#include <cstdint>
#include "field.hpp"


/**
 * @brief Computes the pointwise product of multiple polynomials at a single evaluation point.
 *
 * @details This function is executed on the GPU device. It calculates the product
 * g(x) = p_1(x) * p_2(x) * ... * p_d(x) for a specific evaluation point `x_id`.
 * The evaluations of all polynomials (p_1, ..., p_d) are stored in a single contiguous
 * array `evals`. The evaluations for each polynomial form a block of size `domain_size`.
 *
 * Example layout for evals with num_factors = 3 and domain_size = 8:
 * evals = [p_1(000), p_1(001), p_1(010), p_1(011), p_1(100), p_1(101), p_1(110), p_1(111),
 *          p_2(000), p_2(001), p_2(010), p_2(011), p_2(100), p_2(101), p_2(110), p_2(111),
 *          p_3(000), p_3(001), p_3(010), p_2(011), p_3(100), p_3(101), p_3(110), p_3(111)]
 *
 * Calling `evaluate_pointwise_product(evals, 1, 8, 3)` would compute p_1(001) * p_2(001) * p_3(001).
 *
 * @param evals Pointer to the array containing the evaluations of all polynomial factors.
 * @param x_id The index of the evaluation point (from 0 to domain_size - 1).
 * @param domain_size The size of the evaluation domain for a single polynomial factor.
 * @param num_factors The number of polynomial factors to multiply.
 * @return The product of the evaluations as a finite field element (fr_t).
 */
__device__ fr_t evaluate_pointwise_product(fr_t* evals, unsigned int x_id, unsigned int domain_size, unsigned int num_factors) {
    fr_t result = fr_one();
    for (int i = 0; i < num_factors; i++) result = fr_mul(result, evals[x_id + i * domain_size]);
    return result;
}

/**
 * @brief CUDA kernel to fold multiple polynomial factors into one by computing their pointwise product.
 *
 * @details This kernel parallelizes the `evaluate_pointwise_product` function across the
 * entire evaluation domain. Each thread computes the product for one point in the domain.
 * The result for each point `idx` is written back in-place to `buf[idx]`, effectively
 * replacing the evaluations of the first polynomial factor with the evaluations of the product polynomial.
 * It uses a grid-stride loop to ensure that all points are processed, regardless of the grid size.
 *
 * @param buf Buffer containing the evaluations of the polynomial factors. The results are written to the start of this buffer.
 * @param domain_size The size of the evaluation domain.
 * @param num_factors The number of polynomial factors to fold.
 */
extern "C" __global__ void fold_factors_pointwise(fr_t* buf, unsigned int domain_size, unsigned int num_factors) {
    const int tid = threadIdx.x;
    int idx = blockIdx.x * blockDim.x + tid;
    while (idx < domain_size) {
        buf[idx] = evaluate_pointwise_product(buf, idx, domain_size, num_factors);
        idx += blockDim.x * gridDim.x;
    }
}

/**
 * @brief CUDA kernel to perform one round of folding in the sum-check protocol.
 *
 * @details In each round of the sum-check protocol, a variable is bound to a random
 * challenge `r`. This kernel updates the polynomial evaluations accordingly. If a polynomial
 * `p(X_j, ...)` is defined over variables starting with `X_j`, this computes the evaluations
 * for the new polynomial `p'(...) = p(r, ...)` over a domain half the size.
 * The update rule is `p(r, x) = p(0, x) + r * (p(1, x) - p(0, x))`.
 * This operation is performed in-place on the `evals` buffer.
 *
 * @param num_remaining_vars The number of variables remaining in the polynomial.
 * @param domain_size The total domain size for a single factor polynomial (e.g., 2^l).
 * @param num_blocks_per_factor The number of thread blocks assigned to process a single polynomial factor.
 * @param evals The buffer of polynomial evaluations.
 * @param challenge A pointer to the random challenge `r` for the current round.
 */
extern "C" __global__ void fold_into_half(
    unsigned int num_remaining_vars, unsigned int domain_size, unsigned int num_blocks_per_factor, fr_t* evals, const fr_t* challenge
) {
    // Calculate the thread's ID within the set of threads assigned to this factor.
    int tid = (blockIdx.x % num_blocks_per_factor) * blockDim.x + threadIdx.x;
    // evals[idx] corresponds to p(..., 0, x'), and evals[idx + stride] to p(..., 1, x').
    const int stride = 1 << (num_remaining_vars - 1);
    // Determine the polynomial factor this thread block is working on.
    const int poly_idx = blockIdx.x / num_blocks_per_factor;
    // Calculate the base offset for the current polynomial factor's evaluations.
    const int offset = poly_idx * domain_size;
    while (tid < stride) {
        int idx = offset + tid;
	// p'(x') = p(0, x') + challenge * (p(1, x') - p(0, x'))
	fr_t tmp = fr_sub(evals[idx + stride], evals[idx]);
	tmp = fr_mul(*challenge, tmp);
	evals[idx] = fr_add(tmp, evals[idx]);
        
	// Stride to the next element for this thread.
	tid += blockDim.x * num_blocks_per_factor;
    }
}

/**
 * @brief CUDA kernel to evaluate a partially folded polynomial at a specific point.
 *
 * @details This is similar to `fold_into_half`, but instead of using a random challenge,
 * it uses a specific evaluation point.
 *
 * @param num_vars The number of variables in the current polynomial representation.
 * @param initial_poly_size The domain size of the original, unfolded polynomial.
 * @param num_blocks_per_poly The number of thread blocks assigned to process a single polynomial.
 * @param evals The buffer of polynomial evaluations.
 * @param result Output buffer to store the evaluated polynomial.
 * @param eval_point The point at which to evaluate the first variable of the polynomial.
 */
extern "C" __global__ void eval_folded_poly(
    unsigned int num_vars, unsigned int initial_poly_size, unsigned int num_blocks_per_poly, fr_t* evals, fr_t* result, const fr_t* eval_point
) {
    int tid = (blockIdx.x % num_blocks_per_poly) * blockDim.x + threadIdx.x;
    const int stride = 1 << (num_vars - 1);
    const int buf_offset = (blockIdx.x / num_blocks_per_poly) * stride;
    const int poly_offset = (blockIdx.x / num_blocks_per_poly) * initial_poly_size;
    while (tid < stride) {
        if (fr_eq(*eval_point, fr_zero())) {result[buf_offset + tid] = evals[poly_offset + tid];}
        else if (fr_eq(*eval_point, fr_one())) {result[buf_offset + tid] = evals[poly_offset + tid + stride];}
        else {
	  result[buf_offset + tid] = fr_sub(evals[poly_offset + tid + stride], evals[poly_offset + tid]);
	  result[buf_offset + tid] = fr_mul(*eval_point, result[buf_offset + tid]);
	  result[buf_offset + tid] = fr_add(result[buf_offset + tid], evals[poly_offset + tid]);
	  // result[buf_offset + tid] = (*eval_point) * (evals[poly_offset + tid + stride] - evals[poly_offset + tid]) + evals[poly_offset + tid];
	}
        tid += blockDim.x * num_blocks_per_poly;
    }
}

extern "C" __global__ void sum(fr_t* data, fr_t* result, unsigned int stride, unsigned int index) {
    const int tid = threadIdx.x;
    for (unsigned int s = stride; s > 0; s >>= 1) {
        int idx = tid;
        while (idx < s) {
            data[idx] = fr_add(data[idx], data[idx + s]);
            idx += blockDim.x;
        }
        __syncthreads();
    }
    if (tid == 0) result[index] = data[0];
}
