#include <cstdint>
#include "field.hpp"

extern "C" __global__ void eval_all_binary_combinations(
    fr_t* results,
    unsigned int offset,
    unsigned int el,
    unsigned int num_sub_mpolys,
    const fr_t* coeffs,
    const unsigned int* expo_flat,
    const unsigned int* offsets,
    const unsigned int* lens
) {
    unsigned int tid = blockIdx.x * blockDim.x + threadIdx.x;
    unsigned int num_combinations = 1 << el; // 2^el

    unsigned int stride = gridDim.x * blockDim.x;
    for (unsigned int idx = tid; idx < num_combinations; idx += stride) {
        fr_t point[32];
        for (unsigned int i = 0; i < el; ++i) {
            point[i] = (idx & (1u << (el - 1 - i))) ? fr_one() : fr_zero();
        }
        fr_t val = evaluate_mpolynomial(num_sub_mpolys, coeffs, expo_flat, offsets, lens, point, el);
        results[offset + idx] = val;
    }
}

/**
 * @brief CUDA kernel to fold multiple polynomial factors into one by computing their pointwise product.
 *
 * @details This kernel calculates the product g(x) = p_1(x) * p_2(x) * ... * p_d(x) over the
 * entire evaluation domain. Each thread computes the product for one point in the domain.
 * The result for each point `idx` is written back in-place to `buf[idx]`.
 *
 * Example layout for evals with num_factors = 3 and domain_size = 8:
 * evals = [p_1(000), p_1(001), p_1(010), p_1(011), p_1(100), p_1(101), p_1(110), p_1(111),
 *          p_2(000), p_2(001), p_2(010), p_2(011), p_2(100), p_2(101), p_2(110), p_2(111),
 *          p_3(000), p_3(001), p_3(010), p_2(011), p_3(100), p_3(101), p_3(110), p_3(111)]
 *
 * @param buf Buffer containing the evaluations of the polynomial factors. The results are written to the start of this buffer.
 * @param domain_size The size of the evaluation domain.
 * @param num_factors The number of polynomial factors to fold.
 */
extern "C" __global__ void fold_factors_pointwise(fr_t* buf, unsigned int domain_size, unsigned int num_factors) {
    const int tid = threadIdx.x;
    int idx = blockIdx.x * blockDim.x + tid;
    while (idx < domain_size) {
	fr_t result = buf[idx];
	for (int i = 1; i < num_factors; i++) {
	    result = fr_mul(result, buf[idx + i * domain_size]);
	}
	buf[idx] = result;
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
 * @param evals The buffer of polynomial evaluations.
 * @param domain_size The total domain size for a single factor polynomial (e.g., 2^l).
 * @param num_remaining_vars The number of variables remaining in the polynomial.
 * @param challenge A pointer to the random challenge `r` for the current round.
 * @param num_blocks_per_factor The number of thread blocks assigned to process a single polynomial factor.
 */
extern "C" __global__ void fold_into_half(
    fr_t* evals, unsigned int domain_size, unsigned int num_remaining_vars, const fr_t* challenge, unsigned int num_blocks_per_factor) {
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
 * @param result Output buffer to store the evaluated polynomial.
 * @param evals The buffer of polynomial evaluations.
 * @param eval_point The point at which to evaluate the first variable of the polynomial.
 * @param num_vars The number of variables in the current polynomial representation.
 * @param domain_size The domain size of the original, unfolded polynomial.
 * @param num_blocks_per_poly The number of thread blocks assigned to process a single polynomial.
 */
extern "C" __global__ void eval_folded_poly(
    fr_t* result, fr_t* evals, const fr_t* eval_point, unsigned int num_vars, unsigned int domain_size, unsigned int num_blocks_per_poly
) {
    const int tid = (blockIdx.x % num_blocks_per_poly) * blockDim.x + threadIdx.x;
    const int stride = 1 << (num_vars - 1);
    const int buf_offset = (blockIdx.x / num_blocks_per_poly) * stride;
    const int poly_offset = (blockIdx.x / num_blocks_per_poly) * domain_size;

    __shared__ bool is_zero, is_one;
    if (threadIdx.x == 0) {
        is_zero = fr_eq(*eval_point, fr_zero());
        is_one = fr_eq(*eval_point, fr_one());
    }
    __syncthreads();

    int idx = tid;
    while (idx < stride) {
        fr_t a = evals[poly_offset + idx];
        fr_t b = evals[poly_offset + idx + stride];
        fr_t r;

        if (is_zero) {r = a;}
        else if (is_one) {r = b;}
        else {
	  fr_t diff = fr_sub(b, a);
	  r = fr_add(a, fr_mul(*eval_point, diff));
	}

	result[buf_offset + idx] = r;
        idx += blockDim.x * num_blocks_per_poly;
    }
}

extern "C" __global__ void sum(fr_t* result, fr_t* data, unsigned int stride, unsigned int index) {
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
