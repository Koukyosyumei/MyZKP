#include <cstdint>
#include "field.hpp"


/*
 evals = [p_1(000), p_1(001), p_1(010), p_1(011), p_1(100), p_1(101), p_1(110), p_1(111)
          p_2(000), p_2(001), p_2(010), p_2(011), p_2(100), p_2(101), p_2(110), p_2(111)
          p_3(000), p_3(001), p_3(010), p_2(011), p_3(100), p_3(101), p_3(110), p_3(111)
         ]
 evaluate_pointwise_product(evals, 1, 8, 3) = p_1(001) * p_2(001) * p_3(001)
 */
__device__ fr_t evaluate_pointwise_product(fr_t* evals, unsigned int x_pos, unsigned int domain_size, unsigned int num_factors) {
    fr_t result = fr_one();
    for (int i = 0; i < num_factors; i++) result = fr_mul(result, evals[x_pos + i * domain_size]);
    return result;
}

extern "C" __global__ void fold_factors_pointwise(fr_t* buf, unsigned int domain_size, unsigned int num_factors) {
    const int tid = threadIdx.x;
    int idx = blockIdx.x * blockDim.x + tid;
    while (idx < domain_size) {
        buf[idx] = evaluate_pointwise_product(buf, idx, domain_size, num_factors);
        idx += blockDim.x * gridDim.x;
    }
}

extern "C" __global__ void fold_into_half(
    unsigned int num_remaining_vars, unsigned int domain_size, unsigned int num_blocks_per_factor, fr_t* evals, const fr_t* challenge
) {
    int tid = (blockIdx.x % num_blocks_per_factor) * blockDim.x + threadIdx.x;
    // evals[idx] = p_j(0x), then evals[idx + stride] = p_j(1x)
    const int stride = 1 << (num_remaining_vars - 1);
    // evals[offset : offset + domain_size] = [p_{subpoly_idx}(r_0 r_1, ...), ...., p_{subpoly_idx}(r_0 r_1, ...)]
    const int poly_idx = blockIdx.x / num_blocks_per_factor;
    const int offset = poly_idx * domain_size;
    // blockDim.x * num_blocks_per_poly threads work on folding one factor polynomial p_j.
    while (tid < stride) {
        int idx = offset + tid;
	fr_t tmp = fr_sub(evals[idx + stride], evals[idx]);
	tmp = fr_mul(*challenge, tmp);
	evals[idx] = fr_add(tmp, evals[idx]);
        // evals[idx] = (*challenge) * (evals[idx + stride] - evals[idx]) + evals[idx];
        tid += blockDim.x * num_blocks_per_factor;
    }
}


/*
extern "C" __global__ void fold_into_half(
    unsigned int num_vars, unsigned int initial_poly_size, unsigned int num_blocks_per_poly, fr* evals, fr* result, const fr* eval_point
) {
    int tid = (blockIdx.x % num_blocks_per_poly) * blockDim.x + threadIdx.x;
    const int stride = 1 << (num_vars - 1);
    const int buf_offset = (blockIdx.x / num_blocks_per_poly) * stride;
    const int poly_offset = (blockIdx.x / num_blocks_per_poly) * initial_poly_size;
    while (tid < stride) {
        if (*eval_point == fr::zero()) {result[buf_offset + tid] = evals[poly_offset + tid];}
        else if (*eval_point == fr::one()) {result[buf_offset + tid] = evals[poly_offset + tid + stride];}
        else {
	  result[buf_offset + tid] = fr_sub(evals[poly_offset + tid + stride], evals[poly_offset + tid]);
	  result[buf_offset + tid] = fr_mul(eval_point, result[buf_offset + tid]);
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
}*/
