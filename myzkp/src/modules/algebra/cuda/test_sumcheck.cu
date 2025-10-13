#include <cstdio>
#include <vector>
#include <numeric>
#include "kernels/sumcheck.hpp"

// Macro to check for CUDA errors
#define CUDA_CHECK(err) { \
    cudaError_t e = err; \
    if (e != cudaSuccess) { \
        printf("Cuda error in file '%s' in line %i : %s.\n", __FILE__, __LINE__, cudaGetErrorString(e)); \
        exit(EXIT_FAILURE); \
    } \
}

// Helper to print test results in a consistent format
void check(const char* test_name, bool success) {
    printf("%-30s: %s\n", test_name, success ? "PASS" : "FAIL");
}

// Helper to compare two host arrays of fr_t elements
bool compare_arrays(const std::vector<fr_t>& result, const std::vector<fr_t>& expected) {
    if (result.size() != expected.size()) {
        return false;
    }
    for (size_t i = 0; i < result.size(); ++i) {
        if (!fr_eq(result[i], expected[i])) {
            printf("Mismatch at index %zu: got %llu, expected %llu\n", i, (unsigned long long)result[i].limbs[0], (unsigned long long)expected[i].limbs[0]);
            return false;
        }
    }
    return true;
}

void test_fold_factors_pointwise() {
    const unsigned int domain_size = 256;
    const unsigned int num_factors = 3;
    const size_t total_size = domain_size * num_factors;

    std::vector<fr_t> h_evals(total_size);
    std::vector<fr_t> h_expected(domain_size);

    // Initialize host data and calculate expected result
    for (unsigned int i = 0; i < domain_size; ++i) {
        h_evals[i] = to_fr(i + 1); // p1(i) = i + 1
        h_evals[i + domain_size] = to_fr(i + 2); // p2(i) = i + 2
        h_evals[i + 2 * domain_size] = to_fr(i + 3); // p3(i) = i + 3
        h_expected[i] = to_fr((i + 1) * (i + 2) * (i + 3));
    }

    fr_t* d_evals;
    CUDA_CHECK(cudaMalloc(&d_evals, total_size * sizeof(fr_t)));
    CUDA_CHECK(cudaMemcpy(d_evals, h_evals.data(), total_size * sizeof(fr_t), cudaMemcpyHostToDevice));

    const int block_size = 256;
    const int grid_size = (domain_size + block_size - 1) / block_size;
    fold_factors_pointwise<<<grid_size, block_size>>>(d_evals, domain_size, num_factors);
    CUDA_CHECK(cudaDeviceSynchronize());

    std::vector<fr_t> h_result(domain_size);
    CUDA_CHECK(cudaMemcpy(h_result.data(), d_evals, domain_size * sizeof(fr_t), cudaMemcpyDeviceToHost));
    
    check("fold_factors_pointwise", compare_arrays(h_result, h_expected));

    CUDA_CHECK(cudaFree(d_evals));
}

void test_fold_into_half() {
    const unsigned int num_remaining_vars = 8;
    const unsigned int domain_size = 1 << num_remaining_vars; // 256
    const unsigned int stride = domain_size / 2; // 128

    std::vector<fr_t> h_evals(domain_size);
    std::vector<fr_t> h_expected(stride);
    fr_t challenge = to_fr(3);

    for (unsigned int i = 0; i < domain_size; ++i) {
        h_evals[i] = to_fr(i);
    }
    for (unsigned int i = 0; i < stride; ++i) {
        // p'(i) = p(0,i) + r * (p(1,i) - p(0,i))
        // h_evals[i] + 3 * (h_evals[i+stride] - h_evals[i])
        // i + 3 * ((i+128) - i) = i + 3 * 128 = i + 384
        h_expected[i] = to_fr(i + 384);
    }

    fr_t* d_evals;
    fr_t* d_challenge;
    CUDA_CHECK(cudaMalloc(&d_evals, domain_size * sizeof(fr_t)));
    CUDA_CHECK(cudaMalloc(&d_challenge, sizeof(fr_t)));
    CUDA_CHECK(cudaMemcpy(d_evals, h_evals.data(), domain_size * sizeof(fr_t), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_challenge, &challenge, sizeof(fr_t), cudaMemcpyHostToDevice));

    const int block_size = 256;
    const int grid_size = (stride + block_size - 1) / block_size;
    fold_into_half<<<grid_size, block_size>>>(num_remaining_vars, domain_size, grid_size, d_evals, d_challenge);
    CUDA_CHECK(cudaDeviceSynchronize());

    std::vector<fr_t> h_result(stride);
    CUDA_CHECK(cudaMemcpy(h_result.data(), d_evals, stride * sizeof(fr_t), cudaMemcpyDeviceToHost));

    check("fold_into_half", compare_arrays(h_result, h_expected));

    CUDA_CHECK(cudaFree(d_evals));
    CUDA_CHECK(cudaFree(d_challenge));
}

void test_eval_folded_poly() {
    const unsigned int num_vars = 8;
    const unsigned int initial_poly_size = 1 << num_vars; // 256
    const unsigned int stride = initial_poly_size / 2; // 128

    std::vector<fr_t> h_evals(initial_poly_size);
    for (unsigned int i = 0; i < initial_poly_size; ++i) {
        h_evals[i] = to_fr(i);
    }

    fr_t* d_evals;
    fr_t* d_result;
    fr_t* d_eval_point;
    CUDA_CHECK(cudaMalloc(&d_evals, initial_poly_size * sizeof(fr_t)));
    CUDA_CHECK(cudaMalloc(&d_result, stride * sizeof(fr_t)));
    CUDA_CHECK(cudaMalloc(&d_eval_point, sizeof(fr_t)));
    CUDA_CHECK(cudaMemcpy(d_evals, h_evals.data(), initial_poly_size * sizeof(fr_t), cudaMemcpyHostToDevice));

    const int block_size = 256;
    const int grid_size = (stride + block_size - 1) / block_size;

    // Case 1: eval_point = 0
    fr_t eval_point_0 = fr_zero();
    std::vector<fr_t> h_expected_0(stride);
    for(unsigned int i=0; i<stride; ++i) h_expected_0[i] = to_fr(i);
    CUDA_CHECK(cudaMemcpy(d_eval_point, &eval_point_0, sizeof(fr_t), cudaMemcpyHostToDevice));
    eval_folded_poly<<<grid_size, block_size>>>(num_vars, initial_poly_size, grid_size, d_evals, d_result, d_eval_point);
    CUDA_CHECK(cudaDeviceSynchronize());
    std::vector<fr_t> h_result_0(stride);
    CUDA_CHECK(cudaMemcpy(h_result_0.data(), d_result, stride * sizeof(fr_t), cudaMemcpyDeviceToHost));
    check("eval_folded_poly (point=0)", compare_arrays(h_result_0, h_expected_0));

    // Case 2: eval_point = 1
    fr_t eval_point_1 = fr_one();
    std::vector<fr_t> h_expected_1(stride);
    for(unsigned int i=0; i<stride; ++i) h_expected_1[i] = to_fr(i + stride);
    CUDA_CHECK(cudaMemcpy(d_eval_point, &eval_point_1, sizeof(fr_t), cudaMemcpyHostToDevice));
    eval_folded_poly<<<grid_size, block_size>>>(num_vars, initial_poly_size, grid_size, d_evals, d_result, d_eval_point);
    CUDA_CHECK(cudaDeviceSynchronize());
    std::vector<fr_t> h_result_1(stride);
    CUDA_CHECK(cudaMemcpy(h_result_1.data(), d_result, stride * sizeof(fr_t), cudaMemcpyDeviceToHost));
    check("eval_folded_poly (point=1)", compare_arrays(h_result_1, h_expected_1));

    // Case 3: eval_point = 5
    fr_t eval_point_5 = to_fr(5);
    std::vector<fr_t> h_expected_5(stride);
    for(unsigned int i=0; i<stride; ++i) h_expected_5[i] = to_fr(i + 5 * stride); // i + 5 * 128
    CUDA_CHECK(cudaMemcpy(d_eval_point, &eval_point_5, sizeof(fr_t), cudaMemcpyHostToDevice));
    eval_folded_poly<<<grid_size, block_size>>>(num_vars, initial_poly_size, grid_size, d_evals, d_result, d_eval_point);
    CUDA_CHECK(cudaDeviceSynchronize());
    std::vector<fr_t> h_result_5(stride);
    CUDA_CHECK(cudaMemcpy(h_result_5.data(), d_result, stride * sizeof(fr_t), cudaMemcpyDeviceToHost));
    check("eval_folded_poly (point=5)", compare_arrays(h_result_5, h_expected_5));

    CUDA_CHECK(cudaFree(d_evals));
    CUDA_CHECK(cudaFree(d_result));
    CUDA_CHECK(cudaFree(d_eval_point));
}

void test_sum() {
    const unsigned int stride = 256;
    std::vector<fr_t> h_data(stride);
    uint64_t expected_sum = 0;
    for (unsigned int i = 0; i < stride; ++i) {
        h_data[i] = to_fr(i + 1);
        expected_sum += (i + 1);
    }

    fr_t* d_data;
    fr_t* d_result;
    CUDA_CHECK(cudaMalloc(&d_data, stride * sizeof(fr_t)));
    CUDA_CHECK(cudaMalloc(&d_result, sizeof(fr_t)));
    CUDA_CHECK(cudaMemcpy(d_data, h_data.data(), stride * sizeof(fr_t), cudaMemcpyHostToDevice));

    sum<<<1, stride>>>(d_data, d_result, stride, 0);
    CUDA_CHECK(cudaDeviceSynchronize());

    fr_t h_result;
    CUDA_CHECK(cudaMemcpy(&h_result, d_result, sizeof(fr_t), cudaMemcpyDeviceToHost));

    check("sum reduction", h_result.limbs[0] == expected_sum);

    CUDA_CHECK(cudaFree(d_data));
    CUDA_CHECK(cudaFree(d_result));
}

int main() {
    printf("Running Sum-Check Kernel Unit Tests...\n\n");
    test_fold_factors_pointwise();
    test_fold_into_half();
    test_eval_folded_poly();
    test_sum();
    printf("\n...Tests complete.\n");
    return 0;
}
