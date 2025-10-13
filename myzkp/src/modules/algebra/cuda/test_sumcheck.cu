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

int main() {
    printf("Running Sum-Check Kernel Unit Tests...\n\n");
    test_fold_factors_pointwise();
    // test_fold_into_half();
    // test_eval_folded_poly();
    // test_sum();
    printf("\n...Tests complete.\n");
    return 0;
}
