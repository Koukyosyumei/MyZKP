#include <cstdio>
#include "./kernels/sumcheck.cu"

__global__ void test_kernel(bool *results) {
    fr_t f_z   = fr_zero();
    fr_t f_o   = fr_one();
    fr_t f_mo  = fr_minus_one();
    fr_t f_2   = {{2, 0, 0, 0}};
    fr_t f_3   = {{3, 0, 0, 0}};
    fr_t f_5   = {{5, 0, 0, 0}};
    fr_t f_7   = {{7, 0, 0, 0}};
    fr_t f_12  = {{12, 0, 0, 0}};
    fr_t f_24  = {{24, 0, 0, 0}};
    fr_t f_29  = {{29, 0, 0, 0}};
    fr_t f_35  = {{35, 0, 0, 0}};
    fr_t f_m2  = {{0x43e1f593efffffffULL, 0x2833e84879b97091ULL, 0xb85045b68181585dULL, 0x30644e72e131a029ULL}};
    fr_t f_m12 = {{0x43e1f593effffff5ULL, 0x2833e84879b97091ULL, 0xb85045b68181585dULL, 0x30644e72e131a029ULL}};

    fr_t add_res = fr_add(f_5, f_7);     // 12 = 5 + 7
    results[0] = fr_eq(add_res, f_12);

    fr_t sub_res = fr_sub(f_7, f_5);     // 2 = 7 - 5
    results[1] = fr_eq(sub_res, f_2);

    fr_t sub_neg = fr_sub(f_5, f_7);     // -2 = 5 - 7
    results[2] = fr_eq(sub_neg, f_m2);  // should wrap around mod p

    fr_t mul_res = fr_mul(f_5, f_7);     // 35 = 5 * 7
    results[3] = fr_eq(mul_res, f_35);

    fr_t f_o_add_mo = fr_add(f_o, f_mo); // 0 = 1 + -1
    results[4] = fr_eq(f_o_add_mo, f_z);

    fr_t f_z_s_12 = fr_sub(f_z, f_12);   // -12 = 0 - 12
    results[5] = fr_eq(f_z_s_12, f_m12);

    fr_t f_m2_mul_m12 = fr_mul(f_m2, f_m12); // 24 = -2 * -12
    results[6] = fr_eq(f_m2_mul_m12, f_24);

    fr_t f_m12_mul_2_add_24 = fr_mul(f_m12, f_2);
    f_m12_mul_2_add_24 = fr_add(f_m12_mul_2_add_24, f_24);
    results[7] = fr_eq(f_m12_mul_2_add_24, f_z);

    // Polynomial: x^2 + 2xy + y^3
    const unsigned int num_sub_mpolys = 3;
    fr_t coeffs[3] = {f_o, f_2, f_o};       // 1*x^2, 2*xy, 1*y^3
    unsigned int expo_flat[5] = {2, 1, 1, 0, 3}; // x^2, x^1 y^1, y^3
    unsigned int offsets[3] = {0, 1, 3};
    unsigned int lens[3] = {1, 2, 2};
    fr_t point[2] = {f_3, f_2}; // x=3, y=2
    unsigned int point_len = 2;

    fr_t m_result = evaluate_mpolynomial(num_sub_mpolys, coeffs, expo_flat, offsets, lens, point, point_len);
    results[8] = fr_eq(f_29, m_result);
}

int main() {
    bool host_results[9];
    bool *dev_results;
    cudaMalloc(&dev_results, sizeof(bool) * 9);

    test_kernel<<<1, 1>>>(dev_results);
    cudaMemcpy(host_results, dev_results, sizeof(bool) * 9, cudaMemcpyDeviceToHost);
    cudaFree(dev_results);

    const char *names[9] = {
        "fr_add", "fr_sub", "fr_sub_wrap", "fr_mul", "fr_add_wrap", "fr_sub_wrap", "f_mul_wrap", "f_mul_wrap1", "f_mpoly"
    };

    int pass_count = 0;
    for (int i = 0; i < 9; i++) {
        printf("%-15s : %s\n", names[i], host_results[i] ? "PASS" : "FAIL");
        if (host_results[i]) pass_count++;
    }

    printf("\n%d / 9 tests passed.\n", pass_count);
    return 0;
}

