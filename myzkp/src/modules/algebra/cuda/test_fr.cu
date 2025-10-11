#include <cstdio>
#include "./kernels/field.hpp"

__global__ void test_kernel(bool *results) {
    fr_t f_z   = fr_zero();
    fr_t f_o   = fr_one();
    fr_t f_mo  = fr_minus_one();
    fr_t f_2   = {{2, 0, 0, 0}};
    fr_t f_5   = {{5, 0, 0, 0}};
    fr_t f_7   = {{7, 0, 0, 0}};
    fr_t f_12  = {{12, 0, 0, 0}};
    fr_t f_24  = {{24, 0, 0, 0}};
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

    fr_t f_m2_mul_m12 = fr_mul1(f_m2, f_m12); // 24 = -2 * -12
    fr_print(f_m2_mul_m12);
    results[6] = fr_eq(f_m2_mul_m12, f_24);
}

int main() {
    bool host_results[7];
    bool *dev_results;
    cudaMalloc(&dev_results, sizeof(bool) * 7);

    test_kernel<<<1, 1>>>(dev_results);
    cudaMemcpy(host_results, dev_results, sizeof(bool) * 7, cudaMemcpyDeviceToHost);
    cudaFree(dev_results);

    const char *names[7] = {
        "fr_add", "fr_sub", "fr_sub_wrap", "fr_mul", "fr_add_wrap", "fr_sub_wrap", "f_mul_wrap"
    };

    int pass_count = 0;
    for (int i = 0; i < 7; i++) {
        printf("%-15s : %s\n", names[i], host_results[i] ? "PASS" : "FAIL");
        if (host_results[i]) pass_count++;
    }

    printf("\n%d / 7 tests passed.\n", pass_count);
    return 0;
}

