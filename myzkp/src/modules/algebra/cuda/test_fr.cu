#include <cstdio>
#include "./kernels/field.hpp"

__global__ void test_kernel(bool *results) {
    fr_t a = {{1, 0, 0, 0}};
    fr_t b = {{2, 0, 0, 0}};

    // fr_add
    fr_t add_res = fr_add(a, b);
    results[0] = (add_res.limbs[0] == 3);

    // fr_sub
    fr_t sub_res = fr_sub(b, a);
    results[1] = (sub_res.limbs[0] == 1);

    // fr_sub wrap-around
    fr_t sub_neg = fr_sub(a, b);
    results[2] = fr_gte(sub_neg, a.limbs); // should wrap around mod p

    // fr_mul
    fr_t mul_res = fr_mul(b, b);
    results[3] = (mul_res.limbs[0] == 4);

    // fr_gtr / fr_gte
    //results[4] = fr_gtr(b, a);
    //results[5] = fr_gte(b, b);
}

int main() {
    bool host_results[6];
    bool *dev_results;
    cudaMalloc(&dev_results, sizeof(bool) * 6);

    test_kernel<<<1, 1>>>(dev_results);
    cudaMemcpy(host_results, dev_results, sizeof(bool) * 6, cudaMemcpyDeviceToHost);
    cudaFree(dev_results);

    const char *names[6] = {
        "fr_add", "fr_sub", "fr_sub_wrap", "fr_mul", "fr_gtr", "fr_gte"
    };

    int pass_count = 0;
    for (int i = 0; i < 4; i++) {
        printf("%-15s : %s\n", names[i], host_results[i] ? "PASS" : "FAIL");
        if (host_results[i]) pass_count++;
    }

    printf("\n%d / 6 tests passed.\n", pass_count);
    return 0;
}

