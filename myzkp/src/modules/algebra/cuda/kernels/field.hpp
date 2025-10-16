#include <cstdint>
#include <cstdio>

struct fr_t {
    uint64_t limbs[4]; // 256-bit little-endian
};
	
// Device-only constants (for device compilation)
__device__ __constant__ fr_t FR_MOD_DEV = {{
    0x43e1f593f0000001ULL,
    0x2833e84879b97091ULL,
    0xb85045b68181585dULL,
    0x30644e72e131a029ULL
}};

// R = 2^256
// Montgomery constant: FR_MOD * FR_MOD_INV = -1 (mod R)

__device__ __constant__ fr_t FR_MOD_INV_DEV = {{
    0xc2e1f593efffffffULL,
    0x6586864b4c6911b3ULL,
    0xe39a982899062391ULL,
    0x73f82f1d0d8341b2ULL
}};

__device__ __constant__ fr_t R2_MOD_N_DEV = {{
    0x1bb8e645ae216da7ULL,
    0x53fe3ab1e35c59e3ULL,
    0x8c49833d53bb8085ULL,
    0x0216d0b17f4e44a5ULL
}};

// Host constexpr copies (for host compilation)
static constexpr fr_t FR_MOD_HOST = {{
    0x43e1f593f0000001ULL,
    0x2833e84879b97091ULL,
    0xb85045b68181585dULL,
    0x30644e72e131a029ULL
}};

static constexpr fr_t FR_MOD_INV_HOST = {{
    0xc2e1f593efffffffULL,
    0x6586864b4c6911b3ULL,
    0xe39a982899062391ULL,
    0x73f82f1d0d8341b2ULL
}};

static constexpr fr_t R2_MOD_N_HOST = {{
    0x1bb8e645ae216da7ULL,
    0x53fe3ab1e35c59e3ULL,
    0x8c49833d53bb8085ULL,
    0x0216d0b17f4e44a5ULL
}};

inline __host__ __device__ const fr_t* get_FR_MOD_ptr() {
#ifdef __CUDA_ARCH__
    return &FR_MOD_DEV;
#else
    return &FR_MOD_HOST;
#endif
}

inline __host__ __device__ const fr_t* get_FR_MOD_INV_ptr() {
#ifdef __CUDA_ARCH__
    return &FR_MOD_INV_DEV;
#else
    return &FR_MOD_INV_HOST;
#endif
}

inline __host__ __device__ const fr_t* get_R2_MOD_N_ptr() {
#ifdef __CUDA_ARCH__
    return &R2_MOD_N_DEV;
#else
    return &R2_MOD_N_HOST;
#endif
}

inline fr_t to_fr(uint64_t val) {
    return {{val, 0, 0, 0}};
}

__host__ __device__ bool fr_eq(const fr_t &a, const fr_t &b) {
    for (int i = 3; i >= 0; i--) {
       if (a.limbs[i] != b.limbs[i]) return false;
    }
    return true;
}

__host__ __device__ void fr_print(const fr_t &x) {
    printf("0x");
    for (int i = 3; i >= 0; --i) printf("%016llx", (unsigned long long)x.limbs[i]);
}

__host__ __device__ __forceinline__ fr_t fr_zero() {
    fr_t z = {{0, 0, 0, 0}};
    return z;
}

__host__ __device__ __forceinline__ fr_t fr_one() {
    fr_t o = {{1, 0, 0, 0}};
    return o;
}

__host__ __device__ __forceinline__ fr_t fr_minus_one() {
    fr_t m1 = {
        {
            0x43e1f593f0000000ULL,
            0x2833e84879b97091ULL,
            0xb85045b68181585dULL,
            0x30644e72e131a029ULL
        }
    };
    return m1;
}

__host__ __device__ bool fr_gte(const fr_t &a, const fr_t &b) {
    for (int i = 3; i >= 0; i--) {
        if (a.limbs[i] > b.limbs[i]) return true;
        if (a.limbs[i] < b.limbs[i]) return false;
    }
    return true;
}

__host__ __device__ void fr_reduce(fr_t &a) {
    if (fr_gte(a, *get_FR_MOD_ptr())) {
        uint64_t borrow = 0;
	const fr_t &mod = *get_FR_MOD_ptr();
        for (int i = 0; i < 4; i++) {
            uint64_t tmp = a.limbs[i] - mod.limbs[i] - borrow;
            borrow = (a.limbs[i] < mod.limbs[i] + borrow) ? 1 : 0;
            a.limbs[i] = tmp;
        }
    }
}

__host__ __device__ fr_t fr_add(const fr_t &a, const fr_t &b) {
    fr_t res;
    uint64_t carry = 0;
    for (int i = 0; i < 4; i++) {
        unsigned __int128 tmp = (unsigned __int128)a.limbs[i] + b.limbs[i] + carry;
        res.limbs[i] = (uint64_t)tmp;
        carry = tmp >> 64;
    }
    fr_reduce(res);
    return res;
}

__host__ __device__ fr_t fr_sub(const fr_t &a, const fr_t &b) {
    fr_t res;
    uint64_t borrow = 0;
    for (int i = 0; i < 4; i++) {
        unsigned __int128 tmp = (unsigned __int128)a.limbs[i] - b.limbs[i] - borrow;
        res.limbs[i] = (uint64_t)tmp;
        borrow = (tmp >> 127) & 1; // borrow = 1 if underflow
    }

    const fr_t &mod = *get_FR_MOD_ptr();
    if (borrow) { // add modulus back
        uint64_t carry = 0;
        for (int i = 0; i < 4; i++) {
            unsigned __int128 tmp = (unsigned __int128)res.limbs[i] + mod.limbs[i] + carry;
            res.limbs[i] = (uint64_t)tmp;
            carry = tmp >> 64;
        }
    }

    return res;
}

// R = 2^256
// Montgomery constant: FR_MOD * FR_MOD_INV = -1 (mod R)

// --- Helper Types and Functions ---

// Helper struct for a 512-bit integer
struct u512_t {
    uint64_t limbs[8];
};

// --- New Helper Functions ---

/**
 * @brief Adds two 512-bit numbers. res = a + b.
 */
__host__ __device__ void add_512(u512_t &res, const u512_t &a, const u512_t &b) {
    unsigned __int128 carry = 0;
    for (int i = 0; i < 8; i++) {
        carry += (unsigned __int128)a.limbs[i] + b.limbs[i];
        res.limbs[i] = (uint64_t)carry;
        carry >>= 64;
    }
}

/**
 * @brief Subtracts one 256-bit number from another. res = a - b.
 * Assumes a >= b.
 */
__host__ __device__ void sub_256(fr_t &res, const fr_t &a, const fr_t &b) {
    uint64_t borrow = 0;
    for (int i = 0; i < 4; i++) {
        unsigned __int128 tmp = (unsigned __int128)a.limbs[i] - b.limbs[i] - borrow;
        res.limbs[i] = (uint64_t)tmp;
        borrow = (tmp >> 127) & 1; // 1 if underflow, 0 otherwise
    }
}

/**
 * @brief Multiplies two 256-bit numbers to produce a 512-bit result.
 * @param res The 512-bit result (a * b).
 * @param a The first 256-bit operand.
 * @param b The second 256-bit operand.
 */
__host__ __device__ void mul_512(u512_t &res, const fr_t &a, const fr_t &b) {
    // Initialize result limbs to 0
    for(int i = 0; i < 8; ++i) res.limbs[i] = 0;

    for (int i = 0; i < 4; i++) {
        unsigned __int128 carry = 0;
        for (int j = 0; j < 4; j++) {
            // Add product and existing value in limb
            carry += (unsigned __int128)a.limbs[i] * b.limbs[j] + res.limbs[i + j];
            res.limbs[i + j] = (uint64_t)carry;
            carry >>= 64;
        }
        // Propagate final carry
        res.limbs[i + 4] += (uint64_t)carry;
    }
}

// --- Core Montgomery Reduction and Multiplication Functions ---

/**
 * @brief Reduces a 512-bit number T using Montgomery reduction (REDC).
 * @param t The 512-bit number to reduce.
 * @return The result (T * R^-1) mod N, where R = 2^256.
 */
__host__ __device__ fr_t mont_reduce(const u512_t &t) {
    // 1. m = (T mod R) * N' mod R
    // T mod R is just the lower 4 limbs of T.
    // N' is FR_MOD_INV.
    fr_t t_low;
    for (int i = 0; i < 4; ++i) t_low.limbs[i] = t.limbs[i];
    
    u512_t m_full;
    mul_512(m_full, t_low, *get_FR_MOD_INV_ptr()); // m = (T mod R) * N'
    // We only need the lower 256 bits of m_full for the next step.

    // 2. tmp = (T + m*N)
    u512_t mN;
    mul_512(mN, *(const fr_t*)m_full.limbs, *get_FR_MOD_ptr());

    u512_t sum;
    add_512(sum, t, mN);
    
    // 3. result = tmp / R
    // By construction, the lower 256 bits of `sum` are zero,
    // so division by R (2^256) is just a right shift by 256 bits.
    fr_t res;
    for (int i = 0; i < 4; ++i) res.limbs[i] = sum.limbs[i + 4];

    // 4. Final conditional subtraction
    if (fr_gte(res, *get_FR_MOD_ptr())) {
        sub_256(res, res, *get_FR_MOD_ptr());
    }
    
    return res;
}

/**
 * @brief Multiplies two numbers that are already in the Montgomery domain.
 * @param a_mont A number in Montgomery form (a * R mod N).
 * @param b_mont Another number in Montgomery form (b * R mod N).
 * @return The result (a * b * R mod N).
 */
__host__ __device__ fr_t fr_mul_mont(const fr_t &a_mont, const fr_t &b_mont) {
    u512_t prod;
    mul_512(prod, a_mont, b_mont);
    return mont_reduce(prod);
}

/**
 * @brief Converts a number into the Montgomery domain.
 * @param a A standard number.
 * @return The number in Montgomery form (a * R mod N).
 */
__host__ __device__ fr_t to_mont(const fr_t &a) {
    u512_t prod;
    mul_512(prod, a, *get_R2_MOD_N_ptr());
    return mont_reduce(prod);
}

/**
 * @brief Multiplies two standard numbers. Handles all conversions.
 * @param a The first operand.
 * @param b The second operand.
 * @return The result (a * b) mod N.
 */
__host__ __device__ fr_t fr_mul(const fr_t &a, const fr_t &b) {
    // 1. Convert both operands to Montgomery form
    fr_t a_mont = to_mont(a);
    fr_t b_mont = to_mont(b);
    
    // 2. Multiply in the Montgomery domain
    fr_t res_mont = fr_mul_mont(a_mont, b_mont);
    
    // 3. Convert the result back from Montgomery form
    u512_t zero_ext_res = {{0}};
    for(int i = 0; i < 4; ++i) zero_ext_res.limbs[i] = res_mont.limbs[i];
    
    return mont_reduce(zero_ext_res);
}

__host__ __device__ fr_t fr_pow_uint(fr_t base, unsigned int exp) {
    fr_t res = fr_one();
    fr_t cur = base;
    while (exp > 0) {
        if (exp & 1u) res = fr_mul(res, cur);
        cur = fr_mul(cur, cur);
        exp >>= 1;
    }
    return res;
}

__host__ __device__ fr_t evaluate_mpolynomial(
    unsigned int num_sub_mpolys,
    const fr_t* coeffs,
    const unsigned int* expo_flat,
    const unsigned int* offsets,
    const unsigned int* lens,
    const fr_t* point,
    unsigned int point_len
) {
    fr_t acc = fr_zero();
    for (unsigned int i = 0; i < num_sub_mpolys; ++i) {
        fr_t prod = coeffs[i];
        unsigned int off = offsets[i];
        unsigned int len = lens[i];
        for (unsigned int j = 0; j < len; ++j) {
            if (j >= point_len) continue;
            prod = fr_mul(prod, fr_pow_uint(point[j], expo_flat[off + j]));
        }
        acc = fr_add(acc, prod);
    }
    return acc;
}
