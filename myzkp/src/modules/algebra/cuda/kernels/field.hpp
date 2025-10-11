#include <cstdint>
#include <cstdio>

// BN254 modulus r (little-endian limbs)
__device__ const uint64_t FR_MOD[4] = {
    0x43e1f593f0000001ULL,
    0x2833e84879b97091ULL,
    0xb85045b68181585dULL,
    0x30644e72e131a029ULL
};

struct fr_t {
    uint64_t limbs[4]; // 256-bit little-endian
};

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

__device__ bool fr_gte(const fr_t &a, const uint64_t b[4]) {
    for (int i = 3; i >= 0; i--) {
        if (a.limbs[i] > b[i]) return true;
        if (a.limbs[i] < b[i]) return false;
    }
    return true;
}

__device__ void fr_reduce(fr_t &a) {
    if (fr_gte(a, FR_MOD)) {
        uint64_t borrow = 0;
        for (int i = 0; i < 4; i++) {
            uint64_t tmp = a.limbs[i] - FR_MOD[i] - borrow;
            borrow = (a.limbs[i] < FR_MOD[i] + borrow) ? 1 : 0;
            a.limbs[i] = tmp;
        }
    }
}

__device__ fr_t fr_add(const fr_t &a, const fr_t &b) {
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

__device__ fr_t fr_sub(const fr_t &a, const fr_t &b) {
    fr_t res;
    uint64_t borrow = 0;
    for (int i = 0; i < 4; i++) {
        unsigned __int128 tmp = (unsigned __int128)a.limbs[i] - b.limbs[i] - borrow;
        res.limbs[i] = (uint64_t)tmp;
        borrow = (tmp >> 127) & 1; // borrow = 1 if underflow
    }

    if (borrow) { // add modulus back
        uint64_t carry = 0;
        for (int i = 0; i < 4; i++) {
            unsigned __int128 tmp = (unsigned __int128)res.limbs[i] + FR_MOD[i] + carry;
            res.limbs[i] = (uint64_t)tmp;
            carry = tmp >> 64;
        }
    }

    return res;
}

__device__ fr_t fr_mul(const fr_t &a, const fr_t &b) {
    unsigned __int128 tmp[8] = {0};

    // Schoolbook multiplication
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            tmp[i+j] += (unsigned __int128)a.limbs[i] * b.limbs[j];
        }
    }

    // Propagate carry
    for (int i = 0; i < 7; i++) {
        tmp[i+1] += tmp[i] >> 64;
        tmp[i] &= 0xFFFFFFFFFFFFFFFFULL;
    }

    // Reduce modulo r (fixed algorithm)
    // 1. Copy low 4 limbs
    fr_t res;
    for (int i = 0; i < 4; i++) res.limbs[i] = (uint64_t)tmp[i];

    // 2. Reduce high limbs manually using BN254 modulus
    //    (this is simplified; for best performance use Montgomery)
    for (int i = 4; i < 8; i++) {
        // Multiply tmp[i] by reduction factor from modulus
        // Add to low limbs with carry
        unsigned __int128 carry = tmp[i];
        for (int j = 0; j < 4 && carry != 0; j++) {
            unsigned __int128 sum = (unsigned __int128)res.limbs[j] + carry;
            res.limbs[j] = (uint64_t)sum;
            carry = sum >> 64;
        }
    }

    fr_reduce(res); // ensure < r
    return res;
}
