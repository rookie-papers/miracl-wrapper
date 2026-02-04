#include "../include/Tools.h"
#include "benchmark/benchmark.h"

#include <iostream>
#include <cstring>
#include <cstdlib>
#include <ctime>

using namespace std;

// ==================================================================
// Global Variables & Constants
// ==================================================================

gmp_randstate_t state_BM;
csprng rng;

// BLS12-381 Curve Order
const mpz_class q("0x73EDA753299D7D483339D80809A1D80553BDA402FFFE5BFEFFFFFFFF00000001");

// ==================================================================
// GMP Arithmetic Benchmarks
// ==================================================================

void GMP_add(benchmark::State &state) {
    initState(state_BM);
    mpz_class a = rand_mpz(state_BM);
    mpz_class b = rand_mpz(state_BM);
    for (auto _: state) {
        benchmark::DoNotOptimize(a + b);
    }
}

void GMP_sub(benchmark::State &state) {
    initState(state_BM);
    mpz_class a = rand_mpz(state_BM);
    mpz_class b = rand_mpz(state_BM);
    for (auto _: state) {
        benchmark::DoNotOptimize(a - b);
    }
}

void GMP_mul(benchmark::State &state) {
    initState(state_BM);
    mpz_class a = rand_mpz(state_BM);
    mpz_class b = rand_mpz(state_BM);
    for (auto _: state) {
        benchmark::DoNotOptimize(a * b);
    }
}

void GMP_div(benchmark::State &state) {
    initState(state_BM);
    mpz_class a = rand_mpz(state_BM);
    mpz_class b = rand_mpz(state_BM);
    for (auto _: state) {
        benchmark::DoNotOptimize(a / b);
    }
}

void GMP_modadd(benchmark::State &state) {
    initState(state_BM);
    mpz_class a = rand_mpz(state_BM);
    mpz_class b = rand_mpz(state_BM);
    for (auto _: state) {
        mpz_class res = (a + b) % q;
        benchmark::DoNotOptimize(res);
    }
}

void GMP_modsub(benchmark::State &state) {
    initState(state_BM);
    mpz_class a = rand_mpz(state_BM);
    mpz_class b = rand_mpz(state_BM);
    for (auto _: state) {
        mpz_class res = (a - b) % q;
        benchmark::DoNotOptimize(res);
    }
}

void GMP_modmul(benchmark::State &state) {
    initState(state_BM);
    mpz_class a = rand_mpz(state_BM);
    mpz_class b = rand_mpz(state_BM);
    for (auto _: state) {
        mpz_class res = (a * b) % q;
        benchmark::DoNotOptimize(res);
    }
}

void GMP_moddiv(benchmark::State &state) {
    initState(state_BM);
    mpz_class a = rand_mpz(state_BM);
    mpz_class b = rand_mpz(state_BM);
    for (auto _: state) {
        mpz_class res = (a / b) % q;
        benchmark::DoNotOptimize(res);
    }
}

void GMP_inv(benchmark::State &state) {
    initState(state_BM);
    mpz_class a = rand_mpz(state_BM);
    for (auto _: state) {
        mpz_class res = invert_mpz(a, q);
        benchmark::DoNotOptimize(res);
    }
}

void GMP_pow(benchmark::State &state) {
    initState(state_BM);
    mpz_class a = rand_mpz(state_BM);
    mpz_class b = rand_mpz(state_BM);
    for (auto _: state) {
        mpz_class res = pow_mpz(a, b, q);
        benchmark::DoNotOptimize(res);
    }
}

// ==================================================================
// MIRACL Core BIG Arithmetic Benchmarks
// ==================================================================

void Miracl_add(benchmark::State &state) {
    initRNG(&rng);
    BIG a, b;
    randBig(a, rng);
    randBig(b, rng);
    for (auto _: state) {
        BIG_add(a, a, b);
    }
}

void Miracl_sub(benchmark::State &state) {
    initRNG(&rng);
    BIG a, b;
    randBig(a, rng);
    randBig(b, rng);
    for (auto _: state) {
        BIG_sub(a, a, b);
    }
}

void Miracl_mul(benchmark::State &state) {
    initRNG(&rng);
    BIG a, b;
    DBIG res;
    randBig(a, rng);
    randBig(b, rng);
    for (auto _: state) {
        BIG_mul(res, a, b);
    }
}

void Miracl_modadd(benchmark::State &state) {
    initRNG(&rng);
    BIG a, b, order;
    randBig(a, rng);
    randBig(b, rng);
    BIG_rcopy(order, CURVE_Order);
    for (auto _: state) {
        BIG_modadd(a, a, b, order);
    }
}

void Miracl_modsub(benchmark::State &state) {
    initRNG(&rng);
    BIG a, b, order;
    randBig(a, rng);
    randBig(b, rng);
    BIG_rcopy(order, CURVE_Order);
    for (auto _: state) {
        // Simulation of modsub if native function is missing
        BIG_sub(a, a, b);
        BIG_norm(a);
    }
}

void Miracl_modmul(benchmark::State &state) {
    initRNG(&rng);
    BIG a, b, order;
    randBig(a, rng);
    randBig(b, rng);
    BIG_rcopy(order, CURVE_Order);
    for (auto _: state) {
        BIG_modmul(a, a, b, order);
    }
}

void Miracl_inv(benchmark::State &state) {
    initRNG(&rng);
    BIG a, b, order;
    randBig(a, rng);
    BIG_rcopy(order, CURVE_Order);
    for (auto _: state) {
        BIG_invmodp(a, a, order);
    }
}

// ==================================================================
// ECC (Elliptic Curve Cryptography) Benchmarks
// ==================================================================

void Miracl_ECP_add(benchmark::State &state) {
    initRNG(&rng);
    BIG a, b;
    randBig(a, rng);
    randBig(b, rng);
    ECP P1, P2;
    ECP_generator(&P1);
    ECP_generator(&P2);
    ECP_mul(&P1, a);
    ECP_mul(&P2, b);
    for (auto _: state) {
        ECP_add(&P1, &P2);
    }
}

void Miracl_ECP_mul(benchmark::State &state) {
    initRNG(&rng);
    BIG a;
    randBig(a, rng);
    ECP P1;
    ECP_generator(&P1);
    for (auto _: state) {
        ECP_mul(&P1, a);
    }
}

void Miracl_ECP2_add(benchmark::State &state) {
    initRNG(&rng);
    BIG a, b;
    randBig(a, rng);
    randBig(b, rng);
    ECP2 P1, P2;
    ECP2_generator(&P1);
    ECP2_generator(&P2);
    ECP2_mul(&P1, a);
    ECP2_mul(&P2, b);
    for (auto _: state) {
        ECP2_add(&P1, &P2);
    }
}

void Miracl_ECP2_mul(benchmark::State &state) {
    initRNG(&rng);
    BIG a;
    randBig(a, rng);
    ECP2 P1;
    ECP2_generator(&P1);
    for (auto _: state) {
        ECP2_mul(&P1, a);
    }
}

void Miracl_pair(benchmark::State &state) {
    initRNG(&rng);
    ECP P1;
    ECP2 P2;
    ECP_generator(&P1);
    ECP2_generator(&P2);
    for (auto _: state) {
        FP12 r = e(P1, P2);
        benchmark::DoNotOptimize(r);
    }
}

void Miracl_GT_mul(benchmark::State &state) {
    initRNG(&rng);
    ECP P1;
    ECP2 P2;
    ECP_generator(&P1);
    ECP2_generator(&P2);
    FP12 gt1 = e(P1, P2);
    FP12 gt2; FP12_copy(&gt2, &gt1);

    for (auto _: state) {
        FP12_mul(&gt1, &gt2);
    }
}

void Miracl_GT_pow(benchmark::State &state) {
    initRNG(&rng);
    BIG a;
    randBig(a, rng);
    ECP P1;
    ECP2 P2;
    ECP_generator(&P1);
    ECP2_generator(&P2);
    FP12 gt = e(P1, P2);
    for (auto _: state) {
        FP12_pow(&gt, &gt, a);
    }
}

// ==================================================================
// Utilities (Hash & AES) Benchmarks
// ==================================================================

void Miracl_hash(benchmark::State &state) {
    initState(state_BM);
    mpz_class a = rand_mpz(state_BM);
    BIG order;
    BIG_rcopy(order, CURVE_Order);
    BIG ret;

    for (auto _: state) {
        octet oct_a = mpzToOctet(a);
        hashZp256(ret, &oct_a, order);
        free(oct_a.val);
    }
}

void Miracl_hashToPoint(benchmark::State &state) {
    initState(state_BM);
    mpz_class a = rand_mpz(state_BM);
    BIG order;
    BIG_rcopy(order, CURVE_Order);
    mpz_class q = BIG_to_mpz(order);
    ECP P;

    for (auto _: state) {
        P = hashToPoint(a, q);
    }
}

void Miracl_AES_Enc(benchmark::State &state) {
    int KK = 32; // Key length (256-bit)
    aes a;
    char key[32];
    char block[16];
    char iv[16];

    // Initialize dummy data
    for (int i = 0; i < KK; i++) key[i] = 5;
    for (int i = 0; i < 16; i++) { iv[i] = i; block[i] = i; }

    AES_init(&a, CTR16, KK, key, iv);

    for (auto _: state) {
        AES_encrypt(&a, block);
    }
    AES_end(&a);
}

void Miracl_AES_Dec(benchmark::State &state) {
    int KK = 32;
    aes a;
    char key[32];
    char block[16];
    char iv[16];

    // Initialize dummy data
    for (int i = 0; i < KK; i++) key[i] = 5;
    for (int i = 0; i < 16; i++) { iv[i] = i; block[i] = i; }

    // Setup: Encrypt once to get valid ciphertext for decryption loop
    AES_init(&a, CTR16, KK, key, iv);
    AES_encrypt(&a, block);
    AES_end(&a);

    // Re-init for benchmark loop
    AES_init(&a, CTR16, KK, key, iv);
    for (auto _: state) {
        AES_decrypt(&a, block);
    }
    AES_end(&a);
}

// ==================================================================
// Register Benchmarks
// ==================================================================

// GMP
BENCHMARK(GMP_add);
BENCHMARK(GMP_sub);
BENCHMARK(GMP_mul);
BENCHMARK(GMP_div);
BENCHMARK(GMP_modadd);
BENCHMARK(GMP_modsub);
BENCHMARK(GMP_modmul);
BENCHMARK(GMP_moddiv);
BENCHMARK(GMP_inv);
BENCHMARK(GMP_pow);

// MIRACL Core
BENCHMARK(Miracl_add);
BENCHMARK(Miracl_sub);
BENCHMARK(Miracl_mul);
BENCHMARK(Miracl_modadd);
BENCHMARK(Miracl_modsub);
BENCHMARK(Miracl_modmul);
BENCHMARK(Miracl_inv);

// ECC
BENCHMARK(Miracl_ECP_add);
BENCHMARK(Miracl_ECP_mul);
BENCHMARK(Miracl_ECP2_add);
BENCHMARK(Miracl_ECP2_mul);
BENCHMARK(Miracl_pair);
BENCHMARK(Miracl_GT_mul);
BENCHMARK(Miracl_GT_pow);

// Utils
BENCHMARK(Miracl_hash);
BENCHMARK(Miracl_hashToPoint);
BENCHMARK(Miracl_AES_Enc);
BENCHMARK(Miracl_AES_Dec);

BENCHMARK_MAIN();