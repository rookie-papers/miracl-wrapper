/**
 * @file test_tools.cpp
 * @brief Functional verification for miracl-wrapper Tools library.
 */

#include "../include/Tools.h"
#include <iostream>
#include <cassert>
#include <string>

using namespace std;

#define TEST_PASS(msg) cout << "[PASS] " << msg << endl;
#define TEST_FAIL(msg) { cout << "[FAIL] " << msg << endl; exit(1); }

gmp_randstate_t state_gmp;
csprng rng_tools;

// ==================================================================
// 1. GMP Convenience Test
// ==================================================================
void Test_GMP_Convenience() {
    cout << "\n--- Test 1: GMP vs MIRACL Syntax & Correctness ---" << endl;

    initState(state_gmp);
    initRNG(&rng_tools);

    // Prepare random data
    BIG a_big, b_big, order;
    BIG_rcopy(order, CURVE_Order);
    randBig(a_big, rng_tools);
    randBig(b_big, rng_tools);

    // Convert to mpz for comparison
    mpz_class a_mpz = BIG_to_mpz(a_big);
    mpz_class b_mpz = BIG_to_mpz(b_big);
    mpz_class q_mpz = BIG_to_mpz(order);

    // A. Modular Addition
    {
        // GMP: Concise syntax
        mpz_class res_mpz = (a_mpz + b_mpz) % q_mpz;

        // MIRACL: Verbose function call
        BIG res_big;
        BIG_modadd(res_big, a_big, b_big, order);

        mpz_class res_converted = BIG_to_mpz(res_big);
        if (res_mpz == res_converted) {
            TEST_PASS("Mod Addition (a + b) % q");
        } else {
            TEST_FAIL("Mod Addition mismatch");
        }
    }

    // B. Modular Subtraction
    {
        // GMP: Handle negative modulo result
        mpz_class res_mpz = (a_mpz - b_mpz) % q_mpz;
        if (res_mpz < 0) res_mpz += q_mpz;

        // MIRACL: Must use modneg + modadd (a - b = a + (-b))
        BIG res_big, neg_b;
        BIG_modneg(neg_b, b_big, order);
        BIG_modadd(res_big, a_big, neg_b, order);

        mpz_class res_converted = BIG_to_mpz(res_big);
        if (res_mpz == res_converted) {
            TEST_PASS("Mod Subtraction (a - b) % q");
        } else {
            TEST_FAIL("Mod Subtraction mismatch");
        }
    }

    // C. Modular Multiplication
    {
        // GMP: Concise syntax
        mpz_class res_mpz = (a_mpz * b_mpz) % q_mpz;

        // MIRACL: Handles internal double-precision reduction
        BIG res_big;
        BIG_modmul(res_big, a_big, b_big, order);

        mpz_class res_converted = BIG_to_mpz(res_big);
        if (res_mpz == res_converted) {
            TEST_PASS("Mod Multiplication (a * b) % q");
        } else {
            TEST_FAIL("Mod Multiplication mismatch");
        }
    }
}

// ==================================================================
// 2. Conversion & ECP Adapter Test
// ==================================================================
void Test_Conversion_And_ECP_Adapter() {
    cout << "\n--- Test 2: Conversions & ECP Adapters ---" << endl;

    initRNG(&rng_tools);

    // A. Round-Trip Conversion (BIG -> mpz -> BIG)
    BIG big_original, big_recovered;
    randBig(big_original, rng_tools);

    mpz_class mpz_val = BIG_to_mpz(big_original);
    mpz_to_BIG(mpz_val, big_recovered);

    if (BIG_comp(big_original, big_recovered) == 0) {
        TEST_PASS("Type Conversion (BIG <-> mpz)");
    } else {
        TEST_FAIL("Type Conversion mismatch");
    }

    // B. ECP Adapter (G1 scalar mul using mpz)
    ECP P_raw, P_wrapper;
    ECP_generator(&P_raw);
    ECP_copy(&P_wrapper, &P_raw);

    BIG r_big;
    randBig(r_big, rng_tools);
    mpz_class r_mpz = BIG_to_mpz(r_big);

    // Compare native MIRACL vs Wrapper
    ECP_mul(&P_raw, r_big);      // Native
    ECP_mul(P_wrapper, r_mpz);   // Wrapper

    if (ECP_equals(&P_raw, &P_wrapper)) {
        TEST_PASS("ECP_mul adapter matches native");
    } else {
        TEST_FAIL("ECP_mul adapter logic failed");
    }

    // C. ECP2 Adapter (G2 scalar mul using mpz)
    ECP2 P2_raw, P2_wrapper;
    ECP2_generator(&P2_raw);
    ECP2_copy(&P2_wrapper, &P2_raw);

    ECP2_mul(&P2_raw, r_big);    // Native
    ECP2_mul(P2_wrapper, r_mpz); // Wrapper

    if (ECP2_equals(&P2_raw, &P2_wrapper)) {
        TEST_PASS("ECP2_mul adapter matches native");
    } else {
        TEST_FAIL("ECP2_mul adapter logic failed");
    }
}

// ==================================================================
// 3. FP12 Adapter Test
// ==================================================================
void Test_FP12_Adapter() {
    cout << "\n--- Test 3: FP12 Functions ---" << endl;

    initRNG(&rng_tools);

    // Setup Pairing inputs
    ECP g1; ECP_generator(&g1);
    ECP2 g2; ECP2_generator(&g2);

    // A. Pairing Adapter
    FP12 gt = e(g1, g2);

    if (!FP12_isunity(&gt)) {
        TEST_PASS("Pairing calculation e(g1, g2)");
    } else {
        TEST_FAIL("Pairing failed (result is unity)");
    }

    // B. Power Adapter (using mpz exponent)
    BIG exp_big;
    randBig(exp_big, rng_tools);
    mpz_class exp_mpz = BIG_to_mpz(exp_big);

    FP12 res_raw, res_wrapper;
    FP12_copy(&res_raw, &gt);
    FP12_copy(&res_wrapper, &gt);

    // Compare native MIRACL vs Wrapper
    FP12_pow(&res_raw, &res_raw, exp_big); // Native
    FP12_pow(res_wrapper, exp_mpz);        // Wrapper (in-place)

    if (FP12_equals(&res_raw, &res_wrapper)) {
        TEST_PASS("FP12_pow adapter matches native");
    } else {
        TEST_FAIL("FP12_pow adapter logic failed");
    }
}

int main() {
    cout << "=== Running Wrapper Verification ===" << endl;

    Test_GMP_Convenience();
    Test_Conversion_And_ECP_Adapter();
    Test_FP12_Adapter();

    cout << "\n=== All Tests Passed ===" << endl;
    return 0;
}