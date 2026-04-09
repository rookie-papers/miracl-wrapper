/**
 * @file test_tools.cpp
 * @brief Functional verification for miracl-wrapper Tools library.
 */

#include "../include/Tools.h"
#include <iostream>
#include <string>
#include <vector>

using namespace std;

#define TEST_PASS(msg) do { cout << "[PASS] " << msg << endl; } while(0)
#define TEST_FAIL(msg) do { cout << "[FAIL] " << msg << endl; exit(1); } while(0)

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
        mpz_class res_mpz = (a_mpz + b_mpz) % q_mpz;

        BIG res_big;
        BIG_modadd(res_big, a_big, b_big, order);

        mpz_class res_converted = BIG_to_mpz(res_big);
        if (res_mpz == res_converted) TEST_PASS("Mod Addition (a + b) % q");
        else TEST_FAIL("Mod Addition mismatch");
    }

    // B. Modular Subtraction
    {
        mpz_class res_mpz = (a_mpz - b_mpz) % q_mpz;
        if (res_mpz < 0) res_mpz += q_mpz;

        BIG res_big, neg_b;
        BIG_modneg(neg_b, b_big, order);
        BIG_modadd(res_big, a_big, neg_b, order);

        mpz_class res_converted = BIG_to_mpz(res_big);
        if (res_mpz == res_converted) TEST_PASS("Mod Subtraction (a - b) % q");
        else TEST_FAIL("Mod Subtraction mismatch");
    }

    // C. Modular Multiplication
    {
        mpz_class res_mpz = (a_mpz * b_mpz) % q_mpz;

        BIG res_big;
        BIG_modmul(res_big, a_big, b_big, order);

        mpz_class res_converted = BIG_to_mpz(res_big);
        if (res_mpz == res_converted) TEST_PASS("Mod Multiplication (a * b) % q");
        else TEST_FAIL("Mod Multiplication mismatch");
    }
}

// ==================================================================
// 2. Conversion & ECP Adapter Test
// ==================================================================
void Test_Conversion_And_ECP_Adapter() {
    cout << "\n--- Test 2: Conversions & ECP Adapters ---" << endl;

    initRNG(&rng_tools);

    // A. Round-Trip Conversion
    BIG big_original, big_recovered;
    randBig(big_original, rng_tools);

    mpz_class mpz_val = BIG_to_mpz(big_original);
    mpz_to_BIG(mpz_val, big_recovered);

    if (BIG_comp(big_original, big_recovered) == 0) TEST_PASS("Type Conversion (BIG <-> mpz)");
    else TEST_FAIL("Type Conversion mismatch");

    // B. ECP Adapter (G1 scalar mul)
    ECP P_raw, P_wrapper;
    ECP_generator(&P_raw);
    ECP_copy(&P_wrapper, &P_raw);

    BIG r_big;
    randBig(r_big, rng_tools);
    mpz_class r_mpz = BIG_to_mpz(r_big);

    ECP_mul(&P_raw, r_big);
    ECP_mul(P_wrapper, r_mpz);

    if (ECP_equals(&P_raw, &P_wrapper)) TEST_PASS("ECP_mul adapter matches native");
    else TEST_FAIL("ECP_mul adapter logic failed");

    // C. ECP2 Adapter (G2 scalar mul)
    ECP2 P2_raw, P2_wrapper;
    ECP2_generator(&P2_raw);
    ECP2_copy(&P2_wrapper, &P2_raw);

    ECP2_mul(&P2_raw, r_big);
    ECP2_mul(P2_wrapper, r_mpz);

    if (ECP2_equals(&P2_raw, &P2_wrapper)) TEST_PASS("ECP2_mul adapter matches native");
    else TEST_FAIL("ECP2_mul adapter logic failed");
}

// ==================================================================
// 3. FP12 Adapter Test
// ==================================================================
void Test_FP12_Adapter() {
    cout << "\n--- Test 3: FP12 Functions ---" << endl;

    initRNG(&rng_tools);

    ECP g1; ECP_generator(&g1);
    ECP2 g2; ECP2_generator(&g2);

    // A. Pairing
    FP12 gt = e(g1, g2);
    if (!FP12_isunity(&gt)) TEST_PASS("Pairing calculation e(g1, g2)");
    else TEST_FAIL("Pairing failed");

    // B. Power Adapter
    BIG exp_big;
    randBig(exp_big, rng_tools);
    mpz_class exp_mpz = BIG_to_mpz(exp_big);

    FP12 res_raw, res_wrapper;
    FP12_copy(&res_raw, &gt);
    FP12_copy(&res_wrapper, &gt);

    FP12_pow(&res_raw, &res_raw, exp_big);
    FP12_pow(res_wrapper, exp_mpz);

    if (FP12_equals(&res_raw, &res_wrapper)) TEST_PASS("FP12_pow adapter matches native");
    else TEST_FAIL("FP12_pow adapter logic failed");
}

// ==================================================================
// 4. Variadic Hash Test (HashToZp, HashToG1, HashToG2)
// ==================================================================
void Test_Variadic_Hash() {
    cout << "\n--- Test 4: Variadic Hash Functions (Zp, G1, G2) ---" << endl;

    initRNG(&rng_tools);

    // Prepare elements of various types
    string msg = "Miracl_Wrapper_OpenSource_v2";
    mpz_class num = mpz_class("ABCDEF123456", 16);

    BIG b1; randBig(b1, rng_tools);
    ECP g1; ECP_generator(&g1);
    ECP2 g2; ECP2_generator(&g2);
    FP12 gt = e(g1, g2);

    // Prepare vectors for testing generic template
    vector<ECP> vec_ecp = {g1, g1};
    vector<ECP2> vec_ecp2 = {g2, g2};

    // A. HashToZp: Consistency and Mixed Type Check
    {
        mpz_class h1 = HashToZp(msg, num, b1, g1, g2, gt);
        mpz_class h1_check = HashToZp(msg, num, b1, g1, g2, gt);

        if (h1 == h1_check && h1 > 0) {
            TEST_PASS("HashToZp: Mixed types (Consistency check)");
        } else {
            TEST_FAIL("HashToZp: Consistency check failed");
        }
    }

    // B. HashToG1: Map mixed types to G1 point
    {
        ECP pt1, pt2;
        HashToG1(pt1, msg, g1, vec_ecp);
        HashToG1(pt2, msg, g1, vec_ecp);

        if (ECP_equals(&pt1, &pt2) && !ECP_isinf(&pt1)) {
            TEST_PASS("HashToG1: Mixed types -> G1 (Consistency check)");
        } else {
            TEST_FAIL("HashToG1: Mapping failed or inconsistent");
        }
    }

    // C. HashToG2: Map mixed types to G2 point (New Test)
    {
        ECP2 pt1, pt2;
        // Test with different order and nested-like structure (vector)
        HashToG2(pt1, num, vec_ecp2, gt, "G2_Test_Salt");
        HashToG2(pt2, num, vec_ecp2, gt, "G2_Test_Salt");

        if (ECP2_equals(&pt1, &pt2) && !ECP2_isinf(&pt1)) {
            TEST_PASS("HashToG2: Mixed types -> G2 (Consistency check)");
        } else {
            TEST_FAIL("HashToG2: Mapping failed or inconsistent");
        }
    }

    // D. Cross-Group Independence Check
    {
        // Ensure that hashing same input to Zp, G1, G2 results in "matching" but distinct objects
        mpz_class h_val = HashToZp(msg);
        ECP h_g1; HashToG1(h_g1, msg);
        ECP2 h_g2; HashToG2(h_g2, msg);

        // Verification: If we manually multiply generator by HashToZp result, it should match HashToG1
        ECP g1_gen, g1_manual;
        ECP_generator(&g1_gen);
        ECP_copy(&g1_manual, &g1_gen);
        ECP_mul(g1_manual, h_val);

        if (ECP_equals(&h_g1, &g1_manual)) {
            TEST_PASS("Hash Consistency: HashToG1 matches generator * HashToZp");
        } else {
            TEST_FAIL("Hash Consistency: Cross-group logic mismatch");
        }
    }
}

int main() {
    cout << "=== Running Wrapper Verification ===" << endl;

    mpz_class order = getCurveOrder();
    show_mpz(order.get_mpz_t());
    Test_GMP_Convenience();
    Test_Conversion_And_ECP_Adapter();
    Test_FP12_Adapter();
    Test_Variadic_Hash(); // New Test suite added here

    cout << "\n=== All Tests Passed ===" << endl;
    return 0;
}