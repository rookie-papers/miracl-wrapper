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

// ==================================================================
// 5. Generic Serializer Test
// ==================================================================

struct DummyTestStruct {
    int _int; short _short; long long _ll; bool _bool; std::string _str;
    mpz_class _mpz; ECP _g1; ECP2 _g2; FP12 _gt;
    std::vector<int> _v_int; std::vector<std::string> _v_str; std::vector<mpz_class> _v_mpz;
    std::vector<ECP> _v_g1; std::vector<ECP2> _v_g2; std::vector<FP12> _v_gt;
};

std::string encode(const DummyTestStruct& obj) {
    return GenericSerializer::serialize(
        obj._int, obj._short, obj._ll, obj._bool, obj._str,
        obj._mpz, obj._g1, obj._g2, obj._gt,
        obj._v_int, obj._v_str, obj._v_mpz, obj._v_g1, obj._v_g2, obj._v_gt
    );
}

void decode(const std::string& str, DummyTestStruct& obj) {
    GenericSerializer::deserialize(str,
        obj._int, obj._short, obj._ll, obj._bool, obj._str,
        obj._mpz, obj._g1, obj._g2, obj._gt,
        obj._v_int, obj._v_str, obj._v_mpz, obj._v_g1, obj._v_g2, obj._v_gt
    );
}

// 5.3 Main function for serialization comprehensive testing
void Test_Generic_Serializer() {
    cout << "\n--- Test 5: Generic Serialization (Basic, Crypto, Vector, Struct) ---" << endl;

    // A. Basic Types Test
    {
        int o_int = -12345; short o_short = 255; long long o_ll = 987654321012345LL;
        bool o_bool = true; string o_str = "CryptoSerializationTest";

        string s1 = GenericSerializer::serialize(o_int, o_short, o_ll, o_bool, o_str);
        int r_int; short r_short; long long r_ll; bool r_bool; string r_str;
        GenericSerializer::deserialize(s1, r_int, r_short, r_ll, r_bool, r_str);
        string s2 = GenericSerializer::serialize(r_int, r_short, r_ll, r_bool, r_str);

        if (s1 == s2) { TEST_PASS("Serializer: Basic Types"); } 
        else { TEST_FAIL("Serializer: Basic Types mismatch"); }
    }

    // B. Cryptography Types Test
    {
        mpz_class o_mpz("FFFFFFFFFFFFFFFFFFFFFFFFFFFFFF", 16);
        ECP o_g1; ECP_generator(&o_g1); 
        ECP2 o_g2; ECP2_generator(&o_g2); 
        FP12 o_gt; FP12_one(&o_gt); 

        string s1 = GenericSerializer::serialize(o_mpz, o_g1, o_g2, o_gt);
        mpz_class r_mpz; ECP r_g1; ECP2 r_g2; FP12 r_gt;
        GenericSerializer::deserialize(s1, r_mpz, r_g1, r_g2, r_gt);
        string s2 = GenericSerializer::serialize(r_mpz, r_g1, r_g2, r_gt);

        if (s1 == s2) { TEST_PASS("Serializer: Cryptography Types"); } 
        else { TEST_FAIL("Serializer: Cryptography Types mismatch"); }
    }

    // C. Vector Types Test
    {
        vector<int> o_v_int = {1, 2, 3};
        vector<string> o_v_str = {"Alice", "Bob"};
        ECP g1; ECP_generator(&g1); ECP2 g2; ECP2_generator(&g2); FP12 gt; FP12_one(&gt);

        vector<mpz_class> o_v_mpz = { mpz_class(10), mpz_class("10000000000", 16) };
        vector<ECP> o_v_g1 = {g1, g1};
        vector<ECP2> o_v_g2 = {g2, g2};
        vector<FP12> o_v_gt = {gt, gt};

        string s1 = GenericSerializer::serialize(o_v_int, o_v_str, o_v_mpz, o_v_g1, o_v_g2, o_v_gt);

        vector<int> r_v_int; vector<string> r_v_str; vector<mpz_class> r_v_mpz;
        vector<ECP> r_v_g1; vector<ECP2> r_v_g2; vector<FP12> r_v_gt;
        GenericSerializer::deserialize(s1, r_v_int, r_v_str, r_v_mpz, r_v_g1, r_v_g2, r_v_gt);
        string s2 = GenericSerializer::serialize(r_v_int, r_v_str, r_v_mpz, r_v_g1, r_v_g2, r_v_gt);

        if (s1 == s2) { TEST_PASS("Serializer: Vector Types (Including Crypto)"); } 
        else { TEST_FAIL("Serializer: Vector Types mismatch"); }
    }

    // D. Custom Struct Test
    {
        DummyTestStruct obj1;
        obj1._int = 1001; obj1._short = 255; obj1._ll = 999999999LL;
        obj1._bool = true; obj1._str = "ZKP_Protocol";
        obj1._mpz = mpz_class("ABCDEF123456", 16);
        ECP_generator(&obj1._g1); ECP2_generator(&obj1._g2); FP12_one(&obj1._gt);
        obj1._v_int = {1, 2, 3}; obj1._v_str = {"UAV", "Auth"};
        obj1._v_mpz = {mpz_class(1), mpz_class(2)};
        obj1._v_g1 = {obj1._g1}; obj1._v_g2 = {obj1._g2}; obj1._v_gt = {obj1._gt};

        string s1 = ::encode(obj1);
        DummyTestStruct obj2;
        ::decode(s1, obj2);
        string s2 = ::encode(obj2);

        if (s1 == s2) { TEST_PASS("Serializer: Ultimate Custom Struct (15 Types)"); } 
        else { TEST_FAIL("Serializer: Custom Struct mismatch"); }
    }

    // E. Edge Cases
    {
        string o_empty_str = "";
        vector<int> o_empty_vec;
        bool o_false = false;

        string s1 = GenericSerializer::serialize(o_empty_str, o_empty_vec, o_false);
        string r_empty_str = "dirty_data"; 
        vector<int> r_empty_vec = {1, 2, 3}; 
        bool r_false = true;

        GenericSerializer::deserialize(s1, r_empty_str, r_empty_vec, r_false);
        string s2 = GenericSerializer::serialize(r_empty_str, r_empty_vec, r_false);

        if (s1 == s2) { TEST_PASS("Serializer: Edge Cases (Empty data)"); } 
        else { TEST_FAIL("Serializer: Edge Cases mismatch"); }
    }
}

// ==================================================================
// 6. Nested Struct & Vector of Structs Test (嵌套结构体测试)
// ==================================================================

// 6.1 定义内层结构体
struct UAV_Info {
    int uav_id;
    mpz_class serial_number;
    ECP public_key;
};

// 6.2 为内层结构体提供序列化图纸 (底层库已自动加装甲，不再需要手动转 Hex)
std::string encode(const UAV_Info& uav) {
    return GenericSerializer::serialize(uav.uav_id, uav.serial_number, uav.public_key); 
}

void decode(const std::string& str, UAV_Info& uav) {
    GenericSerializer::deserialize(str, uav.uav_id, uav.serial_number, uav.public_key);
}
// ------------------------------------------------------------------

// 6.3 定义外层结构体 (包含内层结构体，以及内层结构体的 Vector 容器)
struct Fleet_Mission {
    std::string mission_name;
    UAV_Info leader;                 // 嵌套单体结构体
    std::vector<UAV_Info> followers; // 嵌套结构体的容器！(终极测试)
    mpz_class timestamp;
};

// 6.4 为外层结构体提供序列化图纸 (写法与内层完全一致)
std::string encode(const Fleet_Mission& fm) {
    return GenericSerializer::serialize(fm.mission_name, fm.leader, fm.followers, fm.timestamp);
}

void decode(const std::string& str, Fleet_Mission& fm) {
    GenericSerializer::deserialize(str, fm.mission_name, fm.leader, fm.followers, fm.timestamp);
}

// ------------------------------------------------------------------

// 6.5 测试函数 (测试逻辑保持不变，验证底层自动装甲是否生效)
void Test_Nested_Struct() {
    cout << "\n--- Test 6: Nested Structs & Struct Vectors (Auto-Hex Armor Version) ---" << endl;

    // A. 构造复杂的嵌套测试数据
    Fleet_Mission original_mission;
    original_mission.mission_name = "Operation_Nightfall";
    original_mission.timestamp = mpz_class("1688123456789");

    // 赋值 Leader
    original_mission.leader.uav_id = 1;
    original_mission.leader.serial_number = mpz_class("AAA111", 16);
    ECP_generator(&original_mission.leader.public_key);

    // 赋值 Followers (包含两个 UAV_Info)
    UAV_Info f1, f2;
    f1.uav_id = 2; f1.serial_number = mpz_class("BBB222", 16); ECP_generator(&f1.public_key);
    f2.uav_id = 3; f2.serial_number = mpz_class("CCC333", 16); ECP_generator(&f2.public_key);
    original_mission.followers = {f1, f2};

    // B. 一键序列化整个外层结构体
    string serialized_data = ::encode(original_mission);
    
    // 打印出来看看终极装甲的样子 (可选，用于调试验证)
    // cout << "[DEBUG] Serialized Hex Payload: " << serialized_data << endl;
    
    // C. 反序列化到新对象中
    Fleet_Mission recovered_mission;
    ::decode(serialized_data, recovered_mission);

    // D. 重新序列化比对
    string secondary_data = ::encode(recovered_mission);

    if (serialized_data == secondary_data) { 
        TEST_PASS("Serializer: Nested Structs & Struct Vectors"); 
    } else { 
        TEST_FAIL("Serializer: Nested Structs mismatch"); 
    }
}

// ==================================================================
// Main Execution
// ==================================================================
int main() {
    cout << "=== Running Wrapper Verification ===" << endl;

    mpz_class order = getCurveOrder();
    show_mpz(order.get_mpz_t());
    Test_GMP_Convenience();
    Test_Conversion_And_ECP_Adapter();
    Test_FP12_Adapter();
    Test_Variadic_Hash();
    Test_Generic_Serializer(); 
    Test_Nested_Struct();

    cout << "\n=== All Tests Passed ===" << endl;
    return 0;
}