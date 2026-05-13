#pragma once
#include <iostream>
#include <pair_BLS12381.h>
#include <bls_BLS12381.h>
#include <randapi.h>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <vector>
#include <sys/time.h>
#include "gmpxx.h"
#include <sstream>
#include <iomanip>
#include <chrono>
#include <random>
#include <cassert>
#include <string>
#include <stdexcept>

using namespace B384_58;
using namespace BLS12381;
using namespace core;
using namespace std;
using namespace chrono;

/**
 * Initializes the random seed for generating random numbers
 * @param rng Seed to be initialized
 */
void initRNG(csprng *rng);

mpz_class getCurveOrder();

/**
 * Generates a random 256-bit number in the Fq field, where q is the order of the BLS12381 curve.
 * Note: Run initRNG before calling this function to initialize the seed.
 * @param big The variable to store the generated random number
 * @param rng Seed used for random generation (initialized by initRNG)
 */
void randBig(BIG big, csprng &rng);

/**
 * Generates a random point on the G1 group of the BLS12381 curve
 * @param ecp Output G1 point
 */
ECP randECP(csprng &rng);

/**
 * Generates a random point on the G2 group of the BLS12381 curve
 * @param ecp Output G2 point
 */
ECP2 randECP2(csprng &rng);

/**
 * Converts data of type char* to a string
 * @param ch Data to be converted of type char*
 * @return The converted string
 */
string charsToString(char *ch);

/**
 * Converts a BIG integer to an mpz_class integer
 * @param big BIG integer to be converted
 * @return Converted mpz_class integer
 */
mpz_class BIG_to_mpz(BIG big);

/**
 * Converts an mpz_class integer to a BIG integer
 * @param t mpz_class integer to be converted
 * @param big Output BIG integer
 */
void mpz_to_BIG(const mpz_class& t, BIG& big);

/**
 * Converts a string to a BIG integer
 * @param hex_string The string to be converted
 * @param big Output BIG integer
 */
void str_to_BIG(string hex_string, BIG &big);

/**
 * Elliptic curve multiplication with constant using mpz_class
 * @param P1 Elliptic curve point
 * @param t The multiplier
 */
void ECP_mul(ECP& P1, const mpz_class& t);

/**
 * Elliptic curve multiplication with constant using BIG
 * @param P2 Elliptic curve point
 * @param t The multiplier
 */
void ECP2_mul(ECP2& P2, const mpz_class& t);

void FP12_mulMy(FP12 &a, FP12 &b);

void FP12_pow(FP12 &r, const mpz_class &exp);

void FP12_inv(FP12 &r);

/**
 * Initializes the random number generator in GMP
 * @param state Random number generator to be initialized
 */
void initState(gmp_randstate_t& state);

/**
 * Generates a random mpz_class integer
 * @param bits Number of bits for the generated integer
 * @param state Random number generator
 * @return Randomly generated mpz_class integer
 */
mpz_class rand_mpz(gmp_randstate_t state);


/**
 * Computes base^exp % mod, with the result stored in res
 * @param base Base number
 * @param exp Exponent
 * @param mod Modulus
 * @return Result of modular exponentiation
 */
mpz_class pow_mpz(const mpz_class& base, const mpz_class& exp, const mpz_class& mod);

/**
 * Computes the modular multiplicative inverse of a under modulo m, result stored in res
 * @param a Integer to find the inverse of
 * @param m Modulus
 * @return Modular multiplicative inverse, or 0 if it doesn't exist
 */
mpz_class invert_mpz(const mpz_class& a, const mpz_class& m);

/**
 * Computes the Lagrange interpolation coefficients for a given set of points under a modulus.
 * @param x Vector of x-coordinates of the interpolation points.
 * @param y Vector of y-coordinates of the interpolation points.
 * @param modulus Modulus for computations (all operations are performed modulo `modulus`).
 * @return Vector of coefficients of the Lagrange interpolating polynomial under the given modulus.
 *         The coefficients are ordered from the constant term to the highest degree term.
 */
vector<mpz_class> getLagrangeCoffs(const vector<mpz_class> &x, const vector<mpz_class> &y, const mpz_class &modulus);

/**
 * Evaluates a polynomial at a given point under a modulus.
 * @param poly Vector of polynomial coefficients, ordered from the constant term to the highest degree term.
 * @param x Point at which to evaluate the polynomial.
 * @param modulus Modulus for computations (all operations are performed modulo `modulus`).
 * @return The value of the polynomial at the given point `x`, computed modulo `modulus`.
 */
mpz_class computePoly(const vector<mpz_class> &poly, const mpz_class &x, const mpz_class &modulus);

/**
 * Computes the values of the Lagrange basis polynomials at zero for a given set of x-coordinates under a modulus.
 * @param x Vector of x-coordinates of the interpolation points.
 * @param q Modulus for computations (all operations are performed modulo `q`).
 * @return Vector of values, where each value corresponds to the evaluation of a Lagrange basis polynomial
 *         at zero. The result is computed modulo `q`.
 */
vector<mpz_class> getLagrangeBasis(const vector<mpz_class> &x, const mpz_class &q);

/**
 * Outputs a mpz_t integer in hexadecimal format, including a newline
 * @param mpz Integer to view
 */
void show_mpz(mpz_t mpz);

/**
 * Initializes an octet with a buffer of specified maximum length.
 * @param maxLen Maximum buffer size in bytes.
 * @return An initialized octet structure.
 */
octet getOctet(int maxLen);

/**
 * @brief Prints the contents of an octet as a hexadecimal string.
 * @param S The octet to be printed.
 */
void showOctet(const octet *S);

octet mpzToOctet(const mpz_class& num);

mpz_class octetToMpz(const octet& o);

/**
 * @brief Concatenates two octets into a new octet and returns it.
 *
 * @param oc1 The first octet.
 * @param oc2 The second octet.
 * @return octet A new octet containing the concatenated data.
 *               Ensure to free the returned octet's val after use.
 */
octet concat_Octet(const octet *oc1, const octet *oc2);

/**
 * @brief Concatenates the second octet to the first octet, resizing the buffer if needed.
 *
 * This function appends the content of `oc2` to `oc1`. If `oc1`'s buffer is insufficient,
 * it automatically resizes the buffer using `realloc`.
 *
 * @param oc1 Pointer to the octet that stores the result. Its buffer may be resized.
 * @param oc2 Pointer to the octet to be appended.
 * @return true If the concatenation is successful.
 * @return false If inputs are invalid or memory allocation fails.
 *
 * @note Ensure to free `oc1->val` after use if resizing occurs.
 */
bool concatOctet(octet *oc1, const octet *oc2);

/**
 * Hashes an octet to a 256-bit integer, used by hashToZp256
 * @param num Hash result
 * @param ct Octet to hash
 * @param q Order of the elliptic curve to mod the hash result
 */
void hashZp256(BIG res, octet *ct, BIG q);

/**
 * Bilinear pairing
 * @param alpha1 Element on G1
 * @param alpha2 Element on G2
 * @return Result of bilinear pairing, an element on GT
 */
FP12 e(ECP P1, ECP2 P2);

/**
 * Outputs a BIG integer in hexadecimal format, including a newline
 * @param big Integer to view
 */
void showBIG(BIG big);

/**
 * Outputs an FP12 element of the BLS12381 elliptic curve GT in hexadecimal format, including a newline
 * @param fp12 Element on GT to view
 */
void showFP12(FP12 fp12);

/**
 * Outputs a separator line for debugging purposes
 * @param text Text content of the separator line
 */
void printLine(const string& text);




// ============================================================================
// Variadic Template Hash Functions
// ============================================================================


// ------------------------------------------------------------------
// Append string
inline void appendToOctet(octet* hash, const string& str) {
    octet temp = getOctet(1024);
    temp.len = str.length();
    memcpy(temp.val, str.c_str(), temp.len);
    concatOctet(hash, &temp);
    free(temp.val);
}

// Append mpz_class
inline void appendToOctet(octet* hash, const mpz_class& num) {
    appendToOctet(hash, num.get_str(16));
}

// Append MIRACL BIG type
inline void appendToOctet(octet* hash, BIG b) {
    appendToOctet(hash, BIG_to_mpz(b));
}

// Append G1 element (ECP)
inline void appendToOctet(octet* hash, ECP pt) {
    octet temp = getOctet(1024);
    ECP_toOctet(&temp, &pt, true);
    concatOctet(hash, &temp);
    free(temp.val);
}

// Append G2 element (ECP2)
inline void appendToOctet(octet* hash, ECP2 pt) {
    octet temp = getOctet(1024);
    ECP2_toOctet(&temp, &pt, true);
    concatOctet(hash, &temp);
    free(temp.val);
}

// Append GT element (FP12)
inline void appendToOctet(octet* hash, FP12 pt) {
    octet temp = getOctet(1024);
    FP12_toOctet(&temp, &pt);
    concatOctet(hash, &temp);
    free(temp.val);
}

// Generic support for std::vector of any type T
template<typename T>
inline void appendToOctet(octet* hash, const vector<T>& vec) {
    for (const auto& item : vec) {
        // Automatically resolves to the correct appendToOctet based on T
        appendToOctet(hash, item);
    }
}

// ------------------------------------------------------------------


// Base case to terminate recursion
inline void hash_recursive(octet* hash) {}

// Recursive step
template<typename T, typename... Args>
inline void hash_recursive(octet* hash, T first, Args... rest) {
    appendToOctet(hash, first);        // Process the first argument
    hash_recursive(hash, rest...);     // Recursively process the rest
}

// ------------------------------------------------------------------

/**
 * @brief Variadic hash function that maps any number of arguments to a Zp scalar.
 * @tparam Args Template parameter pack for various types (mpz_class, ECP, ECP2, FP12, string, vector, etc.)
 * @param args The sequence of arguments to be hashed.
 * @return The resulting hash as a GMP mpz_class integer modulo the curve order.
 */
template<typename... Args>
inline mpz_class HashToZp(Args... args) {
    octet hash = getOctet(8192);
    hash_recursive(&hash, args...);
    BIG order, ret;
    BIG_rcopy(order, CURVE_Order);
    hashZp256(ret, &hash, order);
    free(hash.val);
    return BIG_to_mpz(ret);
}

/**
 * @brief Variadic hash function that maps any number of arguments to a point on the G1 group.
 * @tparam Args Template parameter pack for various types.
 * @param res Output parameter: the resulting ECP point on G1.
 * @param args The sequence of arguments to be hashed.
 */
template<typename... Args>
inline void HashToG1(ECP& res, Args... args) {
    octet hash = getOctet(8192);
    hash_recursive(&hash, args...);
    BIG order, scalar;
    BIG_rcopy(order, CURVE_Order);
    hashZp256(scalar, &hash, order);

    // Map scalar to G1: res = scalar * G1_generator
    ECP_generator(&res);
    ECP_mul(&res, scalar);

    free(hash.val);
}

/**
 * @brief Variadic hash function that maps any number of arguments to a point on the G2 group.
 * @tparam Args Template parameter pack for various types.
 * @param res Output parameter: the resulting ECP2 point on G2.
 * @param args The sequence of arguments to be hashed.
 */
template<typename... Args>
inline void HashToG2(ECP2& res, Args... args) {
    octet hash = getOctet(8192);
    hash_recursive(&hash, args...);
    BIG order, scalar;
    BIG_rcopy(order, CURVE_Order);
    hashZp256(scalar, &hash, order);

    // Map scalar to G2: res = scalar * G2_generator
    ECP2_generator(&res);
    ECP2_mul(&res, scalar);

    free(hash.val);
}

namespace GenericSerializer {

std::string binToHex(const std::string& input);
std::string hexToBin(const std::string& input);

// ============================================================================
// 1. Base Type Converters Declarations
// Implemented in the .cpp file. Compiler automatically deduces the correct 
// overloaded function based on the parameter type.
// ============================================================================

// --- Cryptographic Types ---
std::string encode(const mpz_class& val);
void decode(const std::string& str, mpz_class& val);

std::string encode(const ECP& val);
void decode(const std::string& str, ECP& val);

std::string encode(const ECP2& val);
void decode(const std::string& str, ECP2& val);

std::string encode(const FP12& val);
void decode(const std::string& str, FP12& val);

// --- Standard C++ Types ---
std::string encode(const std::string& val);
void decode(const std::string& str, std::string& val);

std::string encode(int val);
void decode(const std::string& str, int& val);

std::string encode(short val);
void decode(const std::string& str, short& val);

std::string encode(long long val);
void decode(const std::string& str, long long& val);

std::string encode(bool val);
void decode(const std::string& str, bool& val);

// --- Utility Functions ---
/**
 * @brief Splits a string into a vector of substrings based on a specified delimiter.
 */
std::vector<std::string> split(const std::string& s, char delimiter);


// ============================================================================
// 2. Vector Support (Template Implementation)
// Uses a comma ',' as the delimiter for elements within a vector.
// ============================================================================

template <typename T>
std::string encode(const std::vector<T>& vec) {
    std::ostringstream oss;
    for (size_t i = 0; i < vec.size(); ++i) {
        oss << encode(vec[i]);
        if (i != vec.size() - 1) {
            oss << ",";
        }
    }
    return oss.str();
}

template <typename T>
void decode(const std::string& str, std::vector<T>& vec) {
    vec.clear();
    if (str.empty()) return;
    
    std::vector<std::string> items = split(str, ',');
    for (const auto& item : items) {
        if (item.empty()) continue;
        T val;
        decode(item, val);
        vec.push_back(val);
    }
}

// ============================================================================
// 3. Variadic Template Serialization (Core Logic)
// Uses '#' as the delimiter between different fields.
// ============================================================================

// Base case for serialization recursion
inline void serialize_recursive(std::ostringstream& oss) {}

// Recursive step for serialization
template<typename T, typename... Args>
inline void serialize_recursive(std::ostringstream& oss, const T& first, const Args&... rest) {
    oss << encode(first);
    if (sizeof...(rest) > 0) {
        oss << "#";
    }
    serialize_recursive(oss, rest...);
}

/**
 * @brief Universal serialization function supporting an arbitrary number 
 * of parameters of supported types.
 * @return A single concatenated string representing the serialized data.
 */
template<typename... Args>
std::string serialize(const Args&... args) {
    std::ostringstream oss;
    serialize_recursive(oss, args...);
    return oss.str();
}

// ============================================================================
// 4. Variadic Template Deserialization (Core Logic)
// ============================================================================

// Base case for deserialization recursion
inline void deserialize_recursive(const std::vector<std::string>& fields, size_t& index) {}

// Recursive step for deserialization
template<typename T, typename... Args>
inline void deserialize_recursive(const std::vector<std::string>& fields, size_t& index, T& first, Args&... rest) {
    if (index >= fields.size()) {
        throw std::runtime_error("Deserialization error: Not enough fields in the input string.");
    }
    decode(fields[index++], first);
    deserialize_recursive(fields, index, rest...);
}

/**
 * @brief Universal deserialization function. 
 * Assigns parsed values directly to the passed references.
 * @param data The serialized string.
 * @param args The variables (passed by reference) to populate.
 */
template<typename... Args>
void deserialize(const std::string& data, Args&... args) {
    std::vector<std::string> fields = split(data, '#');
    size_t index = 0;
    deserialize_recursive(fields, index, args...);
}

} // namespace GenericSerializer