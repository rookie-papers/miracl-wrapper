# MIRACL-Wrapper

**MIRACL-Wrapper** is a lightweight, high-performance C++ utility library designed to bridge the gap between **MIRACL Core** (Efficient Elliptic Curve Cryptography) and **GMP** (Arbitrary Precision Arithmetic).

This library encapsulates complex low-level memory management and data conversion, providing a modern C++ interface. It allows developers to use GMP's `mpz_class` directly with MIRACL's Elliptic Curve structures (`ECP`, `FP12`, etc.), significantly accelerating the prototyping and implementation of advanced cryptographic schemes.

---

## Key Features

* **Seamless Integration**: Perform scalar multiplication on Elliptic Curves using GMP integers (`mpz_class`) directly.
* **Automatic Conversion**: Handles bidirectional conversion between MIRACL's `BIG` type and GMP's `mpz_class` transparently.
* **Simplified API**: Provides easy-to-use wrappers for Bilinear Pairings,
* **Dependency Management**: Automatically manages the compilation of MIRACL Core and GMP as static libraries.

---

## ⚙️ Build & Test Instructions

The following instructions guide you through installing dependencies, building the library from source, and running the included unit tests.

### 1️⃣ Install Required Dependencies

Run the following command to install the necessary build tools.

Note: `m4` is strictly required for building the GMP library.

```bash
sudo apt update
sudo apt install -y git cmake python3 build-essential m4
```

### 2️⃣ Clone & Build
This repository uses Git submodules to manage dependencies (MIRACL Core, GMP, Google Benchmark). Recursive cloning is essential.

```bash
# Clone the repository with all submodules
git clone --recurse-submodules https://github.com/rookie-papers/miracl-wrapper.git
cd miracl-wrapper
mkdir build && cd build

# Configure and Build; Release mode is recommended for performance benchmarks
cmake -DCMAKE_BUILD_TYPE=Release ..
make
```
### 3️⃣ Running Verification & Benchmarks
After a successful build, the test executables are generated in the tests/ directory.

```bash
# Functional Verification: Runs unit tests to ensure math logic, data conversion, and wrappers are correct.
./tests/test_tools
Performance Benchmark: Measures the overhead of the wrapper compared to raw MIRACL calls.
./tests/test_benchmark
```

---


## 📦 Integration Guide
You can easily integrate miracl-wrapper into any CMake-based cryptographic project using FetchContent. This approach handles the download, compilation, and linking automatically.

Step 1: Configure CMakeLists.txt
Add the following configuration to your project's CMakeLists.txt:

```bash
include(FetchContent)

FetchContent_Declare(
    WrapperLib
    GIT_REPOSITORY https://github.com/rookie-papers/miracl-wrapper.git
    GIT_TAG        main
)

# Download and populate the library
FetchContent_MakeAvailable(WrapperLib)

# ... (Your executable definitions) ...

# Link your target against WrapperLib
# This automatically pulls in MIRACL Core, GMP, and the Tools headers
target_link_libraries(YourTargetName PRIVATE
    WrapperLib
    # benchmark::benchmark  <-- Optional: Add this line only if you use Google Benchmark
)
```


Step 2: Usage in C++ Code
You can include "Tools.h" and start using the combined features of GMP and MIRACL.

```cpp
#include "Tools.h"
#include <iostream>

void example_usage() {
    // 1. Initialize State
    initState(state_gmp);   // Initialize GMP random state
    initRNG(&rng);          // Initialize MIRACL random state

    // 2. Generate a random scalar using GMP
    mpz_class secret_key = rand_mpz(state_gmp);

    // 3. Perform Elliptic Curve Operations
    ECP P;
    ECP_generator(&P);      // P = Generator of G1

    // Wrapper allows passing mpz_class directly to ECP_mul
    ECP_mul(P, secret_key); 

    std::cout << "Scalar multiplication complete." << std::endl;
}
```