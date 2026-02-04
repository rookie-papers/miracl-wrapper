# =========================================================
# 1. 功能测试 (Correctness Test)
# =========================================================
add_executable(test_tools test_tools.cpp)

# 链接主库 WrapperLib (它自动包含了 MIRACL 和 GMP)
target_link_libraries(test_tools PRIVATE WrapperLib)

# 注册到 CTest (方便用 ctest 命令运行)
add_test(NAME FunctionalTest COMMAND test_tools)


# =========================================================
# 2. 基准测试 (Performance Benchmark)
# =========================================================
add_executable(test_benchmark test_benchmark.cpp)

# 链接主库 + Google Benchmark
target_link_libraries(test_benchmark
        PRIVATE
        WrapperLib
        benchmark::benchmark      # 必须显式链接 benchmark
        benchmark::benchmark_main # (可选) 如果你不想写 main 函数可以用这个，但你写了
)

# Benchmark 通常不注册为 CTest 的一部分，因为它是看输出的，不是单纯 Pass/Fail
# 但如果你想确保它能跑通，也可以加个 test
add_test(NAME PerformanceRun COMMAND test_benchmark)