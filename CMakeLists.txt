cmake_minimum_required(VERSION 3.28)
project(miracl-wrapper VERSION 1.0.0 LANGUAGES C CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# 1. 处理依赖
# ---------------------------------------------------------
# 移除 cmake module path，因为我们要直接编译源码
# list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

# 重要：这里会进入 external 目录，你需要确保 external/CMakeLists.txt 里包含了 gmp
add_subdirectory(external)

# 【已删除】 find_package(GMP REQUIRED)
# 我们不再查找系统库，而是直接使用 add_subdirectory 编译出来的目标

# 2. 定义库
# ---------------------------------------------------------
add_library(WrapperLib STATIC
        src/Tools.cpp
)

# 3. 设置 Include 路径
# ---------------------------------------------------------
target_include_directories(WrapperLib
        PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:include>
        PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/src
)

# 4. 链接依赖
# ---------------------------------------------------------
target_link_libraries(WrapperLib PUBLIC
        miracl_core
        # 这里直接链接子模块生成的 Target 名称。
        # 注意：你需要确认 rookie-papers/GMP 这个库生成的 Target 名字叫什么。
        # 如果它也是标准的 CMake 项目，通常是 gmp 或 gmp::gmp。
        # 假设它生成的静态库 target 叫 gmp
        gmp
        # 如果你的子模块没有生成 gmpxx，这行可能要去掉，或者看子模块是否支持 C++ 绑定
        # gmpxx
)


# =========================================================
# 测试模块配置
# =========================================================
# 开启测试支持 (这样才能用 ctest)
enable_testing()

# 添加 tests 子目录
add_subdirectory(tests)