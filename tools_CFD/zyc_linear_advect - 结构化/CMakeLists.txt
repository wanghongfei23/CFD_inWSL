#=============================================================================
# 项目基本配置
#=============================================================================

# 指定所需的最低 CMake 版本
cmake_minimum_required(VERSION 3.12)

# 定义项目名称和使用的编程语言
project(WCNS_Convection LANGUAGES CXX)

# 设置 C++ 标准版本为 C++23
set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# 设置可执行文件的输出目录为项目根目录下的 bin 目录
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin)


#=============================================================================
# 头文件路径配置
#=============================================================================

# 添加头文件搜索路径，包含 include 目录
include_directories(include)


#=============================================================================
# 第三方库配置
#=============================================================================

# 查找 quadmath 数学库
find_library(QUADMATH_LIBRARY quadmath)
if(NOT QUADMATH_LIBRARY)
    message(FATAL_ERROR "libquadmath未找到 - 请安装gcc libquadmath-dev")
endif()

# 查找并链接 OpenMP 并行计算库
find_package(OpenMP REQUIRED)


#=============================================================================
# 可执行文件配置 - 带时间维度的WCNS对流求解器
#=============================================================================

# 定义可执行文件目标
add_executable(wcns_convection src/with_time.cpp)

# 设置编译选项
target_compile_options(wcns_convection PRIVATE 
    -Wall -Wextra -O0 -march=native
)

# 链接所需库文件
target_link_libraries(wcns_convection PRIVATE 
    ${QUADMATH_LIBRARY} 
    stdc++ 
    m  # 数学库
    OpenMP::OpenMP_CXX
)

# 添加预处理器定义以支持 quadmath 输出
target_compile_definitions(wcns_convection PRIVATE 
    HAVE_QUADMATH_H
)


#=============================================================================
# 可执行文件配置 - 不带时间维度的WCNS对流求解器
#=============================================================================

# 定义可执行文件目标
add_executable(wcns_convection_without_time src/without_time.cpp)

# 设置编译选项
target_compile_options(wcns_convection_without_time PRIVATE 
    -Wall -Wextra -O0 -march=native
)

# 链接所需库文件
target_link_libraries(wcns_convection_without_time PRIVATE 
    ${QUADMATH_LIBRARY} 
    stdc++ 
    m  # 数学库
    OpenMP::OpenMP_CXX
)

# 添加预处理器定义以支持 quadmath 输出
target_compile_definitions(wcns_convection_without_time PRIVATE 
    HAVE_QUADMATH_H
)