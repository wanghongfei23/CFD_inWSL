cmake_minimum_required(VERSION 3.12)
project(WCNS_Convection LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Set output directory for executables
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin)

# Include directories
include_directories(include)

# Quadmath library setup
find_library(QUADMATH_LIBRARY quadmath)
if(NOT QUADMATH_LIBRARY)
    message(FATAL_ERROR "libquadmath not found - install gcc libquadmath-dev")
endif()

# Find the OpenMP package
find_package(OpenMP REQUIRED)

# Executable with time
add_executable(wcns_convection src/with_time.cpp)

target_compile_options(wcns_convection PRIVATE 
    -Wall -Wextra -O0 -march=native
)

target_link_libraries(wcns_convection PRIVATE 
    ${QUADMATH_LIBRARY} 
    stdc++ 
    m  # Math library
    OpenMP::OpenMP_CXX
)

# For quadmath output support
target_compile_definitions(wcns_convection PRIVATE 
    HAVE_QUADMATH_H
)

# Executable without time
add_executable(wcns_convection_without_time src/without_time.cpp)

target_compile_options(wcns_convection_without_time PRIVATE 
    -Wall -Wextra -O0 -march=native
)

target_link_libraries(wcns_convection_without_time PRIVATE 
    ${QUADMATH_LIBRARY} 
    stdc++ 
    m  # Math library
    OpenMP::OpenMP_CXX
)

# For quadmath output support
target_compile_definitions(wcns_convection_without_time PRIVATE 
    HAVE_QUADMATH_H
)