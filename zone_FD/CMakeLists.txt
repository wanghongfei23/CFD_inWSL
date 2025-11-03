cmake_minimum_required(VERSION 3.15)

# set the project name and version
project(congFD VERSION 0.0.0)

configure_file(helloConfig.h.in helloConfig.h)
set(CMAKE_CXX_STANDARD 23)
#string(TIMESTAMP COMPILE_TIME %Y%m%d-%H%M%S)


#add_subdirectory(libMeshb)
#list(APPEND EXTRA_LIBS libMeshb)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -w -O0")
set(CMAKE_C_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -O0")
set(CMAKE_C_FLAGS_Release "${CMAKE_CXX_FLAGS_Release} -O0")

find_package(OpenMP)
#add the executable
add_executable(${PROJECT_NAME} zoneMain.cpp)

target_link_libraries(${PROJECT_NAME} PUBLIC libcongFD OpenMP::OpenMP_CXX)


add_executable(symmetrytest symmetrytest.cpp)
target_link_libraries(symmetrytest PUBLIC libcongFD OpenMP::OpenMP_CXX)



add_subdirectory(src)
add_subdirectory(utils)

include(CMakePrintHelpers)

cmake_print_properties(
TARGETS
OpenMP::OpenMP_CXX
PROPERTIES
INTERFACE_COMPILE_OPTIONS
INTERFACE_INCLUDE_DIRECTORIES
INTERFACE_LINK_LIBRARIES
)




