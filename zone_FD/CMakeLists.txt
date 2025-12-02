# 指定构建此项目所需的最低CMake版本为3.15，必须放在文件开头[1,6,7](@ref)
# 如果当前CMake版本低于此要求，配置过程将报错停止[6](@ref)
cmake_minimum_required(VERSION 3.15)

# 设置项目名称为congFD，版本号为0.0.0[1,5](@ref)
# project命令会定义一些有用的变量，如PROJECT_NAME(值为congFD)和PROJECT_VERSION(值为0.0.0)[5](@ref)
project(congFD VERSION 0.0.0)

# 配置头文件：将模板文件helloConfig.h.in中的CMake变量替换为实际值，生成helloConfig.h[1](@ref)
# 这样可以在C++代码中使用项目配置信息（如版本号）
configure_file(helloConfig.h.in helloConfig.h)

# 设置C++编译标准为C++23[1](@ref)
set(CMAKE_CXX_STANDARD 23)

# 这行被注释掉的代码原本用于生成时间戳字符串并存储在COMPILE_TIME变量中[1](@ref)
# 时间戳格式为年月日-时分秒（如20251201-143055）
#string(TIMESTAMP COMPILE_TIME %Y%m%d-%H%M%S)

# 以下两行被注释的代码涉及添加子目录和库链接[1](@ref)
#add_subdirectory(libMeshb)   # 添加libMeshb子目录，CMake会处理其中的CMakeLists.txt
#list(APPEND EXTRA_LIBS libMeshb)  # 将libMeshb库添加到EXTRA_LIBS变量中

# 这行被注释的代码用于设置编译模式为Release（发布模式）[1](@ref)
# set(CMAKE_BUILD_TYPE Release)  # 这一行段代码，设置编译模式为Release

# 以下三行被注释的代码设置了调试版本的编译选项（无优化-O0）[1](@ref)
# set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -w -O0")
# set(CMAKE_C_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -O0")
# set(CMAKE_C_FLAGS_Release "${CMAKE_CXX_FLAGS_Release} -O0")

# 设置C++编译器标志：添加-g（调试信息）、-w（抑制所有警告）、-O2（二级优化）[1](@ref)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -w -O2")

# 设置Debug模式的C编译器标志：添加-O2优化级别[1](@ref)
# set(CMAKE_C_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -O2")
set(CMAKE_C_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -O0")

# 设置Release模式的C编译器标志：添加-O2优化级别[1](@ref)
set(CMAKE_C_FLAGS_Release "${CMAKE_CXX_FLAGS_Release} -O2")

# 查找OpenMP包，OpenMP用于支持并行编程[1](@ref)
# 如果找到OpenMP，会定义OpenMP::OpenMP_CXX等目标供后续链接使用
find_package(OpenMP)

# 创建名为${PROJECT_NAME}（即congFD）的可执行文件，源文件为zoneMain.cpp[1,4](@ref)
add_executable(${PROJECT_NAME} zoneMain.cpp)

# 为congFD目标链接库：libcongFD和OpenMP::OpenMP_CXX[1](@ref)
# PUBLIC关键字表示这些库的依赖关系会传递给链接congFD的其他目标
target_link_libraries(${PROJECT_NAME} PUBLIC libcongFD OpenMP::OpenMP_CXX)

# 创建另一个名为symmetrytest的可执行文件，源文件为symmetrytest.cpp[1,4](@ref)
add_executable(symmetrytest symmetrytest.cpp)

# 为symmetrytest目标链接相同的库[1](@ref)
target_link_libraries(symmetrytest PUBLIC libcongFD OpenMP::OpenMP_CXX)

# 添加src和utils子目录，CMake会处理这些子目录中的CMakeLists.txt文件[1](@ref)
# 通常这些子目录中定义了libcongFD库的构建规则
add_subdirectory(src)
add_subdirectory(utils)

# 包含CMakePrintHelpers模块，该模块提供了cmake_print_properties等调试辅助命令[3](@ref)
include(CMakePrintHelpers)

# 打印OpenMP::OpenMP_CXX目标的属性，用于调试和验证OpenMP配置[3](@ref)
# 会显示编译选项、头文件目录和链接库信息
cmake_print_properties(
TARGETS
OpenMP::OpenMP_CXX
PROPERTIES
INTERFACE_COMPILE_OPTIONS      # 接口编译选项
INTERFACE_INCLUDE_DIRECTORIES  # 接口头文件目录
INTERFACE_LINK_LIBRARIES       # 接口链接库
)
