# 设置CMake最低版本要求
cmake_minimum_required(VERSION 3.20)

# 定义项目名称
project(CPP_QOE_MAX)

# 设置C++标准为C++2a，因为GCC 9.4.0使用这个标志来支持C++20特性
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF) # 确保不使用编译器扩展

# 添加项目的源文件和头文件
add_executable(CPP_QOE_MAX
        min.h
        min.cpp
        max.h
        qoe.h
        max.cpp
        main.cpp
        calculate.cpp
        calculate.h
        resource.h
        dijkstra.h
        dijkstra.cpp
        qoe.cpp)

# 链接pthread库，因为您的项目使用了std::thread
find_package(Threads REQUIRED)
target_link_libraries(CPP_QOE_MAX PRIVATE Threads::Threads)