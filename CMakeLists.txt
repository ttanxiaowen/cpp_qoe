cmake_minimum_required(VERSION 3.30)
project(c_new_qoe_1)

set(CMAKE_CXX_STANDARD 20)

add_executable(c_new_qoe
        calculate.cpp
        calculate.h
        resource.h
        dijkstra.h
        dijkstra.cpp
        q_max_qoe.cpp
        main.cpp)
