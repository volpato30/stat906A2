cmake_minimum_required(VERSION 3.12)
project(Stat906A2)


set(CMAKE_CXX_STANDARD 14)
include_directories(./include)

add_executable(Stat906A2 main.cpp include/q1.h q1.cpp include/q2.h q2.cpp include/helper.h)
