cmake_minimum_required(VERSION 3.21)
project(HLLL)

set(CMAKE_CXX_STANDARD 11)

aux_source_directory(. SRC_LIST)
add_executable(HLLL ${SRC_LIST})
add_compile_options(-std=c++11 -O3 -m64)
