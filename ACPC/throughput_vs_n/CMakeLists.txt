cmake_minimum_required(VERSION 3.10)
project(Apache_CPC)

set(CMAKE_CXX_STANDARD 11)

aux_source_directory(. SRC_LIST)
add_executable(Apache_CPC_th_n ${SRC_LIST})
add_compile_options(-std=c++11 -O3 -m64)
