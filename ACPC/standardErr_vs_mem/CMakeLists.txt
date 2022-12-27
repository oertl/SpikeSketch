cmake_minimum_required(VERSION 3.21)
project(Apache_CPC)

set(CMAKE_CXX_STANDARD 11)

aux_source_directory(. SRC_LIST)
add_executable(ACPC_stdError_mem ${SRC_LIST})
add_compile_options(-std=c++11 -O3 -m64)
