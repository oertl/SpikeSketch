cmake_minimum_required(VERSION 3.21)
project(H3L_HLL)

set(CMAKE_CXX_STANDARD 11)

aux_source_directory(. SRC_LIST)
add_executable(H3L_HLL ${SRC_LIST})
add_compile_options(-std=c++11 -O3 -m64)
