cmake_minimum_required(VERSION 3.20)
project(cspike)

set(CMAKE_CXX_STANDARD 11)

add_compile_options(-std=c++11 -O3 -m64)
aux_source_directory(. SRC_LIST)
add_executable(cspike_th_mem ${SRC_LIST})






