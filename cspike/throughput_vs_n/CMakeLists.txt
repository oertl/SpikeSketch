cmake_minimum_required(VERSION 3.10)
project(cspike)

set(CMAKE_CXX_STANDARD 11)

add_compile_options(-std=c++11 -O3 -m64)
aux_source_directory(. SRC_LIST)
add_executable(cspike_th_n ${SRC_LIST})






