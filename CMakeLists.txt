cmake_minimum_required(VERSION 3.25)
project(segment_spmm)

set(CMAKE_CXX_STANDARD 20)

# Verbose CMake
set(CMAKE_VERBOSE_MAKEFILE ON)

# Matrix Market format support
add_library(matrix_market extern/matrix_market/src/mmio.c)
target_include_directories(matrix_market PUBLIC extern/matrix_market/include)

# Sparse matrix library
add_library(sparse_matrix lib/src/sparse_matrix.cpp)
target_include_directories(sparse_matrix PUBLIC lib/include)

# Dense matrix library
add_library(dense_matrix lib/src/dense_matrix.cpp)
target_include_directories(dense_matrix PUBLIC lib/include)

# IO library
add_library(mt_io lib/src/mt_io.cpp)
target_link_libraries(mt_io PUBLIC matrix_market sparse_matrix dense_matrix)

# SpMM library
add_library(spmm_csr src/spmm_csr.cpp)
target_include_directories(spmm_csr PUBLIC include)
target_link_libraries(spmm_csr PUBLIC sparse_matrix dense_matrix)

add_executable(segment_spmm main.cpp)

add_executable(run_dram_spmm_csr app/run_dram_spmm_csr.cpp)
target_link_libraries(run_dram_spmm_csr PRIVATE mt_io spmm_csr)
