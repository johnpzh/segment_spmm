//
// Created by Zhen Peng on 6/24/25.
//

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <chrono>
#include <vector>
#include <string>
#include "mt_io.h"
#include "utils.h"
#include "spmm_csr.h"

std::vector<std::string> MATRICES = {
    "bcsstk17",
    "cop20k_A",
    "scircuit"
};

int main(int argc, char *argv[])
{
  if (argc < 2) {
    fprintf(stderr, "Usage: %s <dir>\n", argv[0]);
    exit(EXIT_FAILURE);
  }

  auto master_tt_start = std::chrono::high_resolution_clock::now();

//  bool debugging = false;
//  char *filename = nullptr;
//  if (strcmp(argv[argc - 1], "debug") == 0) {
//    debugging = true;
//    filename = argv[1];
//  }

  uint64_t num_repeats = 4;

  std::vector<std::string> matrix_list;
  std::vector<double> spmm_no_tiling_avg_time_list;
  std::vector<double> spmm_tiling_avg_time_list;

  for (std::string &mtx : MATRICES) {
    std::string file_name = std::string{argv[1]} + "/" + mtx + "/" + mtx + ".mtx";

    matrix_list.push_back(mtx);
    std::vector<double> no_tiling_exe_time_list;
    std::vector<double> tiling_exe_time_list;

    /// No Tiling
    for (uint64_t r_i = 0; r_i < num_repeats; ++r_i) {
      /// Sparse matrix A
      MT::DRAMCSRMatrix csr_matrix;
      MT::read_csr_matrix(file_name.c_str(), &csr_matrix);

      /// Dense matrix B
      int64_t A1 = csr_matrix.num_rows_;
      int64_t A2 = csr_matrix.num_cols_;
      int64_t B2 = MAX_B2;
      MT::DenseMatrix dense_matrix;
      MT::create_random_dense_matrix(A2, B2, &dense_matrix);

      /// Dense matrix C
      MT::DenseMatrix output_matrix(A1, B2);
      output_matrix.reset();

      /// Kernel
      auto tt_start = std::chrono::high_resolution_clock::now();
      MT::spmm_csr_seq(/*A=*/&csr_matrix, /*B=*/&dense_matrix, /*C=*/&output_matrix);
      auto tt_end = std::chrono::high_resolution_clock::now();
      std::chrono::duration<double> tt_duration = tt_end - tt_start;
      std::cout << "DRAM no-tiling, matrix: " << mtx << ", time_exe(s): " << tt_duration.count() << "\n";

      no_tiling_exe_time_list.push_back(tt_duration.count());
    }
    /// Save no-tiling all exe time
    std::string no_tiling_filename;
    {
      std::stringstream ss;
      ss << "output.spmm.dram.no-tiling.mtx-" << mtx << ".log";
      no_tiling_filename = ss.str();
    }
    save_exe_times_into_file(no_tiling_filename, no_tiling_exe_time_list);
    /// Calculate no-tiling average
    double no_tiling_avg_time = calculate_average_in_vector(no_tiling_exe_time_list);
    spmm_no_tiling_avg_time_list.push_back(no_tiling_avg_time);

    /// Tiling
    for (uint64_t r_i = 0; r_i < num_repeats; ++r_i) {
      /// Sparse matrix A
      MT::DRAMCSRMatrix csr_matrix;
      MT::read_csr_matrix(file_name.c_str(), &csr_matrix);

      /// Dense matrix B
      int64_t A1 = csr_matrix.num_rows_;
      int64_t A2 = csr_matrix.num_cols_;
      int64_t B2 = MAX_B2;
      MT::DenseMatrix dense_matrix;
      MT::create_random_dense_matrix(A2, B2, &dense_matrix);

      /// Dense matrix C
      MT::DenseMatrix output_matrix(A1, B2);
      output_matrix.reset();

      /// Tile sizes
      int64_t tile_dim_size = 512;
      int64_t A1_tile = tile_dim_size;
      int64_t A2_tile = tile_dim_size;
      int64_t B1_tile = A2_tile;
      int64_t B2_tile = tile_dim_size;

      /// Kernel
      auto tt_start = std::chrono::high_resolution_clock::now();
      MT::spmm_csr_seq_column_segment_v1(
          /*A=*/&csr_matrix, A1_tile, A2_tile,
          /*B=*/&dense_matrix, B1_tile, B2_tile,
          /*C=*/&output_matrix);
      auto tt_end = std::chrono::high_resolution_clock::now();
      std::chrono::duration<double> tt_duration = tt_end - tt_start;
      std::cout << "DRAM tiling, matrix: " << mtx << ", time_exe(s): " << tt_duration.count() << "\n";

      tiling_exe_time_list.push_back(tt_duration.count());
    }
    /// Save all exe times
    std::string tiling_filename;
    {
      std::stringstream ss;
      ss << "output.spmm.dram.tiling.mtx-" << mtx << ".log";
      tiling_filename = ss.str();
    }
    save_exe_times_into_file(tiling_filename, tiling_exe_time_list);

    /// Calculate average
    double tiling_avg_time = calculate_average_in_vector(tiling_exe_time_list);
    spmm_tiling_avg_time_list.push_back(tiling_avg_time);
  }

  /// Save results to a collection file
  std::string collect_filename("output.spmm.dram.collection.csv");
  std::ofstream fout;
  fout.open(collect_filename);
  if (fout.is_open()) {
    std::string header("Matrix,DRAM.No-Tiling(s),DRAM.Tiling(s)");
    fout << header << "\n";
    for (uint64_t row_i = 0; row_i < matrix_list.size(); ++row_i) {
      fout << matrix_list[row_i] << "," << spmm_no_tiling_avg_time_list[row_i] << "," << spmm_tiling_avg_time_list[row_i] << "\n";
    }

    std::cout << "Saved to file " << collect_filename << "\n";
  } else {
    std::cerr << "Error: cannot open file " << collect_filename << "\n";
  }

  auto master_tt_end = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> master_tt_duration = master_tt_end - master_tt_start;
  std::cout << "TOTAL_EXE_TIME(S): " << master_tt_duration.count() << "\n";

  return EXIT_SUCCESS;
}


//int main(int argc, char *argv[])
//{
//  if (argc < 2) {
//    fprintf(stderr, "Usage: %s <dir> [debug]\n", argv[0]);
//    exit(EXIT_FAILURE);
//  }
//
//  bool debugging = false;
//  char *filename = nullptr;
//  if (strcmp(argv[argc - 1], "debug") == 0) {
//    debugging = true;
//    filename = argv[1];
//  }
//
//  /// Sparse matrix A
//  MT::DRAMCSRMatrix csr_matrix;
//  MT::read_csr_matrix(filename, &csr_matrix);
//
//  /// Print matrix A
//  if (debugging) {
//    MT::print_csr_matrix(&csr_matrix);
//  }
//
//  /// Dense matrix B
//  int64_t A1 = csr_matrix.num_rows_;
//  int64_t A2 = csr_matrix.num_cols_;
//  int64_t B2 = MAX_B2;
//  if (debugging) {
////    B2 = A2;
//    B2 = 4;
//  }
//  MT::DenseMatrix dense_matrix;
//  MT::create_random_dense_matrix(A2, B2, &dense_matrix);
//
//  /// Dense matrix C
//  MT::DenseMatrix output_matrix(A1, B2);
//  output_matrix.reset();
//
//  /// Tile sizes
//  int64_t tile_dim_size = 3;
//  int64_t A1_tile = tile_dim_size;
//  int64_t A2_tile = tile_dim_size;
//  int64_t B1_tile = A2_tile;
//  int64_t B2_tile = tile_dim_size;
//
//  /// Kernel
//  auto tt_start = std::chrono::high_resolution_clock::now();
////  MT::spmm_csr_seq(/*A=*/&csr_matrix, /*B=*/&dense_matrix, /*C=*/&output_matrix);
//  MT::spmm_csr_seq_column_segment(
//      /*A=*/&csr_matrix, A1_tile, A2_tile,
//      /*B=*/&dense_matrix, B1_tile, B2_tile,
//      /*C=*/&output_matrix);
//  auto tt_end = std::chrono::high_resolution_clock::now();
//  std::chrono::duration<double> tt_duration = tt_end - tt_start;
//  std::cout << "time_exe(s): " << tt_duration.count() << "\n";
//
//
//  if (debugging) {
//    /// Print output
//    printf("\nPrinting output ...\n");
//    for (int64_t row = 0; row < A1; ++row) {
//      for (int64_t col = 0; col < B2; ++col) {
//        printf("%f ", output_matrix.values_[row * B2 + col]);
//      }
//      printf("\n");
//    }
//  }
//
//  return EXIT_SUCCESS;
//}