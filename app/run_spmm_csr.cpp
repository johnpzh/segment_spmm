//
// Created by Zhen Peng on 6/24/25.
//

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <chrono>
#include "mt_io.h"
#include "utils.h"
#include "spmm_csr.h"


int main(int argc, char *argv[])
{
  if (argc < 2) {
    fprintf(stderr, "Usage: %s <dir> [debug]\n", argv[0]);
    exit(EXIT_FAILURE);
  }

  bool debugging = false;
  char *filename = nullptr;
  if (strcmp(argv[argc - 1], "debug") == 0) {
    debugging = true;
    filename = argv[1];
  }

  /// Sparse matrix A
  MT::CSRMatrix csr_matrix;
  MT::read_csr_matrix(filename, &csr_matrix);

  /// Print matrix A
  if (debugging) {
    MT::print_csr_matrix(&csr_matrix);
  }

  /// Dense matrix B
  int64_t A1 = csr_matrix.num_rows_;
  int64_t A2 = csr_matrix.num_cols_;
  int64_t B2 = MAX_B2;
  if (debugging) {
//    B2 = A2;
    B2 = 4;
  }
  MT::DenseMatrix dense_matrix;
  MT::create_random_dense_matrix(A2, B2, &dense_matrix);

  /// Dense matrix C
  MT::DenseMatrix output_matrix(A1, B2);
  output_matrix.reset();

  /// Tile sizes
  int64_t tile_dim_size = 3;
  int64_t A1_tile = tile_dim_size;
  int64_t A2_tile = tile_dim_size;
  int64_t B1_tile = A2_tile;
  int64_t B2_tile = tile_dim_size;

  /// Kernel
  auto tt_start = std::chrono::high_resolution_clock::now();
//  MT::spmm_csr_seq(/*A=*/&csr_matrix, /*B=*/&dense_matrix, /*C=*/&output_matrix);
  MT::spmm_csr_seq_column_segment(
      /*A=*/&csr_matrix, A1_tile, A2_tile,
      /*B=*/&dense_matrix, B1_tile, B2_tile,
      /*C=*/&output_matrix);
  auto tt_end = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> tt_duration = tt_end - tt_start;
  std::cout << "time_exe(s): " << tt_duration.count() << "\n";


  if (debugging) {
    /// Print output
    printf("\nPrinting output ...\n");
    for (int64_t row = 0; row < A1; ++row) {
      for (int64_t col = 0; col < B2; ++col) {
        printf("%f ", output_matrix.values_[row * B2 + col]);
      }
      printf("\n");
    }
  }

  return EXIT_SUCCESS;
}