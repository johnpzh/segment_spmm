//
// Created by Zhen Peng on 6/24/2025.
//

#include "spmm_csr.h"
//#include <omp.h>

/// C[i,j] = A[i,k] * B[k,j]
void MT::spmm_csr_seq(const MT::DRAMCSRMatrix *A,
                      const MT::DenseMatrix *B,
                      MT::DenseMatrix *C/*out*/)
{
  uint64_t A1 = A->num_rows_;
  uint64_t B2 = B->num_cols_;

  for (uint64_t i = 0; i < A1; ++i) {
    for (uint64_t k_loc = A->row_offsets_[i]; k_loc < A->row_offsets_[i + 1]; ++k_loc) {
      uint64_t k = A->col_ind_[k_loc];
      double A_val = A->values_[k_loc];

      for (uint64_t j = 0; j < B2; ++j) {
        double B_val = B->values_[k * B2 + j];
        C->values_[i * B2 + j] += A_val * B_val;
      }
    }
  }
}

/// C[i,j] = A[i,k] * B[k,j]
void MT::spmm_csr_seq_column_segment_v0(const MT::DRAMCSRMatrix *A, uint64_t A1_tile, uint64_t A2_tile,
                                     const MT::DenseMatrix *B, uint64_t B1_tile, uint64_t B2_tile,
                                     MT::DenseMatrix *C/*out*/)
{
  uint64_t A1 = A->num_rows_;
  uint64_t A2 = A->num_cols_;
  uint64_t B2 = B->num_cols_;

  for (uint64_t ii = 0; ii < A1; ii += A1_tile) {
    uint64_t i_bound = std::min(ii + A1_tile, A1);
    for (uint64_t kk = 0; kk < A2; kk += A2_tile) {
      uint64_t k_bound = std::min(kk + A2_tile, A2);
      for (uint64_t jj = 00; jj < B2; jj += B2_tile) {
        uint64_t j_bound = std::min(jj + B2_tile, B2);
        /// Tile
        for (uint64_t i = ii; i < i_bound; ++i) {
          for (uint64_t k_loc = A->row_offsets_[i]; k_loc < A->row_offsets_[i + 1]; ++k_loc) {
            uint64_t k = A->col_ind_[k_loc];
            if (k < kk) {
              continue;
            }
            if (k >= k_bound) {
              break;
            }
            double A_val = A->values_[k_loc];
            for (uint64_t j = jj; j < j_bound; ++j) {
              double B_val = B->values_[k * B2 + j];
              C->values_[i * B2 + j] += A_val * B_val;
            }
          }
        }
      }
    }
  }
}

/// C[i,j] = A[i,k] * B[k,j]
void MT::spmm_csr_seq_column_segment_v1(const MT::DRAMCSRMatrix *A, uint64_t A1_tile, uint64_t A2_tile,
                                     const MT::DenseMatrix *B, uint64_t B1_tile, uint64_t B2_tile,
                                     MT::DenseMatrix *C/*out*/)
{
  uint64_t A1 = A->num_rows_;
  uint64_t A2 = A->num_cols_;
  uint64_t B2 = B->num_cols_;

  for (uint64_t ii = 0; ii < A1; ii += A1_tile) {
    uint64_t i_bound = std::min(ii + A1_tile, A1);
    for (uint64_t i = ii; i < i_bound; ++i) {
      for (uint64_t kk = 0; kk < A2; kk += A2_tile) {
        uint64_t k_bound = std::min(kk + A2_tile, A2);
        for (uint64_t k_loc = A->row_offsets_[i]; k_loc < A->row_offsets_[i + 1]; ++k_loc) {
          uint64_t k = A->col_ind_[k_loc];
          if (k < kk) {
            continue;
          }
          if (k >= k_bound) {
            break;
          }
          double A_val = A->values_[k_loc];
          for (uint64_t jj = 0; jj < B2; jj += B2_tile) {
            uint64_t j_bound = std::min(jj + B2_tile, B2);
            for (uint64_t j = jj; j < j_bound; ++j) {
              double B_val = B->values_[k * B2 + j];
              C->values_[i * B2 + j] += A_val * B_val;
            }
          }
        }
      }
    }
  }
}

/// C[i,j] = A[i,k] * B[k,j]
void MT::spmm_csr_seq_column_segment_v2(const MT::DRAMCSRMatrix *A, uint64_t A1_tile, uint64_t A2_tile,
                                     const MT::DenseMatrix *B, uint64_t B1_tile, uint64_t B2_tile,
                                     MT::DenseMatrix *C/*out*/)
{
  uint64_t A1 = A->num_rows_;
  uint64_t A2 = A->num_cols_;
  uint64_t B2 = B->num_cols_;

  for (uint64_t ii = 0; ii < A1; ii += A1_tile) {
    uint64_t i_bound = std::min(ii + A1_tile, A1);
    uint64_t num_tiles;
    std::vector<std::vector<uint64_t>> table;
    MT::get_tiling_table_for_pos(A,
                                 /*row_idx_start=*/ii,
                                 /*row_idx_bound=*/i_bound,
                                 /*A2_tile=*/A2_tile,
                                 num_tiles/*out*/,
                                 table/*out*/);
    for (uint64_t kk_tile_i = 0; kk_tile_i < num_tiles; ++kk_tile_i) {
      for (uint64_t jj = 0; jj < B2; jj += B2_tile) {
        uint64_t j_bound = std::min(jj + B2_tile, B2);
        /// Tile
        for (uint64_t i = ii; i < i_bound; ++i) {
//          /// test
//          printf("table[%llu][%llu]: %llu\n", i - ii, kk_tile_i, table[i - ii][kk_tile_i]);
//          printf("-> table[%llu][%llu]: %llu\n", i - ii, kk_tile_i + 1, table[i - ii][kk_tile_i + 1]);
//          /// end
          for (uint64_t k_loc = table[i - ii][kk_tile_i]; k_loc < table[i - ii][kk_tile_i + 1]; ++k_loc) {
            uint64_t k = A->col_ind_[k_loc];
            double A_val = A->values_[k_loc];
            for (uint64_t j = jj; j < j_bound; ++j) {
              C->values_[i * B2 + j] += A_val * B->values_[k * B2 + j];
            }
          }
        }
      }
    }
  }
}