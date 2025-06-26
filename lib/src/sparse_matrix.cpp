//
// Created by Zhen Peng on 1/26/2024.
//

#include <iostream>
#include "sparse_matrix.h"


/******* DRAMCSRMatrix **********/

void MT::DRAMCSRMatrix::alloc_row_offsets_(uint64_t size) {
  free(row_offsets_);
  row_offsets_ = (uint64_t *) malloc(size * sizeof(*row_offsets_));
  if (row_offsets_ == nullptr) {
    fprintf(stderr, "Error: %s:%d MT::DRAMCSRMatrix::alloc_row_b_offsets_() failed.\n", __FILE__, __LINE__);
    return;
  }
}

void MT::DRAMCSRMatrix::alloc_col_ind(uint64_t size) {
  free(col_ind_);
  col_ind_ = (uint64_t *) malloc(size * sizeof(*col_ind_));
  if (col_ind_ == nullptr) {
    fprintf(stderr, "Error: %s:%d MT::DRAMCSRMatrix::alloc_col_b_ind() failed.\n", __FILE__, __LINE__);
    return;
  }
}

void MT::DRAMCSRMatrix::alloc_values(uint64_t size) {
  free(values_);
  values_ = (double *) malloc(size * sizeof(*values_));
  if (values_ == nullptr) {
    fprintf(stderr, "Error: %s:%d MT::DRAMCSRMatrix::alloc_values() failed.\n", __FILE__, __LINE__);
    return;
  }
}

void MT::DRAMCSRMatrix::alloc_matrix(uint64_t num_rows, uint64_t num_cols, uint64_t nnz) {
  num_rows_ = num_rows;
  num_cols_ = num_cols;
  nnz_ = nnz;
  alloc_row_offsets_(num_rows + 1);
  alloc_col_ind(nnz);
  alloc_values(nnz);
}

MT::DRAMCSRMatrix::DRAMCSRMatrix(uint64_t num_rows, uint64_t num_cols, uint64_t nnz) {
  alloc_matrix(num_rows, num_cols, nnz);
}

MT::DRAMCSRMatrix::~DRAMCSRMatrix() {
  free(row_offsets_);
  free(col_ind_);
  free(values_);
}

void MT::get_tiling_table_for_pos(
    const DRAMCSRMatrix *A,
    uint64_t row_idx_start,
    uint64_t row_idx_bound,
    uint64_t A2_tile,
    uint64_t &num_tiles/*out*/,
    std::vector<std::vector<uint64_t>> &table/*out*/)
{
  uint64_t left_col_id = UINT64_MAX;
  uint64_t right_col_id = 0;
  for (uint64_t i = row_idx_start; i < row_idx_bound; ++i) {
    for (uint64_t k_loc = A->row_offsets_[i]; k_loc < A->row_offsets_[i + 1]; ++k_loc) {
      uint64_t k = A->col_ind_[k_loc];
      if (k < left_col_id) {
        left_col_id = k;
      }
      if (k > right_col_id) {
        right_col_id = k;
      }
    }
  }
  num_tiles = ((right_col_id - left_col_id + 1) + (A2_tile - 1)) / A2_tile;
  table = std::vector<std::vector<uint64_t>>(row_idx_bound - row_idx_start, std::vector<uint64_t>(num_tiles + 1));

  for (uint64_t row_idx = row_idx_start; row_idx < row_idx_bound; ++row_idx) {
    uint64_t k_loc_start = A->row_offsets_[row_idx];
    uint64_t k_loc_bound = A->row_offsets_[row_idx + 1];
    uint64_t curr_k_loc = k_loc_start;
    table[row_idx - row_idx_start][0] = k_loc_start;

    uint64_t col_id_start = left_col_id;
    for (uint64_t tile_idx = 1; tile_idx < num_tiles + 1; ++tile_idx) {
      uint64_t col_id_bound = std::min(col_id_start + A2_tile, right_col_id + 1);
      while (curr_k_loc < k_loc_bound) {
        uint64_t k = A->col_ind_[curr_k_loc];
        if (k >= col_id_bound || k < col_id_start) {
          break;
        }
        ++curr_k_loc;
      }
      table[row_idx - row_idx_start][tile_idx] = curr_k_loc;
      col_id_start += A2_tile;
    }
  }

//  /// test
//  std::cout << __FILE__ << ":" << __LINE__ << "\n";
//  std::cout << "row_idx_start: " << row_idx_start << "\n";
//  std::cout << "row_idx_bound: " << row_idx_bound << "\n";
//  std::cout << "A2_tile: " << A2_tile << "\n";
//  std::cout << "num_tiles: " << num_tiles << "\n";
//  printf("print the table...\n");
//  for (auto &row : table) {
//    for (auto e : row) {
//      std::cout << e << " ";
//    }
//    std::cout << "\n";
//  }
//  std::cout << "\n";
//  /// end test
}

/******* End DRAMCSRMatrix **********/