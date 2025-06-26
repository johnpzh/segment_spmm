//
// Created by Zhen Peng on 1/26/2024.
//

#ifndef SPMM_FORMAT_SPARSE_MATRIX_H
#define SPMM_FORMAT_SPARSE_MATRIX_H

#include <iostream>
#include <vector>
#include <limits.h>

namespace MT {

/// CSR Matrix
struct DRAMCSRMatrix {
  uint64_t num_rows_ = 0;
  uint64_t num_cols_ = 0;
  uint64_t nnz_ = 0;
  uint64_t *row_offsets_ = nullptr;
  uint64_t *col_ind_ = nullptr;
  double *values_ = nullptr;

  DRAMCSRMatrix() = default;
  DRAMCSRMatrix(uint64_t num_rows,
                uint64_t num_cols,
                uint64_t nnz);
  virtual ~DRAMCSRMatrix();

  void alloc_row_offsets_(uint64_t size);
  void alloc_col_ind(uint64_t size);
  void alloc_values(uint64_t size);
  void alloc_matrix(uint64_t num_rows,
                    uint64_t num_cols,
                    uint64_t nnz);
};

void get_tiling_table_for_pos(
    const DRAMCSRMatrix *A,
    uint64_t row_idx_start,
    uint64_t row_idx_bound,
    uint64_t A2_tile,
    uint64_t &num_tiles/*out*/,
    std::vector<std::vector<uint64_t>> &table/*out*/);

} /// namespace MT

#endif //SPMM_FORMAT_SPARSE_MATRIX_H
