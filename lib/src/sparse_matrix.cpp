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

/******* End DRAMCSRMatrix **********/