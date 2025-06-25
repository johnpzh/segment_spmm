//
// Created by Zhen Peng on 1/26/2024.
//

#include <string.h>
#include "dense_matrix.h"

void MT::DenseMatrix::alloc_values(uint64_t size) {
  free(values_);
  values_ = (double *) malloc(size * sizeof(double));
  if (values_ == nullptr) {
    fprintf(stderr, "Error: %s:%d MT::DenseMatrix::alloc_values() failed.\n", __FILE__, __LINE__);
    return;
  }
}

void MT::DenseMatrix::alloc_matrix(uint64_t num_rows, uint64_t num_cols) {
  num_rows_ = num_rows;
  num_cols_ = num_cols;
  alloc_values(num_rows * num_cols);
}

MT::DenseMatrix::DenseMatrix(uint64_t num_rows, uint64_t num_cols) {
  alloc_matrix(num_rows, num_cols);
}

void MT::DenseMatrix::reset() {
  memset(values_, 0, num_rows_ * num_cols_ * sizeof(*values_));
}

MT::DenseMatrix::~DenseMatrix() {
  free(values_);
}