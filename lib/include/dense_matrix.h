//
// Created by Zhen Peng on 1/26/2024.
//

#ifndef SPMM_FORMAT_DENSE_MATRIX_H
#define SPMM_FORMAT_DENSE_MATRIX_H

#include <iostream>

namespace MT {

struct DenseMatrix {
  int64_t num_rows_ = 0;
  int64_t num_cols_ = 0;
  double *values_ = nullptr;

  DenseMatrix() = default;
  DenseMatrix(int64_t num_rows,
              int64_t num_cols);
  virtual ~DenseMatrix();

  void alloc_values(int64_t size);
  void alloc_matrix(int64_t num_rows,
                    int64_t num_cols);

  void reset();

};

} /// namespace MT


#endif //SPMM_FORMAT_DENSE_MATRIX_H
