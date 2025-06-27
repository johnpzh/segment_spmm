//
// Created by Zhen Peng on 1/26/2024.
//

#ifndef SPMM_FORMAT_DENSE_MATRIX_H
#define SPMM_FORMAT_DENSE_MATRIX_H

#include <iostream>

namespace MT {

struct DRAMDenseMatrix {
  uint64_t num_rows_ = 0;
  uint64_t num_cols_ = 0;
  double *values_ = nullptr;

  DRAMDenseMatrix() = default;
  DRAMDenseMatrix(uint64_t num_rows,
                  uint64_t num_cols);
  virtual ~DRAMDenseMatrix();

  void alloc_values(uint64_t size);
  void alloc_matrix(uint64_t num_rows,
                    uint64_t num_cols);

  void reset();

};

} /// namespace MT


#endif //SPMM_FORMAT_DENSE_MATRIX_H
