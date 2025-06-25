//
// Created by Zhen Peng on 6/24/2025.
//

#ifndef SPMM_FORMAT_SPMM_SCR_H
#define SPMM_FORMAT_SPMM_SCR_H

#include "sparse_matrix.h"
#include "dense_matrix.h"

namespace MT {

void spmm_csr_seq(const MT::CSRMatrix *A,
                  const MT::DenseMatrix *B,
                  MT::DenseMatrix *C/*out*/);

void spmm_csr_seq_column_segment(const MT::CSRMatrix *A, int64_t A1_tile, int64_t A2_tile,
                                 const MT::DenseMatrix *B, int64_t B1_tile, int64_t B2_tile,
                                 MT::DenseMatrix *C/*out*/);

}  /// namespace MT

#endif //SPMM_FORMAT_SPMM_SCR_H
