//
// Created by Zhen Peng on 6/24/2025.
//

#ifndef SPMM_FORMAT_SPMM_SCR_H
#define SPMM_FORMAT_SPMM_SCR_H

#include "sparse_matrix.h"
#include "dense_matrix.h"

namespace MT {

void spmm_csr_seq(const MT::DRAMCSRMatrix *A,
                  const MT::DRAMDenseMatrix *B,
                  MT::DRAMDenseMatrix *C/*out*/);

void spmm_csr_seq_column_segment_v0(const MT::DRAMCSRMatrix *A, uint64_t A1_tile, uint64_t A2_tile,
                                    const MT::DRAMDenseMatrix *B, uint64_t B1_tile, uint64_t B2_tile,
                                    MT::DRAMDenseMatrix *C/*out*/);
void spmm_csr_seq_column_segment_v1(const MT::DRAMCSRMatrix *A, uint64_t A1_tile, uint64_t A2_tile,
                                    const MT::DRAMDenseMatrix *B, uint64_t B1_tile, uint64_t B2_tile,
                                    MT::DRAMDenseMatrix *C/*out*/);
void spmm_csr_seq_column_segment_v2(const MT::DRAMCSRMatrix *A, uint64_t A1_tile, uint64_t A2_tile,
                                    const MT::DRAMDenseMatrix *B, uint64_t B1_tile, uint64_t B2_tile,
                                    MT::DRAMDenseMatrix *C/*out*/);

}  /// namespace MT

#endif //SPMM_FORMAT_SPMM_SCR_H
