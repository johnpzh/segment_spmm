//
// Created by Zhen Peng on 1/27/2024.
//

#ifndef SPMM_FORMAT_MT_IO_H
#define SPMM_FORMAT_MT_IO_H

#include <iostream>
#include <vector>
#include <string.h>
#include "mmio.h"
#include "sparse_matrix.h"
#include "dense_matrix.h"

namespace MT {

/// Read Matrix Market file
int read_mm(const char *filename,
            MM_typecode &matcode,
            int &M, /// number of rows
            int &N, /// number of cols
            int &nz,  /// number of non-zeros
            std::vector< std::vector<std::pair<int, double> > > &edgelist);

/// Read CSR Matrix
void read_csr_matrix(const char *filename, DRAMCSRMatrix *csr_matrix);

void print_csr_matrix(const DRAMCSRMatrix *csr_matrix);

/// Create Dense Matrix
void create_random_dense_matrix(uint64_t num_rows,
                                uint64_t num_cols,
                                DenseMatrix *dense_matrix);

}  /// namespace MT

#endif //SPMM_FORMAT_MT_IO_H
