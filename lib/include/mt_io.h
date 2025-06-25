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

/// Read COO Matrix
void read_coo_matrix(const char *filename, COOMatrix *coo_matrix);

void print_coo_matrix(const COOMatrix *coo_matrix);

/// Read CSR Matrix
void read_csr_matrix(const char *filename, CSRMatrix *csr_matrix);

void print_csr_matrix(const CSRMatrix *csr_matrix);

/// Read BSR Matrix
void read_bsr_matrix(const char *filename,
                     int64_t row_b_size,  /// number of rows in each block
                     int64_t col_b_size,  /// number of columns in each block
                     BSRMatrix *bsr_matrix  /* output */);

void print_bsr_matrix(const BSRMatrix *matrix);

/// Read Blocked Ellpack Matrix
void read_blocked_ell_matrix(const char *filename,
                             int64_t bsize,
                             BlockedEllMatrix *blocked_ell_matrix);

void print_blocked_ell_matrix(const BlockedEllMatrix *blocked_ell_matrix);

/// Read Sliced Ellpack Matrix
void read_sliced_ell_matrix(const char *filename,
                            int64_t slice_size,
                            SlicedEllMatrix *sliced_ell_matrix);

void print_sliced_ell_matrix(const SlicedEllMatrix *sliced_ell_matrix);

/// Create Dense Matrix
void create_random_dense_matrix(int64_t num_rows,
                                int64_t num_cols,
                                DenseMatrix *dense_matrix);

}  /// namespace MT

#endif //SPMM_FORMAT_MT_IO_H
