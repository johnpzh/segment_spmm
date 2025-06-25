//
// Created by Zhen Peng on 1/26/2024.
//

#ifndef SPMM_FORMAT_SPARSE_MATRIX_H
#define SPMM_FORMAT_SPARSE_MATRIX_H

#include <iostream>

namespace MT {

/// COO Matrix
struct COOMatrix {
  int64_t num_rows_ = 0;
  int64_t num_cols_ = 0;
  int64_t nnz_ = 0;
  int64_t *row_ind_ = nullptr;
  int64_t *col_ind_ = nullptr;
  double *values_ = nullptr;

  COOMatrix() = default;
  COOMatrix(int64_t num_rows,
            int64_t num_cols,
            int64_t nnz);
  virtual ~COOMatrix();

  void alloc_row_ind(int64_t size);
  void alloc_col_ind(int64_t size);
  void alloc_values(int64_t size);
  void alloc_matrix(int64_t num_rows,
                    int64_t num_cols,
                    int64_t nnz);
};


/// CSR Matrix
struct CSRMatrix {
  int64_t num_rows_ = 0;
  int64_t num_cols_ = 0;
  int64_t nnz_ = 0;
  int64_t *row_offsets_ = nullptr;
  int64_t *col_ind_ = nullptr;
  double *values_ = nullptr;

  CSRMatrix() = default;
  CSRMatrix(int64_t num_rows,
            int64_t num_cols,
            int64_t nnz);
  virtual ~CSRMatrix();

  void alloc_row_offsets_(int64_t size);
  void alloc_col_ind(int64_t size);
  void alloc_values(int64_t size);
  void alloc_matrix(int64_t num_rows,
                    int64_t num_cols,
                    int64_t nnz);
};


/// BSR Matrix
struct BSRMatrix {
  int64_t num_b_rows_ = 0;  /// number of block rows
  int64_t num_b_cols_ = 0;  /// number of block columns
  int64_t bnnz_ = 0;  /// number of blocks
  int64_t row_b_size_ = 0;  /// number of rows of each block
  int64_t col_b_size_ = 0;  /// number of columns of each block
  int64_t nnz_ = 0;  /// number of non-zeros
  int64_t *row_b_offsets_ = nullptr;  /// Block row offsets. Array of size `num_b_rows_ + 1`.
  int64_t *col_b_ind_ = nullptr;  /// Block column indices. Array of size `bnnz_`.
  double *values_ = nullptr;  /// Values. Array of size `bnnz_ * row_b_size_ * col_b_size_`.

  BSRMatrix() = default;
  BSRMatrix(int64_t num_brows,
            int64_t num_bcols,
            int64_t bnnz,
            int64_t row_bsize,
            int64_t col_bsize,
            int64_t nnz);
  virtual ~BSRMatrix();

  void alloc_row_b_offsets_(int64_t size);
  void alloc_col_b_ind(int64_t size);
  void alloc_values(int64_t size);
  void alloc_matrix(int64_t num_brows,
                    int64_t num_bcols,
                    int64_t bnnz,
                    int64_t row_bsize,
                    int64_t col_bsize,
                    int64_t nnz);
};


/// Blocked-Ellpack Matrix
struct BlockedEllMatrix {
  int64_t num_rows_ = 0;  /// number of rows
  int64_t num_cols_ = 0;  /// number of columns
  int64_t bsize_ = 0;  /// number of rows or columns of each Ell-block (that is a square)
  int64_t ell_cols_ = 0;  /// Actual number of columns of the Blocked-Ellpack format (for the whole matrix)
  int64_t nnz_ = 0;  /// number of non-zeros
  int64_t *col_b_ind_ = nullptr;  /// Blocked-ELL Column indices. Array with `[num_rows_ / bsize_][ell_cols_ / bsize_]` elements
  double *values_ = nullptr;  /// Values. Array with `num_rows_ * ell_cols_` elements

  BlockedEllMatrix() = default;
  BlockedEllMatrix(int64_t num_rows,
                   int64_t num_cols,
                   int64_t bsize,
                   int64_t ell_cols,
                   int64_t nnz);
  virtual ~BlockedEllMatrix();

  void alloc_col_b_ind(int64_t size);
  void alloc_values(int64_t size);
  void alloc_matrix(int64_t num_rows,
                    int64_t num_cols,
                    int64_t bsize,
                    int64_t ell_cols,
                    int64_t nnz);
};


/// Sliced-Ellpack Matrix
struct SlicedEllMatrix {
  int64_t num_rows_ = 0;  /// number of rows
  int64_t num_cols_ = 0;  /// number of columns
  int64_t nnz_ = 0;  /// number of nonzeros
  int64_t values_size_ = 0;  /// total number of elements (nonzero and padding)
  int64_t slice_size_ = 0;  /// number of rows per slice
  int64_t *slice_offsets_ = nullptr;  /// slice offsets. Array of size `ceiling(num_rows_ / slice_size_) + 1`.
  int64_t *col_ind_ = nullptr;  /// column indices. Array of size `values_size_`.
  double *values_ = nullptr;  /// values. Array of size `values_size_`.

  SlicedEllMatrix() = default;
  SlicedEllMatrix(int64_t num_rows,
                  int64_t num_cols,
                  int64_t nnz,
                  int64_t values_size,
                  int64_t slice_size);
  virtual ~SlicedEllMatrix();

  void alloc_slice_offsets(int64_t size);
  void alloc_col_ind(int64_t size);
  void alloc_values(int64_t size);
  void alloc_matrix(int64_t num_rows,
                    int64_t num_cols,
                    int64_t nnz,
                    int64_t values_size,
                    int64_t slice_size);
};

} /// namespace MT

#endif //SPMM_FORMAT_SPARSE_MATRIX_H
