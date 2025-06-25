//
// Created by Zhen Peng on 1/26/2024.
//

#include <iostream>
#include "sparse_matrix.h"


/******* COOMatrix **********/
void MT::COOMatrix::alloc_row_ind(int64_t size) {
  free(row_ind_);
  row_ind_ = (int64_t *) malloc(size * sizeof(*row_ind_));
  if (row_ind_ == nullptr) {
    fprintf(stderr, "Error: %s:%d MT::COOMatrix::alloc_row_ind() failed.\n", __FILE__, __LINE__);
    return;
  }
}

void MT::COOMatrix::alloc_col_ind(int64_t size) {
  free(col_ind_);
  col_ind_ = (int64_t *) malloc(size * sizeof(*col_ind_));
  if (col_ind_ == nullptr) {
    fprintf(stderr, "Error: %s:%d MT::COOMatrix::alloc_col_b_ind() failed.\n", __FILE__, __LINE__);
    return;
  }
}

void MT::COOMatrix::alloc_values(int64_t size) {
  free(values_);
  values_ = (double *) malloc(size * sizeof(*values_));
  if (values_ == nullptr) {
    fprintf(stderr, "Error: %s:%d MT::COOMatrix::alloc_values() failed.\n", __FILE__, __LINE__);
    return;
  }
}

void MT::COOMatrix::alloc_matrix(int64_t num_rows,
                                 int64_t num_cols,
                                 int64_t nnz) {
  num_rows_ = num_rows;
  num_cols_ = num_cols;
  nnz_ = nnz;
  alloc_row_ind(nnz);
  alloc_col_ind(nnz);
  alloc_values(nnz);
}

MT::COOMatrix::COOMatrix(int64_t num_rows,
                         int64_t num_cols,
                         int64_t nnz) {
  alloc_matrix(num_rows,
               num_cols,
               nnz);
}

MT::COOMatrix::~COOMatrix() {
  free(row_ind_);
  free(col_ind_);
  free(values_);
}

/******* End COOMatrix **********/


/******* CSRMatrix **********/

void MT::CSRMatrix::alloc_row_offsets_(int64_t size) {
  free(row_offsets_);
  row_offsets_ = (int64_t *) malloc(size * sizeof(*row_offsets_));
  if (row_offsets_ == nullptr) {
    fprintf(stderr, "Error: %s:%d MT::CSRMatrix::alloc_row_b_offsets_() failed.\n", __FILE__, __LINE__);
    return;
  }
}

void MT::CSRMatrix::alloc_col_ind(int64_t size) {
  free(col_ind_);
  col_ind_ = (int64_t *) malloc(size * sizeof(*col_ind_));
  if (col_ind_ == nullptr) {
    fprintf(stderr, "Error: %s:%d MT::CSRMatrix::alloc_col_b_ind() failed.\n", __FILE__, __LINE__);
    return;
  }
}

void MT::CSRMatrix::alloc_values(int64_t size) {
  free(values_);
  values_ = (double *) malloc(size * sizeof(*values_));
  if (values_ == nullptr) {
    fprintf(stderr, "Error: %s:%d MT::CSRMatrix::alloc_values() failed.\n", __FILE__, __LINE__);
    return;
  }
}

void MT::CSRMatrix::alloc_matrix(int64_t num_rows, int64_t num_cols, int64_t nnz) {
  num_rows_ = num_rows;
  num_cols_ = num_cols;
  nnz_ = nnz;
  alloc_row_offsets_(num_rows + 1);
  alloc_col_ind(nnz);
  alloc_values(nnz);
}

MT::CSRMatrix::CSRMatrix(int64_t num_rows, int64_t num_cols, int64_t nnz) {
  alloc_matrix(num_rows, num_cols, nnz);
}

MT::CSRMatrix::~CSRMatrix() {
  free(row_offsets_);
  free(col_ind_);
  free(values_);
}

/******* End CSRMatrix **********/


/******* BSRMatrix **********/

void MT::BSRMatrix::alloc_row_b_offsets_(int64_t size) {
  free(row_b_offsets_);
  row_b_offsets_ = (int64_t *) malloc(size * sizeof(*row_b_offsets_));
  if (row_b_offsets_ == nullptr) {
    fprintf(stderr, "Error: %s:%d MT::BSRMatrix::alloc_row_b_offsets_() failed.\n", __FILE__, __LINE__);
    return;
  }
}

void MT::BSRMatrix::alloc_col_b_ind(int64_t size) {
  free(col_b_ind_);
  col_b_ind_ = (int64_t *) malloc(size * sizeof(*col_b_ind_));
  if (col_b_ind_ == nullptr) {
    fprintf(stderr, "Error: %s:%d MT::BSRMatrix::alloc_col_b_ind() failed.\n", __FILE__, __LINE__);
    return;
  }
}

void MT::BSRMatrix::alloc_values(int64_t size) {
  free(values_);
  values_ = (double *) malloc(size * sizeof(*values_));
  if (values_ == nullptr) {
    fprintf(stderr, "Error: %s:%d MT::BSRMatrix::alloc_values() failed.\n", __FILE__, __LINE__);
    return;
  }
}

void MT::BSRMatrix::alloc_matrix(int64_t num_brows,
                                 int64_t num_bcols,
                                 int64_t bnnz,
                                 int64_t row_bsize,
                                 int64_t col_bsize,
                                 int64_t nnz) {
  num_b_rows_ = num_brows;
  num_b_cols_ = num_bcols;
  bnnz_ = bnnz;
  row_b_size_ = row_bsize;
  col_b_size_ = col_bsize;
  nnz_ = nnz;

  alloc_row_b_offsets_(num_brows + 1);
  alloc_col_b_ind(bnnz);
  alloc_values(bnnz * row_bsize * col_bsize);
}

MT::BSRMatrix::BSRMatrix(int64_t num_brows, int64_t num_bcols, int64_t bnnz, int64_t row_bsize,
                              int64_t col_bsize, int64_t nnz) {
  alloc_matrix(num_brows,
               num_bcols,
               bnnz,
               row_bsize,
               col_bsize,
               nnz);
}

MT::BSRMatrix::~BSRMatrix() {
  free(row_b_offsets_);
  free(col_b_ind_);
  free(values_);
}

/******* End BSRMatrix **********/


/******* BlockedEllMatrix **********/

void MT::BlockedEllMatrix::alloc_col_b_ind(int64_t size) {
  free(col_b_ind_);
  col_b_ind_ = (int64_t *) malloc(size * sizeof(*col_b_ind_));
  if (col_b_ind_ == nullptr) {
    fprintf(stderr, "Error: %s:%d MT::BlockedEllMatrix::alloc_col_b_ind() failed.\n", __FILE__, __LINE__);
    return;
  }
}

void MT::BlockedEllMatrix::alloc_values(int64_t size) {
  free(values_);
  values_ = (double *) malloc(size * sizeof(*values_));
  if (values_ == nullptr) {
    fprintf(stderr, "Error: %s:%d MT::BlockedEllMatrix::alloc_values() failed.\n", __FILE__, __LINE__);
    return;
  }
}

void MT::BlockedEllMatrix::alloc_matrix(int64_t num_rows, int64_t num_cols, int64_t bsize, int64_t ell_cols, int64_t nnz) {
  num_rows_ = num_rows;
  num_cols_ = num_cols;
  bsize_ = bsize;
  ell_cols_ = ell_cols;
  nnz_ = nnz;

  int64_t num_brows = (num_rows + bsize - 1) / bsize;
  int64_t num_bcols = (ell_cols + bsize - 1) / bsize;
  alloc_col_b_ind(num_brows * num_bcols);
  alloc_values(num_rows * ell_cols);
}

MT::BlockedEllMatrix::BlockedEllMatrix(int64_t num_rows, int64_t num_cols, int64_t bsize, int64_t ell_cols, int64_t nnz) {
  alloc_matrix(num_rows,
               num_cols,
               bsize,
               ell_cols,
               nnz);
}

MT::BlockedEllMatrix::~BlockedEllMatrix() {
  free(col_b_ind_);
  free(values_);
}

/******* End BlockedEllMatrix **********/


/******* SlicedEllMatrix **********/

void MT::SlicedEllMatrix::alloc_slice_offsets(int64_t size) {
  free(slice_offsets_);
  slice_offsets_ = (int64_t *) malloc(size * sizeof(*slice_offsets_));
  if (slice_offsets_ == nullptr) {
    fprintf(stderr, "Error: %s:%d MT::SlicedEllMatrix::alloc_slice_offsets() failed.\n", __FILE__, __LINE__);
    return;
  }
}

void MT::SlicedEllMatrix::alloc_col_ind(int64_t size) {
  free(col_ind_);
  col_ind_ = (int64_t *) malloc(size * sizeof(*col_ind_));
  if (col_ind_ == nullptr) {
    fprintf(stderr, "Error: %s:%d MT::SlicedEllMatrix::alloc_col_ind() failed.\n", __FILE__, __LINE__);
    return;
  }
}

void MT::SlicedEllMatrix::alloc_values(int64_t size) {
  free(values_);
  values_ = (double *) malloc(size * sizeof(*values_));
  if (values_ == nullptr) {
    fprintf(stderr, "Error: %s:%d MT::SlicedEllMatrix::alloc_values() failed.\n", __FILE__, __LINE__);
    return;
  }
}

void MT::SlicedEllMatrix::alloc_matrix(int64_t num_rows, int64_t num_cols, int64_t nnz, int64_t values_size,
                                       int64_t slice_size) {
  num_rows_ = num_rows;
  num_cols_ = num_cols;
  nnz_ = nnz;
  values_size_ = values_size;
  slice_size_ = slice_size;

  alloc_slice_offsets((num_rows + slice_size - 1) / slice_size + 1);
  alloc_col_ind(values_size);
  alloc_values(values_size);
}

MT::SlicedEllMatrix::SlicedEllMatrix(int64_t num_rows, int64_t num_cols, int64_t nnz, int64_t values_size,
                                          int64_t slice_size) {
  alloc_matrix(num_rows,
               num_cols,
               nnz,
               values_size,
               slice_size);
}

MT::SlicedEllMatrix::~SlicedEllMatrix() {
  free(slice_offsets_);
  free(col_ind_);
  free(values_);
}

/******* End SlicedEllMatrix **********/