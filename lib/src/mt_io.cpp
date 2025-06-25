//
// Created by Zhen peng on 1/27/2024.
//

#include <algorithm>
#include "mt_io.h"
#include "utils.h"


/// Read Matrix Market file
int MT::read_mm(const char *filename,
                MM_typecode &matcode,
                int &M, /// number of rows
                int &N, /// number of cols
                int &nz,  /// number of non-zeros
                std::vector< std::vector< std::pair<int, double> > > &edgelist) {
  int ret_code;
//  MM_typecode matcode;
  FILE *f;

  if ((f = fopen(filename, "r")) == nullptr) {
    fprintf(stderr, "Error: can't open file %s .\n", filename);
    exit(1);
  }

  printf("\nLoading %s ...\n", filename);

  if (mm_read_banner(f, &matcode) != 0)
  {
    printf("Could not process Matrix Market banner.\n");
    exit(1);
  }


  /*  This is how one can screen matrix types if their application */
  /*  only supports a subset of the Matrix Market data types.      */

  if (mm_is_complex(matcode) && mm_is_matrix(matcode) &&
      mm_is_sparse(matcode) )
  {
    printf("Sorry, this application does not support ");
    printf("Market Market type: [%s]\n", mm_typecode_to_str(matcode));
    exit(1);
  }

//  /// Check if the matrix is symmetric
  bool is_symmetric = mm_is_symmetric(matcode);

  /* find out size of sparse matrix .... */
  if ((ret_code = mm_read_mtx_crd_size(f, &M, &N, &nz)) !=0)
    exit(ret_code);

  /* NOTE: when reading in doubles, ANSI C requires the use of the "l"  */
  /*   specifier as in "%lg", "%lf", "%le", otherwise errors will occur */
  /*  (ANSI C X3.159-1989, Sec. 4.9.6.2, p. 136 lines 13-15)            */

  mm_write_banner(stdout, matcode);
  mm_write_mtx_crd_size(stdout, M, N, nz);
  int num_symmetric_edge = 0;
  edgelist.resize(M);
//  for (int i = 0; i < nz; ++i) {
//    int from;
//    int to;
//    double v;
//    fscanf(f, "%d %d %lg\n", &from, &to, &v);
//    --from;  // adjust from 1-based to 0-based
//    --to;
////    printf("i: %d from: %d to: %d v: %lf\n", i, from+1, to+1, v);
////    edgelist[from].emplace_back(to, v);
//    edgelist[from].push_back(std::make_pair(to, v));  /// Notes: emplace_back would get errors for the shipset1 matrix.
//
//    if (is_symmetric && from != to) {
//      edgelist[to].emplace_back(from, v);
//      ++num_symmetric_edge;
//    }
//
//    if (i < 9) {
//      printf("%d %d %20.19g\n", from + 1, to + 1, v);
//    } else if (i == 9) {
//      printf("( ... )\n");
//    }
//  }

  char buffer[65536];
  int count = 0;
  while (fgets(buffer, sizeof(buffer), f) != nullptr) {
    int from;
    int to;
    double v;
    int num_nums = sscanf(buffer, "%d %d %lg\n", &from, &to, &v);
    if (num_nums == 2) {
      /// Some matrices do not contains values of edges.
      v = 17.0;
    } else if (num_nums != 3) {
      fprintf(stderr, "Error: %s:%d the line (%s) in the matrix (%s) is broken, does not contain 2 or 3 numbers.\n",
              __FILE__, __LINE__, buffer, filename);
      exit(EXIT_FAILURE);
    }
    --from;  // adjust from 1-based to 0-based
    --to;
//    edgelist[from].emplace_back(to, v);
    edgelist[from].push_back(std::make_pair(to, v));  /// Notes: emplace_back would get errors for the shipset1 matrix.

    if (is_symmetric && from != to) {
      edgelist[to].emplace_back(from, v);
      ++num_symmetric_edge;
    }

    if (count < 9) {
      printf("%d %d %20.19g\n", from + 1, to + 1, v);
    } else if (count == 9) {
      printf("( ... )\n");
    }

    ++count;
  }
  if (count != nz) {
    fprintf(stderr, "Error: %s:%d matrix %s only contains %d edges, not %d as claimed in the header.\n",
            __FILE__, __LINE__, filename, count, nz);
    exit(EXIT_FAILURE);
  }

  if (is_symmetric) {
    nz += num_symmetric_edge;
  }

  /// Sort the edges
  for (int row_ind = 0; row_ind < M; ++row_ind) {
    std::sort(edgelist[row_ind].begin(), edgelist[row_ind].end());
  }

  /// Close the file
  fclose(f);

  return ret_code;
}

/// Read COO Matrix
void MT::read_coo_matrix(const char *filename,
                         COOMatrix *coo_matrix) {
  /// Read Matrix Market file
  MM_typecode matcode;
  int M;  // number of rows
  int N;  // number of columns
  int nz;  // number of non-zeros
//  std::vector< std::vector<int> > edgelist;
//  std::vector< std::vector<double> > value_list;
  std::vector< std::vector< std::pair<int, double> > > edgelist;
  read_mm(filename,
          matcode,
          M,
          N,
          nz,
          edgelist);

  /// Out-of-memory check
  {
    int64_t memory_size = 3 * nz * sizeof(*coo_matrix->values_) / 1000000000;
    if (memory_size > MAX_MEMORY_GB) {
      fprintf(stderr, "Error: %s:%d COO matrix requires %lld GB memory, larger than the limit %lld BG.\n",
              __FILE__, __LINE__, memory_size, MAX_MEMORY_GB);
      return;
    }
  }
  /// Allocate coo_matrix
  coo_matrix->alloc_matrix(M, N, nz);

  /// Copy edgelist to coo_matrix
  int offset = 0;
  for (int row = 0; row < M; ++row) {
    for (int col_i = 0; col_i < edgelist[row].size(); ++col_i) {
      int col = edgelist[row][col_i].first;
      double val = edgelist[row][col_i].second;
//      double val = value_list[row][col_i];
      coo_matrix->row_ind_[offset] = row;
      coo_matrix->col_ind_[offset] = col;
      coo_matrix->values_[offset] = val;
      ++offset;
    }
  }

//  /************************/
//  /* now write out matrix */
//  /************************/
//
//  printf("\nPrinting COO matrix...\n");
//  mm_write_banner(stdout, matcode);
//  mm_write_mtx_crd_size(stdout, M, N, nz);
//  for (int i = 0; i < nz; ++i)
//    fprintf(stdout, "%lld %lld %20.19g\n",
//            coo_matrix->row_ind_[i] + 1,
//            coo_matrix->col_ind_[i] + 1,
//            coo_matrix->values_[i]);
}


/// Read CSR Matrix
void MT::read_csr_matrix(const char *filename,
                         CSRMatrix *csr_matrix) {
  /// Read Matrix Market file
  MM_typecode matcode;
  int M;  // number of rows
  int N;  // number of columns
  int nz;  // number of non-zeros
//  std::vector< std::vector<int> > edgelist;
//  std::vector< std::vector<double> > value_list;
  std::vector< std::vector< std::pair<int, double> > > edgelist;
  read_mm(filename,
          matcode,
          M,
          N,
          nz,
          edgelist);

  /// Out-of-memory check
  {
    int64_t memory_size = nz * sizeof(*csr_matrix->values_) / 1000000000;
    if (memory_size > MAX_MEMORY_GB) {
      fprintf(stderr, "Error: %s:%d CSR matrix requires %lld GB memory, larger than the limit %lld BG.\n",
              __FILE__, __LINE__, memory_size, MAX_MEMORY_GB);
      return;
    }
  }
  /// Allocate csr_matrix
  csr_matrix->alloc_matrix(M, N, nz);

  /// Copy edgelist to csr_matrix
  int offset = 0;
  for (int row = 0; row < M; ++row) {
    csr_matrix->row_offsets_[row] = offset;
    for (int col_i = 0; col_i < edgelist[row].size(); ++col_i) {
      int col = edgelist[row][col_i].first;
      double val = edgelist[row][col_i].second;
//      double val = value_list[row][col_i];
      csr_matrix->col_ind_[offset] = col;
      csr_matrix->values_[offset] = val;
      ++offset;
    }
  }
  csr_matrix->row_offsets_[M] = offset;

//  /************************/
//  /* now write out matrix */
//  /************************/
//
//  printf("\nPrinting CSR matrix...\n");
//  mm_write_banner(stdout, matcode);
//  mm_write_mtx_crd_size(stdout, M, N, nz);
//  int print_count = 0;
//  int print_limit = 9;
//  for (int64_t row = 0; row < csr_matrix->num_rows_ && print_count < print_limit; ++row) {
//    for (int64_t col_loc = csr_matrix->row_offsets_[row]; col_loc < csr_matrix->row_offsets_[row + 1]; ++col_loc) {
//      int64_t col = csr_matrix->col_ind_[col_loc];
//      double val = csr_matrix->values_[col_loc];
//      fprintf(stdout, "%lld %lld %20.19g\n",
//              row + 1,
//              col + 1,
//              val);
//      ++print_count;
//      if (print_count == print_limit) {
//        printf("( ... )\n");
//        break;
//      }
//    }
//  }
}


/// Read BSR (blocked compressed sparse row) Matrix
void MT::read_bsr_matrix(const char *filename,
                         int64_t row_b_size,  /// number of rows in each block
                         int64_t col_b_size,  /// number of columns in each block
                         BSRMatrix *bsr_matrix  /* output */) {
  /// Read Matrix Market file
  MM_typecode matcode;
  int M;  // number of rows
  int N;  // number of columns
  int nz;  // number of non-zeros
  std::vector< std::vector< std::pair<int, double> > > edgelist;
  read_mm(filename,
          matcode,
          M,
          N,
          nz,
          edgelist);

//  int64_t row_bsize = 8;  /// number of rows of each block
//  int64_t col_bsize = 8;  /// number of columns of each block
  int64_t num_brows = (M + row_b_size - 1) / row_b_size;  /// number of block rows
  int64_t num_bcols = (N + col_b_size - 1) / col_b_size;  /// number of block columns

  /// Find out non-zero blocks
  std::vector< std::vector<bool> > is_non_zero_block(num_brows, std::vector<bool>(num_bcols, false));
  int64_t bnnz = 0;
  for (int64_t row_ind = 0; row_ind < M; ++row_ind) {
    int64_t brow_ind = row_ind / row_b_size;
    for (auto &edge : edgelist[row_ind]) {
      int64_t col_ind = edge.first;
      int64_t bcol_ind = col_ind / col_b_size;
      if (!is_non_zero_block[brow_ind][bcol_ind]) {
        is_non_zero_block[brow_ind][bcol_ind] = true;
        ++bnnz;
      }
    }
  }

  /// Out-of-memory check
  {
    int64_t memory_size = bnnz * row_b_size * col_b_size * sizeof(*bsr_matrix->values_) / 1000000000;
    if (memory_size > MAX_MEMORY_GB) {
      fprintf(stderr, "Error: %s:%d BSR matrix requires %lld GB memory, larger than the limit %lld BG.\n",
              __FILE__, __LINE__, memory_size, MAX_MEMORY_GB);
      return;
    }
  }

  /// Allocate the sparse matrix
  bsr_matrix->alloc_matrix(num_brows, num_bcols, bnnz, row_b_size, col_b_size, nz);

  /// Populate block_offsets_ and col_b_ind_
  int64_t block_loc = 0;
  std::vector< std::vector<int64_t> > non_zero_block_ind(num_brows, std::vector<int64_t>(num_bcols, 0));  /// /// There are non_zero_block_ind[brow_ind][bcol_ind] non-zero blocks before the blocks[brow_ind][bcol_ind];
  for (int64_t brow_ind = 0; brow_ind < num_brows; ++brow_ind) {
    bsr_matrix->row_b_offsets_[brow_ind] = block_loc;
    for (int64_t bcol_ind = 0; bcol_ind < num_bcols; ++bcol_ind) {
      if (is_non_zero_block[brow_ind][bcol_ind]) {
        bsr_matrix->col_b_ind_[block_loc] = bcol_ind;
        non_zero_block_ind[brow_ind][bcol_ind] = block_loc;
        ++block_loc;
      }
    }
  }
  bsr_matrix->row_b_offsets_[num_brows] = block_loc;

  /// Populate values
  memset(bsr_matrix->values_, 0, bnnz * row_b_size * col_b_size * sizeof(*bsr_matrix->values_));
  for (int64_t row_ind = 0; row_ind < M; ++row_ind) {
    int64_t brow_ind = row_ind / row_b_size;
    int64_t ele_row_ind = row_ind % row_b_size;
    for (auto &edge : edgelist[row_ind]) {
      int64_t col_ind = edge.first;
      double value = edge.second;
      int64_t bcol_ind = col_ind / col_b_size;
      int64_t ele_col_ind = col_ind % col_b_size;
      int64_t block_ind = non_zero_block_ind[brow_ind][bcol_ind];
      int64_t val_loc = (ele_row_ind * col_b_size + ele_col_ind) + block_ind * row_b_size * col_b_size;
      bsr_matrix->values_[val_loc] = value;
    }
  }
}

/// Backup
///// Read BSR (blocked compressed sparse row) Matrix
//void MT::read_bsr_matrix(const char *filename,
//                         int64_t row_b_size,  /// number of rows in each block
//                         int64_t col_b_size,  /// number of columns in each block
//                         BSRMatrix *bsr_matrix  /* output */) {
//  /// Read Matrix Market file
//  MM_typecode matcode;
//  int M;  // number of rows
//  int N;  // number of columns
//  int nz;  // number of non-zeros
//  std::vector< std::vector< std::pair<int, double> > > edgelist;
//  read_mm(filename,
//          matcode,
//          M,
//          N,
//          nz,
//          edgelist);
//
//  /// Allocate csr_matrix
////  int64_t row_bsize = 8;  /// number of rows of each block
////  int64_t col_bsize = 8;  /// number of columns of each block
//  int64_t num_brows = (M + row_b_size - 1) / row_b_size;  /// number of block rows
//  int64_t num_bcols = (N + col_b_size - 1) / col_b_size;  /// number of block columns
//
//  /// Read from edgelist to buffer. The buffer stores each value into its block
//  int64_t block_capacity = row_b_size * col_b_size;
//  std::vector< std::vector< std::vector<double> > > block_value_matrix(num_brows, std::vector< std::vector<double> >(num_bcols, std::vector<double>()));
//  int64_t bnnz = 0;  /// number of non-zero blocks
//  for (int64_t row_ind = 0; row_ind < M; ++row_ind) {
//    int64_t brow_ind = row_ind / row_b_size;
//    int64_t ele_row_ind = row_ind % row_b_size;
//    for (auto &edge : edgelist[row_ind]) {
//      int col_ind = edge.first;
//      double value = edge.second;
//      int64_t bcol_ind = col_ind / col_b_size;
//      int64_t ele_col_ind = col_ind % col_b_size;
//
//      /// Add the value to its block
//      auto &block = block_value_matrix[brow_ind][bcol_ind];
//      if (block.empty()) {
//        block.resize(block_capacity, 0.0);
//        ++bnnz;
//      }
//      block[ele_row_ind * col_b_size + ele_col_ind] = value;
//    }
//  }
//  edgelist.clear();
//
//  /// Allocate the sparse matrix
//  bsr_matrix->alloc_matrix(num_brows, num_bcols, bnnz, row_b_size, col_b_size, nz);
//
//  /// Store to the BSR matrix
//  int64_t val_loc = 0;
//  int64_t bcol_loc = 0;
//  int brow_offset = 0;
//  for (int64_t brow_ind = 0; brow_ind < num_brows; ++brow_ind) {
//    bsr_matrix->row_b_offsets_[brow_ind] = brow_offset;
//    for (int64_t bcol_ind = 0; bcol_ind < num_bcols; ++bcol_ind) {
//      auto &block = block_value_matrix[brow_ind][bcol_ind];
//      if (block.empty()) {
//        continue;
//      }
//      bsr_matrix->col_b_ind_[bcol_loc++] = bcol_ind;
//      for (double value : block) {
//        bsr_matrix->values_[val_loc++] = value;
//      }
//      ++brow_offset;
//    }
//  }
//  bsr_matrix->row_b_offsets_[num_brows] = brow_offset;
//}
/// End backup

/// Read Blocked Ellpack matrix
void MT::read_blocked_ell_matrix(const char *filename,
                                 int64_t bsize,
                                 MT::BlockedEllMatrix *blocked_ell_matrix) {
  /// Read Matrix Market file
  MM_typecode matcode;
  int M;  // number of rows
  int N;  // number of columns
  int nz;  // number of non-zeros
  std::vector< std::vector< std::pair<int, double> > > edgelist;
  read_mm(filename,
          matcode,
          M,
          N,
          nz,
          edgelist);

  /// Read from edgelist to buffer. The buffer stores each value into its block
//  { /// test
//    printf("#### %s:%d Reading block_value_matrix...\n", __FILE__, __LINE__);
//  }
  int64_t num_b_rows = (M + bsize - 1) / bsize;
  int64_t num_b_cols = (N + bsize - 1) / bsize;  /// number of block columns of the sparse matrix

  /// Get non-zero-block indicators
  std::vector< std::vector<bool> > is_non_zero_block(num_b_rows, std::vector<bool>(num_b_cols, false));
  for (int64_t row_ind = 0; row_ind < M; ++row_ind) {
    int64_t brow_ind = row_ind / bsize;
    for (auto &edge : edgelist[row_ind]) {
      int64_t col_ind = edge.first;
      int64_t bcol_ind = col_ind / bsize;
      if (!is_non_zero_block[brow_ind][bcol_ind]) {
        is_non_zero_block[brow_ind][bcol_ind] = true;
      }
    }
  }

  /// Calculate the ell_cols
//  { /// test
//    printf("#### %s:%d Find out the value of col_ells_...\n", __FILE__, __LINE__);
//  }
  int64_t max_non_zero_bcol = 0;
  for (int64_t brow_ind = 0; brow_ind < num_b_rows; ++brow_ind) {
    int64_t non_zero_bcol = 0;
    for (int64_t bcol_ind = 0; bcol_ind < num_b_cols; ++bcol_ind) {
      if (is_non_zero_block[brow_ind][bcol_ind]) {
        ++non_zero_bcol;
      }
    }
    if (non_zero_bcol > max_non_zero_bcol) {
      max_non_zero_bcol = non_zero_bcol;
    }
  }
  int64_t ell_cols = max_non_zero_bcol * bsize;

  /// Allocate
//  { /// test
//    printf("#### %s:%d Allocating BLOCKED-ELL matrix...\n", __FILE__, __LINE__);
//  }
  /// Size after Padding
  int64_t num_rows = num_b_rows * bsize;
  int64_t num_cols = num_b_cols * bsize;
  int64_t num_ell_b_cols = max_non_zero_bcol;  /// number of block columns of the Ellpack matrix
//  int64_t num_ell_b_cols = ell_cols / bsize;  /// number of block columns of the Ellpack matrix
  /// Out-of-memory check
  {
    int64_t memory_size = num_rows * ell_cols * sizeof(*blocked_ell_matrix->values_) / 1000000000;
    if (memory_size > MAX_MEMORY_GB) {
      fprintf(stderr, "Error: %s:%d BLOCKED-ELL requires %lld GB memory, larger than the limit %lld BG.\n",
      __FILE__, __LINE__, memory_size, MAX_MEMORY_GB);
      return;
    }
  }
  blocked_ell_matrix->alloc_matrix(num_rows, num_cols, bsize, ell_cols, nz);

  /// Store to the Blocked Ellpack matrix's column block indices (col_b_ind_)
//  { /// test
//    printf("#### %s:%d Storing col_b_ind_ ...\n", __FILE__, __LINE__);
//  }
  std::vector< std::vector<int64_t> > non_zero_block_col_ind(num_b_rows, std::vector<int64_t>(num_b_cols, 0)); /// In row brow_ind, there are non_zero_block_col_ind[brow_ind][bcol_ind] non-zero blocks before the blocks[brow_ind][bcol_ind];
  for (int64_t brow_ind = 0; brow_ind < num_b_rows; ++brow_ind) {
    int64_t bcol_loc_base = brow_ind * num_ell_b_cols;
    int64_t bcol_loc = 0;
    int64_t non_zero_bcol_ind = 0;
    for (int64_t bcol_ind = 0; bcol_ind < num_b_cols; ++bcol_ind) {
      if (is_non_zero_block[brow_ind][bcol_ind]) {
        blocked_ell_matrix->col_b_ind_[bcol_loc_base + bcol_loc++] = bcol_ind;
        non_zero_block_col_ind[brow_ind][bcol_ind] = non_zero_bcol_ind++;
      }
    }
    while (bcol_loc < num_ell_b_cols) {
      blocked_ell_matrix->col_b_ind_[bcol_loc_base + bcol_loc++] = -1;
    }
  }

  /// Store to the Blocked Ellpack matrix's values (values_)
//  { /// test
//    printf("#### %s:%d Storing values_ ...\n", __FILE__, __LINE__);
//    printf("num_rows_: %lld "
//           "ell_cols_: %lld "
//           "values_size: %lld\n",
//           blocked_ell_matrix->num_rows_,
//           blocked_ell_matrix->ell_cols_,
//           blocked_ell_matrix->num_rows_ * blocked_ell_matrix->ell_cols_);
//  }
  /// All values set 0 initially.
  memset(blocked_ell_matrix->values_, 0, num_rows * ell_cols * sizeof(*blocked_ell_matrix->values_));
  for (int64_t row_ind = 0; row_ind < M; ++row_ind) {
    int64_t brow_ind = row_ind / bsize;
    int64_t row_ind_base = row_ind * ell_cols;
    for (auto &edge : edgelist[row_ind]) {
      int64_t col_ind = edge.first;
      double value = edge.second;
      int64_t bcol_ind = col_ind / bsize;
      int64_t ele_col_ind = col_ind % bsize;
      int64_t val_loc = row_ind_base + non_zero_block_col_ind[brow_ind][bcol_ind] * bsize + ele_col_ind;
      blocked_ell_matrix->values_[val_loc] = value;
    }
  }
}

/// Backup
///// Read Blocked Ellpack matrix
//void MT::read_blocked_ell_matrix(const char *filename,
//                                 int64_t bsize,
//                                 MT::BlockedEllMatrix *blocked_ell_matrix) {
//  /// Read Matrix Market file
//  MM_typecode matcode;
//  int M;  // number of rows
//  int N;  // number of columns
//  int nz;  // number of non-zeros
//  std::vector< std::vector< std::pair<int, double> > > edgelist;
//  read_mm(filename,
//          matcode,
//          M,
//          N,
//          nz,
//          edgelist);
//
//  /// Read from edgelist to buffer. The buffer stores each value into its block
//  { /// test
//    printf("#### %s:%d Reading block_value_matrix...\n", __FILE__, __LINE__);
//  }
//  int64_t num_b_rows = (M + bsize - 1) / bsize;
//  int64_t num_b_cols = (N + bsize - 1) / bsize;  /// number of block columns of the sparse matrix
//  std::vector< std::vector< std::vector<double> > > block_value_matrix(num_b_rows, std::vector< std::vector<double> >(num_b_cols, std::vector<double>()));
//
////  int64_t ell_cols = 0;
//  for (int64_t row_ind = 0; row_ind < M; ++row_ind) {
////    if (ell_cols < edgelist[row_ind].size()) {
////      ell_cols = edgelist[row_ind].size();
////    }
//    int64_t brow_ind = row_ind / bsize;
//    int64_t ele_row_ind = row_ind % bsize;
//    for (auto &edge : edgelist[row_ind]) {
//      int64_t col_ind = edge.first;
//      double value = edge.second;
//      int64_t bcol_ind = col_ind / bsize;
//      int64_t ele_col_ind = col_ind % bsize;
//
//      auto &block = block_value_matrix[brow_ind][bcol_ind];
//      if (block.empty()) {
//        block.resize(bsize * bsize, 0.0);
//      }
//      block[ele_row_ind * bsize + ele_col_ind] = value;
//    }
//  }
//  edgelist.clear();
//
//  /// Find out the value of col_ells_
//  { /// test
//    printf("#### %s:%d Find out the value of col_ells_...\n", __FILE__, __LINE__);
//  }
//  int64_t max_nonzero_blocks_row = 0;
//  for (int64_t brow_ind = 0; brow_ind < num_b_rows; ++brow_ind) {
//    int non_zero_blocks = 0;
//    for (int64_t bcol_ind = 0; bcol_ind < num_b_cols; ++bcol_ind) {
//      if (!block_value_matrix[brow_ind][bcol_ind].empty()) {
//        ++non_zero_blocks;
//      }
//    }
//    if (max_nonzero_blocks_row < non_zero_blocks) {
//      max_nonzero_blocks_row = non_zero_blocks;
//    }
//  }
//  int64_t ell_cols = max_nonzero_blocks_row * bsize;
//
//  /// Allocate the sparse matrix
//  { /// test
//    printf("#### %s:%d Allocating BLOCKED-ELL matrix...\n", __FILE__, __LINE__);
//  }
//  /// Size after Padding
//  int64_t num_rows = num_b_rows * bsize;
//  int64_t num_cols = num_b_cols * bsize;
//  int64_t num_ell_b_cols = ell_cols / bsize;  /// number of block columns of the Ellpack matrix
////  ell_cols = num_ell_b_cols * bsize;
//  blocked_ell_matrix->alloc_matrix(num_rows, num_cols, bsize, ell_cols, nz);
//
//  /// Store to the Blocked Ellpack matrix's values (values_)
//  { /// test
//    printf("#### %s:%d Storing values_ ...\n", __FILE__, __LINE__);
//    printf("num_rows_: %lld "
//           "ell_cols_: %lld "
//           "values_size: %lld\n",
//           blocked_ell_matrix->num_rows_,
//           blocked_ell_matrix->ell_cols_,
//           blocked_ell_matrix->num_rows_ * blocked_ell_matrix->ell_cols_);
//  }
//  for (int64_t brow_ind = 0; brow_ind < num_b_rows; ++brow_ind) {
//    { /// test
//      printf("brow_ind: %lld\n", brow_ind);
//    }
//    int64_t row_ind_base = brow_ind * bsize;
//    for (int64_t ele_row_ind = 0; ele_row_ind < bsize; ++ele_row_ind) {
//      int64_t row_ind = ele_row_ind + row_ind_base;
//      int64_t value_loc = row_ind * ell_cols;
//      for (int64_t bcol_ind = 0; bcol_ind < num_b_cols; ++bcol_ind) {
//        auto &block = block_value_matrix[brow_ind][bcol_ind];
//        if (block.empty()) {
//          continue;
//        }
//        for (int64_t ele_col_ind = 0; ele_col_ind < bsize; ++ele_col_ind) {
//          blocked_ell_matrix->values_[value_loc++] = block[ele_row_ind * bsize + ele_col_ind];
//        }
//      }
////      /// Padding at the end of the row
////      while (value_loc < (row_ind + 1) * ell_cols) {
////        blocked_ell_matrix->values_[value_loc++] = 0.0;
////      }
//    }
//  }
//
//  /// Store to the Blocked Ellpack matrix's column block indices (col_b_ind_)
//  { /// test
//    printf("#### %s:%d Storing col_b_ind_ ...\n", __FILE__, __LINE__);
//  }
//  for (int64_t brow_ind = 0; brow_ind < num_b_rows; ++brow_ind) {
//    int64_t bcol_loc_base = brow_ind * num_ell_b_cols;
//    int64_t bcol_loc = 0;
//    for (int64_t bcol_ind = 0; bcol_ind < num_b_cols; ++bcol_ind) {
//      if (block_value_matrix[brow_ind][bcol_ind].empty()) {
//        continue;
//      }
//      blocked_ell_matrix->col_b_ind_[bcol_loc_base + bcol_loc++] = bcol_ind;
//    }
//    while (bcol_loc < num_ell_b_cols) {
//      blocked_ell_matrix->col_b_ind_[bcol_loc_base + bcol_loc++] = -1;
//    }
//  }
//}
/// End backup


/// Read Sliced Ellpack matrix
void MT::read_sliced_ell_matrix(const char *filename,
                                int64_t slice_size,
                                MT::SlicedEllMatrix *sliced_ell_matrix) {
  /// Read Matrix Market file
  MM_typecode matcode;
  int M;  // number of rows
  int N;  // number of columns
  int nz;  // number of non-zeros
  std::vector< std::vector< std::pair<int, double> > > edgelist;
  read_mm(filename,
          matcode,
          M,
          N,
          nz,
          edgelist);

  /// Find out the number of ellpack columns for each slice
  int64_t num_slices = (M + slice_size - 1) / slice_size;
  int64_t num_rows = num_slices * slice_size;
  std::vector<int64_t> slice_col_sizes(num_slices, 0);
  for (int64_t row_ind = 0; row_ind < M; ++row_ind) {
    int64_t slice_ind = row_ind / slice_size;
    if (edgelist[row_ind].size() > slice_col_sizes[slice_ind]) {
      slice_col_sizes[slice_ind] = edgelist[row_ind].size();
    }
  }

  /// Find the size of column indices and values
  int64_t values_size = 0;
  for (int64_t slice_ind = 0; slice_ind < num_slices; ++slice_ind) {
    values_size += slice_col_sizes[slice_ind] * slice_size;
  }

  /// Out-of-memory check
  {
    int64_t memory_size = values_size * sizeof(*sliced_ell_matrix->values_) / 1000000000;
    if (memory_size > MAX_MEMORY_GB) {
      fprintf(stderr, "Error: %s:%d SLICED-ELL matrix requires %lld GB memory, larger than the limit %lld BG.\n",
              __FILE__, __LINE__, memory_size, MAX_MEMORY_GB);
      return;
    }
  }

  /// Allocate the matrix
  sliced_ell_matrix->alloc_matrix(num_rows, N, nz, values_size, slice_size);

  /// Populate elements into the matrix
  int64_t val_loc = 0;
  for (int64_t slice_ind = 0; slice_ind < num_slices; ++slice_ind) {
    sliced_ell_matrix->slice_offsets_[slice_ind] = val_loc;
    int64_t slice_col_size = slice_col_sizes[slice_ind];
    int64_t row_ind_base = slice_ind * slice_size;
    for (int64_t col_loc = 0; col_loc < slice_col_size; ++col_loc) {
      for (int64_t ele_row_ind = 0; ele_row_ind < slice_size; ++ele_row_ind) {
        int64_t row_ind = ele_row_ind + row_ind_base;
        if (row_ind < edgelist.size()) {
          if (col_loc < edgelist[row_ind].size()) {
            auto &edge = edgelist[row_ind][col_loc];
            sliced_ell_matrix->col_ind_[val_loc] = edge.first;
            sliced_ell_matrix->values_[val_loc] = edge.second;
          } else {
            /// Padding the row ending
            sliced_ell_matrix->col_ind_[val_loc] = -1;
            sliced_ell_matrix->values_[val_loc] = 0.0;
          }
        } else {
          /// Padding the empty row
          sliced_ell_matrix->col_ind_[val_loc] = -1;
          sliced_ell_matrix->values_[val_loc] = 0.0;
        }
        ++val_loc;
      }
    }
  }
  sliced_ell_matrix->slice_offsets_[num_slices] = val_loc;

}

/// Backup
///// Read Sliced Ellpack matrix
//void MT::read_sliced_ell_matrix(const char *filename,
//                                int64_t slice_size,
//                                MT::SlicedEllMatrix *sliced_ell_matrix) {
//  /// Read Matrix Market file
//  MM_typecode matcode;
//  int M;  // number of rows
//  int N;  // number of columns
//  int nz;  // number of non-zeros
//  std::vector< std::vector< std::pair<int, double> > > edgelist;
//  read_mm(filename,
//          matcode,
//          M,
//          N,
//          nz,
//          edgelist);
//
//  /// Find out the number of ellpack columns for each slice
//  int64_t num_slices = (M + slice_size - 1) / slice_size;
//  int64_t num_rows = num_slices * slice_size;
//  std::vector<int64_t> slice_col_sizes(num_slices, 0);
////  int64_t nnz = 0;  /// number of non-zeros
//  for (int64_t row_ind = 0; row_ind < M; ++row_ind) {
//    int64_t slice_ind = row_ind / slice_size;
//    if (edgelist[row_ind].size() > slice_col_sizes[slice_ind]) {
//      slice_col_sizes[slice_ind] = edgelist[row_ind].size();
//    }
////    nnz += edgelist[row_ind].size();
//  }
//
//
////  std::vector< std::vector< std::pair<int, double> > > ellpack_edgelist(num_rows, std::vector< std::pair<int, double> >());
//  /// Padding the edgelist
//  int64_t row_ind = 0;
//  int64_t values_size = 0; /// total number of elements (nonzero and padding)
//  for ( ; row_ind < M; ++row_ind) {
//    int64_t slice_ind = row_ind / slice_size;
//    int64_t ell_col_size = slice_col_sizes[slice_ind];
//    int64_t row_size = edgelist[row_ind].size();
//    for (int64_t count = 0; count < ell_col_size - row_size; ++count) {
////      edgelist[row_ind].emplace_back(-1, 0.0);
//      edgelist[row_ind].push_back(std::make_pair(-1, 0.0));
//    }
//    values_size += edgelist[row_ind].size();
//  }
//  /// Padding extra rows
//  for ( ; row_ind < num_rows; ++row_ind) {
//    int64_t slice_ind = row_ind / slice_size;
//    int64_t ell_col_size = slice_col_sizes[slice_ind];
////    edgelist.push_back(std::vector< std::pair<int, double> >(ell_col_size, std::make_pair(-1, 0.0)));
//    edgelist.emplace_back(ell_col_size, std::make_pair(-1, 0.0));
//    values_size += ell_col_size;
//  }
//
//  /// Allocate the matrix
//  sliced_ell_matrix->alloc_matrix(num_rows, N, nz, values_size, slice_size);
//
//  /// Copy to Sliced Ellpack matrix
//  int64_t val_loc = 0;
//  for (int64_t slice_ind = 0; slice_ind < num_slices; ++slice_ind) {
//    sliced_ell_matrix->slice_offsets_[slice_ind] = val_loc;
//    int64_t slice_row_base = slice_ind * slice_size;
//    int64_t ell_col_size = slice_col_sizes[slice_ind];
//    for (int64_t ele_col_ind = 0; ele_col_ind < ell_col_size; ++ele_col_ind) {
//      for (int64_t ele_row_ind = 0; ele_row_ind < slice_size; ++ele_row_ind) {
//        /// column-major in a slice
//        row_ind = ele_row_ind + slice_row_base;
//        auto &edge = edgelist[row_ind][ele_col_ind];
//        sliced_ell_matrix->col_ind_[val_loc] = edge.first;
//        sliced_ell_matrix->values_[val_loc] = edge.second;
//        ++val_loc;
//      }
//    }
//  }
//  sliced_ell_matrix->slice_offsets_[num_slices] = val_loc;
//}
/// End backup

/// Create a dense matrix
void MT::create_random_dense_matrix(int64_t num_rows,
                                    int64_t num_cols,
                                    DenseMatrix *dense_matrix) {
  srand(17);
  dense_matrix->alloc_matrix(num_rows, num_cols);
  int64_t nnz = num_rows * num_cols;
  for (int64_t num_i = 0; num_i < nnz; ++num_i) {
    dense_matrix->values_[num_i] = 1.7;
//    dense_matrix.values_[num_i] = rand();
  }
}

/// Print COO Matrix
void MT::print_coo_matrix(const MT::COOMatrix *coo_matrix) {
  printf("\nPrinting COO matrix...\n");
  printf("num_rows_: %lld\n", coo_matrix->num_rows_);
  printf("num_cols_: %lld\n", coo_matrix->num_cols_);
  int64_t nnz = coo_matrix->nnz_;
  printf("nnz_: %lld\n", nnz);

  printf("row_ind_:\n");
  for (int64_t i = 0; i < nnz; ++i) {
    printf("%lld ", coo_matrix->row_ind_[i]);
  }
  printf("\n");

  printf("col_ind_:\n");
  for (int64_t i = 0; i < nnz; ++i) {
    printf("%lld ", coo_matrix->col_ind_[i]);
  }
  printf("\n");

  printf("values_:\n");
  for (int64_t i = 0; i < nnz; ++i) {
    printf("%lf ", coo_matrix->values_[i]);
  }
  printf("\n");
}

/// Print CSR Matrix
void MT::print_csr_matrix(const MT::CSRMatrix *csr_matrix) {
  printf("\nPrinting CSR Matrix...\n");
  printf("num_rows_: %lld\n", csr_matrix->num_rows_);
  printf("num_cols_: %lld\n", csr_matrix->num_cols_);
  printf("nnz_: %lld\n", csr_matrix->nnz_);
  printf("row_offsets_:\n");
  for (int64_t i = 0; i < csr_matrix->num_rows_ + 1; ++i) {
    printf("%lld ", csr_matrix->row_offsets_[i]);
  }
  printf("\n");

  printf("col_ind_:\n");
  for (int64_t i = 0; i < csr_matrix->nnz_; ++i) {
    printf("%lld ", csr_matrix->col_ind_[i]);
  }
  printf("\n");

  printf("values_:\n");
  for (int64_t i = 0; i < csr_matrix->nnz_; ++i) {
    printf("%lf ", csr_matrix->values_[i]);
  }
  printf("\n");
}

/// Print BSR Matrix
void MT::print_bsr_matrix(const BSRMatrix *matrix) {
  printf("\nBSR Matrix:\n");
  printf("num_b_rows_ : %lld\n", matrix->num_b_rows_);
  printf("num_b_cols_: %lld\n", matrix->num_b_cols_);
  printf("bnnz_: %lld\n", matrix->bnnz_);
  printf("row_b_size_: %lld\n", matrix->row_b_size_);
  printf("col_b_size_: %lld\n", matrix->col_b_size_);
  printf("matrix->row_b_offsets_:\n");
  for (int64_t i = 0; i < matrix->num_b_rows_ + 1; ++i) {
    printf("%lld ", matrix->row_b_offsets_[i]);
  }
  printf("\n");

  printf("matrix->col_b_ind_:\n");
  for (int64_t i = 0; i < matrix->bnnz_; ++i) {
    printf("%lld ", matrix->col_b_ind_[i]);
  }
  printf("\n");

  printf("matrix->values_:\n");
  for (int64_t i = 0; i < matrix->bnnz_ * matrix->row_b_size_ * matrix->col_b_size_; ++i) {
    printf("%lf ", matrix->values_[i]);
  }
  printf("\n");
}


/// Print Blocked ELlpack matrix
void MT::print_blocked_ell_matrix(const MT::BlockedEllMatrix *blocked_ell_matrix) {
  /// Print column block indices (col_b_ind_)
  printf("\nPrint Blocked Ellpack Matrix...\n");
  printf("num_rows_: %lld\n", blocked_ell_matrix->num_rows_);
  printf("num_cols_: %lld\n", blocked_ell_matrix->num_cols_);
  printf("bsize_: %lld\n", blocked_ell_matrix->bsize_);
  printf("ell_cols_: %lld\n", blocked_ell_matrix->ell_cols_);
  printf("col_b_ind_:\n");

  int64_t bsize = blocked_ell_matrix->bsize_;
  int64_t ell_cols = blocked_ell_matrix->ell_cols_;
  int64_t num_ell_b_cols = ell_cols / bsize;
  int64_t num_rows = blocked_ell_matrix->num_rows_;
  int64_t num_cols = blocked_ell_matrix->num_cols_;
  int64_t num_b_rows = num_rows / bsize;
  for (int64_t brow_ind = 0; brow_ind < num_b_rows; ++brow_ind) {
    for (int64_t bcol_ind = 0; bcol_ind < num_ell_b_cols; ++bcol_ind) {
      printf("%lld ", blocked_ell_matrix->col_b_ind_[brow_ind * num_ell_b_cols + bcol_ind]);
    }
    printf("\n");
  }
//  printf("\n");

  printf("values_:\n");
  for (int64_t row_ind = 0; row_ind < num_rows; ++row_ind) {
    for (int64_t col_ind = 0; col_ind < ell_cols; ++col_ind) {
      printf("%lf ", blocked_ell_matrix->values_[row_ind * ell_cols + col_ind]);
    }
    printf("\n");
  }
}


/// Print Sliced Ellpack Matrix
void MT::print_sliced_ell_matrix(const MT::SlicedEllMatrix *sliced_ell_matrix) {
  printf("\nPrinting Sliced Ellpack matrix...\n");
  printf("num_rows_: %lld\n", sliced_ell_matrix->num_rows_);
  printf("num_cols_: %lld\n", sliced_ell_matrix->num_cols_);
  printf("nnz_: %lld\n", sliced_ell_matrix->nnz_);
  printf("values_size_: %lld\n", sliced_ell_matrix->values_size_);
  printf("slice_size_: %lld\n", sliced_ell_matrix->slice_size_);

  printf("slice_offsets_:\n");
  int64_t num_slice = (sliced_ell_matrix->num_rows_ + sliced_ell_matrix->slice_size_ - 1) / sliced_ell_matrix->slice_size_;
  for (int64_t i = 0; i < num_slice + 1; ++i) {
    printf("%lld ", sliced_ell_matrix->slice_offsets_[i]);
  }
  printf("\n");

  printf("col_ind_:\n");
  for (int64_t i = 0; i < sliced_ell_matrix->values_size_; ++i) {
    printf("%lld ", sliced_ell_matrix->col_ind_[i]);
  }
  printf("\n");

  printf("values_:\n");
  for (int64_t i = 0; i < sliced_ell_matrix->values_size_; ++i) {
    printf("%lf ", sliced_ell_matrix->values_[i]);
  }
  printf("\n");
}


