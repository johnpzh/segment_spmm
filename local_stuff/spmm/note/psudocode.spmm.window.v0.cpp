
void fun()
{
  // C[i,j] = A[i,k] * B[k,j]
  for (int ii = 0; ii < A1; ii += A1_tile) {
    /// Get the boundary
    int i_bound = min(ii + A1_tile, A1);

    int num_tiles;
    table = vector<vector>(i_bound - ii, vector(num_tiles + 1));  /// dims: A1_tile * (num_tiles + 1);
    GetTable(A=A,
             row_idx_start=ii,
             row_idx_bound=i_bound,
             A2_tile=A2_tile,
             num_tiles=num_tiles/*out*/,
             table=table/*out*/);

    for (int kk_tile_i = 0; kk_tile_i < num_tiles; ++kk_tile_i) {
      for (int jj = 0; jj < B2; jj += B2_tile) {
        int j_bound = min(jj + B2_tile, B2);
        /// Tile
        for (int i = ii; i < i_bound; ++i) {
          for (int k_loc = table[i - ii][kk_tile_i]; k_loc < table[i - ii][kk_tile_i + 1]; ++k_loc) {
            int k = A.crds[k_loc];
            for (int j = jj; j < j_bound; ++j) {
              C[i][j] = A.vals[k_loc] * B[k][j];
            }
          }
        }
      }
    }
  }
}


void GetTable(matrix A,
              int row_idx_start,
              int row_idx_bound,
              int A2_tile,
              int &num_tiles,
              vector &table/*out*/)
{
  int left_col_id = INT_MAX;
  int right_col_id = 0;
  for (int i = row_idx_start; i < row_idx_bound; ++i) {
    for (int k_loc = A.pos[i]; k_loc < A.pos[i + 1]; ++k_loc) {
      int k = A.crds[k_loc];
      if (k < left_col_id) {
        left_col_id = k;
      }
      if (k > right_col_id) {
        right_col_id = k;
      }
    }
  }
  num_tiles = ((right_col_id - left_col_id + 1) + (A2_tile - 1)) / A2_tile;

  for (int row_idx = row_idx_start; row_idx < row_idx_bound; ++row_idx) {
    int k_loc_start = A.pos[row_idx];
    int k_loc_bound = A.pos[row_idx + 1];
    int curr_k_loc = k_loc_start;
    table[row_idx - row_idx_start][0] = k_loc_start;

    int col_id_start = left_col_id;
    for (int tile_idx = 1; tile_idx < num_tiles + 1; ++tile_idx) {
      int col_id_bound = min(col_id_start + A2_tile, right_col_id + 1);
      while (curr_k_loc < k_loc_bound) {
        int k = A.crds[curr_k_loc];
        if (k >= col_id_bound || k < col_id_start) {
          break;
        }
        ++curr_k_loc;
      }
      table[row_idx - row_idx_start][tile_idx] = curr_k_loc;
      col_id_start += A2_tile;
    }
  }
}