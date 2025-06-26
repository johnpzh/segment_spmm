//
// Created by Zhen Peng on 6/25/25.
//

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <chrono>
#include <vector>
#include <string>
#include "mt_io.h"
#include "utils.h"
#include "spmm_csr.h"

int main()
{
  const char *filename = "/Users/peng599/pppp/comet-amais-memory/COMET-memAnalysis2-dev-new/test/integration/data/test_rank2.mtx";
  MT::DRAMCSRMatrix csr_matrix;
  MT::read_csr_matrix(filename, &csr_matrix);
  int64_t A1 = csr_matrix.num_rows_;
  int64_t A2 = csr_matrix.num_cols_;

  uint64_t num_tiles;
  std::vector<std::vector<uint64_t>> table;
  MT::get_tiling_table_for_pos(
      &csr_matrix,
      /*row_idx_start=*/0,
      /*row_idx_bound=*/5,
      /*A2_tile=*/2,
      num_tiles/*out*/,
      table/*out*/);

  /// print
  MT::print_csr_matrix(&csr_matrix);
  printf("\nprint the table...\n");
  for (auto &row : table) {
    for (auto e : row) {
      std::cout << e << " ";
    }
    std::cout << "\n";
  }
  std::cout << "\n";

  return EXIT_SUCCESS;
}