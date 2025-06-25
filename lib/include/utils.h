//
// Created by Zhen Peng on 2/15/2024.
//

#ifndef SPMM_FORMAT_UTILS_H
#define SPMM_FORMAT_UTILS_H

#include <vector>
#include <fstream>

const uint64_t MAX_B2 = 4096;
const uint64_t MAX_MEMORY_GB = 16;  /// NVIDIA P100 GPU has 16 GB memory.

template<class T>
double calculate_average_in_vector(std::vector<T> exe_times)
{
  if (exe_times.empty()) {
    return 0;
  }

  T avg_time = 0;
  for (auto time : exe_times) {
    avg_time += time;
  }
  avg_time /= exe_times.size();

  return avg_time;
}

template<class T>
void save_exe_times_into_file(std::string filename, std::vector<T> exe_times)
{
  std::ofstream fout;
  fout.open(filename);
  if (fout.is_open()) {
    for (auto time : exe_times) {
      fout << time << "\n";
    }
  } else {
    std::cerr << "Error: cannot create file " << filename << "\n";
  }
}

#endif //SPMM_FORMAT_UTILS_H
