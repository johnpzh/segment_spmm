//
// Created by Zhen Peng on 1/26/2024.
//

#ifndef SPMM_FORMAT_RUNHELPER_H
#define SPMM_FORMAT_RUNHELPER_H

#include <iostream>
#include <algorithm>
#include <vector>
#include <chrono>
#include <math.h>

namespace MT {

template<class Func>
class RunHelper {
private:
  int warmup_ = 3;
  double iter_time_ = 1.0;  /// run benchmark for at least iter_time_ second.
  bool clear_cache_ = false;
  int force_repeat_ = 0;  /// run benchmark only force_repaet times, no warmup.
  int min_repeat_ = 50;  /// minimum times of repeat
  int max_repeat_ = 1000;  /// maximum times of repeat

  std::vector<double> runtimes_;
  double mean_ = 0.0;  /// average runtime
  double std_ = 0.0;  /// standard deviation
  double cv_ = 0.0;  /// coefficient of variation
  Func func_;  /// the benchmark function

  const double cv_threshold_ = 0.05;  ///  coefficient of variation threshold 2%
//  const double cv_threshold_ = 0.02;  ///  coefficient of variation threshold 2%

public:
  explicit RunHelper(Func func) : func_(func) { };
  virtual ~RunHelper() = default;
  void run();

  void set_warmup(int warmup) {
    warmup_ = warmup;
  }

  void set_iter_time(double iter_time) {
    iter_time_ = iter_time;
  }

  void set_clear_cache(bool clear_cache) {
    clear_cache_ = clear_cache;
  }

  void set_force_repeat(int force_repeat) {
    force_repeat_ = force_repeat;
  }

  void set_min_repeat(int min_repeat) {
    min_repeat_ = min_repeat;
  }

  void set_max_repeat(int max_repeat) {
    max_repeat_ = max_repeat;
  }

  void set_func(Func func) {
    func_ = func;
  }

  double get_mean() {
    return mean_;
  }

  double get_std() {
    return std_;
  }

  double get_cv() {
    return cv_;
  }
};

template<class Func>
void RunHelper<Func>::run() {
//  func_();
  int repeat;
  if (force_repeat_) {
    repeat = force_repeat_;
    warmup_ = 0;
  } else {
    /// Warmup
    while (true) {
      printf("\n#### Warm up %d runs in seconds ####\n", warmup_);
      double total_time = 0.0;
      std::vector<double> runtimes;
      for (int i = 0; i < warmup_; ++i) {
        auto start = std::chrono::high_resolution_clock::now();
        func_();
        auto end = std::chrono::high_resolution_clock::now();
        double onetime_sec = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() / 1000000.0;
        printf("%f\n", onetime_sec);
        total_time += onetime_sec;
        runtimes.push_back(onetime_sec);
      }

      repeat = iter_time_ / (total_time / warmup_);
      if (repeat < min_repeat_) {
        repeat = min_repeat_;
      } else if (repeat > max_repeat_) {
        repeat = max_repeat_;
      }

      double mean = total_time / warmup_;
      double std = 0.0;
      for (double t : runtimes) {
        std += (t - mean) * (t - mean);
      }
      std /= warmup_;
      std = sqrt(std);
      double cv = std / mean;
      if (cv < cv_threshold_) {
//      if (cv < 0.10) {
        break;
      } else {
        printf("Coefficient of variation %f%% is too high, warm up again ...\n", cv * 100);
      }
    }
  }

  /// Benchmarking
  printf("\n#### Repeat %d runs in seconds (after %d warm-ups) ####\n", repeat, warmup_);
  runtimes_.clear();

  for (int i = 0; i < repeat; ++i) {
    if (clear_cache_) {
      volatile std::vector<int64_t> dumb(0, 1024 * 1024 * 20);
    }
    auto start = std::chrono::high_resolution_clock::now();
    func_();
    auto end = std::chrono::high_resolution_clock::now();
    double onetime_sec = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() / 1000000.0;
    runtimes_.push_back(onetime_sec);
    printf("%f\n", onetime_sec);
  }

  /// Statistics
  std::sort(runtimes_.begin(), runtimes_.end());

  double mean = 0.0;
  for (double t : runtimes_) {
    mean += t;
  }
  mean /= repeat;
  mean_ = mean;

  double std = 0.0;
  for (double t : runtimes_) {
    std += (t - mean) * (t - mean);
  }
  std /= repeat;
  std = sqrt(std);
  std_ = std;

  cv_ = std / mean;

  if (cv_ > cv_threshold_) {
//  if (cv_ > 0.10) {
    printf("Coefficient of variation %f%% is too high, rerunning...\n", cv_ * 100.0);
    run();  /// Re-run
    return;
  }

  printf("#### mean: %f std: %f cv: %f%% repeat: %d warmup: %d ####\n", mean_, std_, cv_ * 100, repeat, warmup_);
}

} /// namespace MT



#endif //SPMM_FORMAT_RUNHELPER_H
