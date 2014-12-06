#ifndef CS207_JB_PARALLEL_HPP
#define CS207_JB_PARALLEL_HPP

#include <iostream>
#include <vector>
#include <cassert>
#include <cmath>
#include <omp.h>
#include <algorithm>

#include <boost/iterator/transform_iterator.hpp>
#include <boost/iterator/filter_iterator.hpp>

#include "CS207/Util.hpp"


namespace jb_parallel {

  // RAII helper class for timing code
  struct Timer {
    std::string msg;
    CS207::Clock clock;
    Timer(const std::string& s) : msg(s) {
      clock.start();
    }
    ~Timer() {
      double elapsed = clock.seconds();
      std::cout << msg << ": " << elapsed << "s" << std::endl;
    }
  };

  // sorts the motherfucker
  void parallel_merge_sort(std::vector<int>& array) {
    int sz = array.size();
    #pragma omp parallel
    {
      int id = omp_get_thread_num();
      int nthreads = omp_get_num_threads();

      auto start = array.begin() + id * sz / nthreads;
      auto end = array.end();
      if (id != nthreads - 1){
        end = array.begin() + (id + 1) * sz / nthreads;
      }
      std::sort(start, end);
    }
  }

  template<class Iter, class UnaryFunction>
  void for_each(Iter first, Iter last, UnaryFunction f) {
    int dist = last - first;
    // Parallelize blockwise
#pragma omp parallel
    {
      int id = omp_get_thread_num();
      int nthreads = omp_get_num_threads();

      auto start = first + id * dist / nthreads;
      auto end = last;
      if (id != nthreads - 1){
        end = first + (id + 1) * dist / nthreads;
      }
      for (auto it = start; it != end; ++it) {
        f(*it);
      }
    }

  }

  // void parallel_transform

  // 



} // namespace jb_parallel

#endif