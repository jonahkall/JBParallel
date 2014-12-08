#ifndef CS207_JB_PARALLEL_HPP
#define CS207_JB_PARALLEL_HPP

#include <iostream>
#include <vector>
#include <cassert>
#include <cmath>
#include <omp.h>
#include <algorithm>
#include <climits>

#include <boost/iterator/transform_iterator.hpp>
#include <boost/iterator/filter_iterator.hpp>

#include "CS207/Util.hpp"


namespace jb_parallel {

  // RAII helper class for timing code.
  // Cris wrote this.
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

  // Simple function which tells the user how many
  // threads OMP is making available.
  int threads_available() {
    int p;
    #pragma omp parallel
    {
      p = omp_get_num_threads();
    }
    return p;
  }

  // Sorts an array of ints in O(n log n /p)
  void parallel_sort(std::vector<int>& array) {
    int sz = array.size();
    int nt = threads_available();
    // currently only works for 4 cores.
    assert(nt == 4);
    std::vector<std::vector<int>::iterator> bounds(nt*2);
    #pragma omp parallel
    {
      int id = omp_get_thread_num();
      int nthreads = omp_get_num_threads();
      //std::cout << nthreads << std::endl;

      auto start = array.begin() + id * sz / nthreads;
      //std::cout << "start bound: " << id * sz / nthreads << std::endl;
      auto end = array.end();
      if (id != nthreads - 1){
        end = array.begin() + (id + 1) * sz / nthreads;
      }
      bounds[2*id] = start;
      bounds[2*id + 1] = end;
      std::sort(start, end);
    }
    std::sort(bounds.begin(), bounds.end());
    // Time to do some in_place merges
    // currently hardcoded for 4 cores
    std::inplace_merge(bounds[0],bounds[2], bounds[3]);
    std::inplace_merge(bounds[4],bounds[6], bounds[7]);
    std::inplace_merge(bounds[0], bounds[4], bounds[7]);
  }

  template<typename T>
  T parallel_min(std::vector<T>& x) {
    int nt = threads_available();
    std::vector<int> mins(nt);
    auto sz = x.size();
    #pragma omp parallel
    {
      int id = omp_get_thread_num();
      int nthreads = omp_get_num_threads();
      auto start = id * sz / nthreads;

      int end = sz;
      if (id != nthreads - 1){
        end = (id + 1) * sz / nthreads;
      }
      int min_seen = LONG_MAX;
      for (int i = start; i < end; ++i) {
        min_seen = std::min(min_seen, x[i]);
      }
      mins[id] = min_seen;
    }
    int min_total = LONG_MAX;
    for (int i = 0; i < mins.size(); ++i) {
      min_total = std::min(min_total, mins[i]);
    }
    return min_total;
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

  template<class Iter, class UnaryFunction>
  void parallel_transform(Iter first, Iter last, UnaryFunction f) {
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
        *it = f(*it);
      }
    }
  }

} // namespace jb_parallel

#endif
