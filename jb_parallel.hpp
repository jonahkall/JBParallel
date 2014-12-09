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
  template<typename IteratorType>
  void parallel_sort(IteratorType begin, IteratorType end) {
    int sz = end - begin;
    int nt = threads_available();
    // currently only works for 4 cores.
    assert(nt == 4);
    std::vector<IteratorType> bounds(nt*2);
    #pragma omp parallel
    {
      int id = omp_get_thread_num();
      int nthreads = omp_get_num_threads();
      //std::cout << nthreads << std::endl;

      auto start = begin + id * sz / nthreads;
      //std::cout << "start bound: " << id * sz / nthreads << std::endl;
      auto finish = end;
			if (id != nthreads - 1){
        finish = begin + (id + 1) * sz / nthreads;
      }
      bounds[2*id] = start;
      bounds[2*id + 1] = finish;
      std::sort(start, finish);
    }
    std::sort(bounds.begin(), bounds.end());
    // Time to do some in_place merges
    // currently hardcoded for 4 cores
    std::inplace_merge(bounds[0],bounds[2], bounds[3]);
    std::inplace_merge(bounds[4],bounds[6], bounds[7]);
    std::inplace_merge(bounds[0], bounds[4], bounds[7]);
  }

	template <typename IteratorType>
	// function for finding minimum of standard vector
	// in O(n/p). Works for any number of cores.
	// Requires a random access iterator
	typename std::iterator_traits<IteratorType>::value_type
	parallel_min(IteratorType begin, IteratorType end) {
    int nt = threads_available();
    std::vector<int> mins(nt);
    auto sz = end - begin;
    #pragma omp parallel
    {
      int id = omp_get_thread_num();
      int nthreads = omp_get_num_threads();
      auto start = id * sz / nthreads;

      int finish = sz;
      if (id != nthreads - 1){
        finish = (id + 1) * sz / nthreads;
      }
			auto min_seen = *(begin + start);
      for (auto i = begin + start; i <  begin + finish; ++i) {
        min_seen = std::min(min_seen, *i);
      }
      mins[id] = min_seen;
    }
		auto min_total = mins[0];
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
