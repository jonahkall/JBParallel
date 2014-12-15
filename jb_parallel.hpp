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

/** @file jb_parallel.hpp
 * @brief Define the key functions built on top of
 * OpenMP to be used to optimize code which iterates
 * over a random access data structure.
 */

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
    // This compiler directive instructs OpenMP
    // to make this block of code parallel.
    #pragma omp parallel
    {
      p = omp_get_num_threads();
    }
    return p;
  }

  /** A function that sorts a range in parallel
   * @param begin, end the range of what is to be sored
	 * @tparam IteratorType must meet the requirements of RandomAccessIterator
   * @pre begin < end and begin and end are in the same range
	 * @pre IteratorType::value_type supports <, + int (looser than random access!) 
	 * @pre Machine has exactly 1, 2, or 4 available cores
   * @post For all i in [begin, end] !(i+1 < i)
   * Runtime O((nlog(n))/p)
	 * std::sort used in this function uses a quicksort/heapsort hybrid known as intro sort
   */
  template<typename IteratorType>
  void parallel_sort(IteratorType begin, IteratorType end) {
    int sz = end - begin;
    int nt = threads_available();
    if (nt == 1) {
      std::sort(begin, end);
      return;
    }
    std::vector<IteratorType> bounds(nt*2);
    #pragma omp parallel
    {
      int id = omp_get_thread_num();
      int nthreads = omp_get_num_threads();

      auto start = begin + id * sz / nthreads;
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
    // This section hardcoded for four cores.
    if (nt == 4) {
      std::inplace_merge(bounds[0],bounds[2], bounds[3]);
      std::inplace_merge(bounds[4],bounds[6], bounds[7]);
      std::inplace_merge(bounds[0], bounds[4], bounds[7]);
    }

    // This sectino hardcoded for 2 cores: only a single inplace
    // merge necessary.
    if (nt == 2) {
      std::inplace_merge(bounds[0],bounds[2], bounds[3]);
    }
  }

 /** A function that finds the minimum in a range in parallel
  * @param begin, end the range of what is to be sored
  * @tparam IteratorType must meet the requirements of RandomAccessIterator
  * @pre begin < end and begin and end are in the same range
  * @pre IteratorType::value_type supports the < function
  * @post For all i in [begin, end] !(i+1 < i)
  * Runtime O((nlog(n))/4  + 2n)
  * std::sort used in this function uses a quicksort/heapsort hybrid known as intro sort
  */

	template <typename IteratorType>

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

  /** A function that applies 
   * @param first, last the range to apply the function to
	 * @param f the unary function object to be applied
	 * @tparam Iter must meet the requirements of RandomAccessIterator
	 * @tparam UnaryFunction must meet the requirements of MoveConstructible
   * @pre first < last and first and last are in same range
	 * @pre 
   * @post for all i in [first,last] *(new first) = f(* old first)
   * Runtime O(O(f)*(last - first)/nthreads)
  */
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

 /** A function that applies 
  * @param first, last the range to apply the function to
 * @param f the unary function object to be applied
 * @tparam Iter must meet the requirements of RandomAccessIterator
 * @tparam UnaryFunction must meet the requirements of MoveConstructible
  * @pre first < last and first and last are in same range
  * @return for all i in [first,last] return[i]= f[range[i]]
  * Runtime O(O(f)*(last - first)/nthreads)
 */
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

  /** BRIEF DESCRIPTION
   * @param NAME DESCRIPTION
   * @return DESCRIPTION OF RETURN VALUE
   * @pre PRECONDITION
   * @post POSTCONDITION
   *
   * LONGER DESCRIPTION, RUNTIME, EXAMPLES, ETC 
  */
  template<class Iter, class UnaryFunction, class T>
  void parallel_reduction(Iter first, Iter last, UnaryFunction f, T& counter) {
    int dist = last - first;
    int nt = threads_available();
    std::vector<T> counts(nt, 0);
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
        counts[id] += f(*it);
      }
    }
    for (int i = 0; i < nt; ++i) {
      counter += counts[i];
    }
  }

} // namespace jb_parallel

#endif
