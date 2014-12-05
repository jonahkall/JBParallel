#ifndef CS207_JB_PARALLEL_HPP
#define CS207_JB_PARALLEL_HPP

#include <iostream>
#include <vector>
#include <cassert>
#include <cmath>
#include <omp.h>

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
    (void) array;
  }

  template<class Iter, class UnaryFunction>
  void for_each(Iter first, Iter last, UnaryFunction f) {
#pragma omp parallel for [cyclic, block, dynamic, guided]
    for (auto it = first; it != last; ++it) {
      f(*it);
    }
  }

}

/* 

#pragma omp parallel
{
  int id = ?
  int 
}

*/

#endif