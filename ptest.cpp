// Performance benchmarks and test for functions in jb_parallel.hpp
// Jonah Kallenbach and Brendan Bozorgmir, CS207 final project 2014.

#include <iostream>
#include <vector>
#include <cassert>
#include <cmath>
#include <omp.h>
#include <algorithm>
#include <random>
#include <chrono>
#include <climits>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include <boost/iterator/transform_iterator.hpp>
#include <boost/iterator/filter_iterator.hpp>

#include "CS207/Util.hpp"
#include "jb_parallel.hpp"

using namespace jb_parallel;

int main () {
  std::vector<int> x;
  std::default_random_engine generator(
      std::chrono::system_clock::now().time_since_epoch().count());
  std::uniform_int_distribution<int> distribution(1,12000000);
  auto gen = std::bind(distribution, generator);
  (void) gen;
  // allocate 100 million random integers
  srand(time(NULL));
  for (int i = 0; i < 30000000; ++i) {
    x.push_back(rand() % 20000000);
  }

  std::vector<int> min_test;
  for (int i = 0; i < 100000400; ++i) {
    min_test.push_back(rand() % 20000000);
  }

  { Timer timer("Serial Min");
    std::cout << "Min element found by serial: " 
              << *(std::min_element(min_test.begin(), min_test.end()))
              << std::endl;
  }

  { Timer timer("jb_parallel parallel_min");
    std::cout << "Min element found by parallel: "
              << parallel_min(min_test.begin(), min_test.end()) << std::endl;
  }

  unsigned N = 40000000;
  std::vector<double> a;

   // Normal loop
  a = std::vector<double>(N,4);
  { Timer timer("Serial Parallel Loop");
    for (auto it = a.begin(); it < a.end(); ++it)
      *it = std::exp(std::sqrt(*it));
  }

    // Parallel loop
  a = std::vector<double>(N,4);
  { Timer timer("OMP Parallel Loop");
#pragma omp parallel for
    for (auto it = a.begin(); it < a.end(); ++it) {
      *it = std::exp(std::sqrt(*it));
    }
  }

  // Wrapped in a function
  a = std::vector<double>(N,4);
  { Timer timer("JB Parallel_Transform");
    // Create a lambda function and call parallel_transform
    auto func = [](double ai) { return std::exp(std::sqrt(ai)); };
    jb_parallel::parallel_transform(a.begin(), a.end(), func);
  }

  std::vector<int> y = x;

  { Timer timer("Serial Sort");
    std::sort(y.begin(), y.end());
  }


  { Timer timer("Parallel Sort");
    parallel_sort(x.begin(),x.end());
  }

  // Check that parallel sort is sorting correctly.
  for (size_t i = 0; i < x.size() - 1; ++i) {
    assert(x[i] <= x[i + 1]);
  }

  return 0;
}