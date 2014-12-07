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
  std::default_random_engine generator(std::chrono::system_clock::now().time_since_epoch().count());
  std::uniform_int_distribution<int> distribution(1,12000000);
  auto gen = std::bind(distribution, generator);
  (void)gen;
  // allocate 100 million random integers
  srand(time(NULL));
  for (int i = 0; i < 110000000; ++i) {
    x.push_back(rand() % 200000000);
  }

  /*for (int i = 0; i < 20; ++i) {
    std::cout << x[i] << " ";
  }
  std::cout << std::endl << std::endl << std::endl;

  parallel_sort(x);*/

  { Timer timer("Serial Min");
    std::cout << *(std::min_element(x.begin(), x.end())) << std::endl;
  }

  { Timer timer("jb_parallel parallel_min");
    std::cout << parallel_min(x) << std::endl;
  }

  return 0;
}