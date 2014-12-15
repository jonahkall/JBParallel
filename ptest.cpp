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

#include "thrust/sort.h"
#include "thrust/for_each.h"
#include "thrust/reduce.h"

#include "CS207/Util.hpp"
#include "jb_parallel.hpp"

using namespace jb_parallel;

int main () {
  std::default_random_engine generator(
      std::chrono::system_clock::now().time_since_epoch().count());
  std::uniform_int_distribution<int> distribution(1,12000000);
  auto gen = std::bind(distribution, generator);
  (void) gen;
  // allocate 100 million random integers
  // srand(time(NULL));
  // for (int i = 0; i < 30000000; ++i) {
  //   x.push_back(rand() % 20000000);
  // }
//
// for (int j = 1; j < pow(10,10); j *= 10){
//   std::vector<int> min_test;
//   for (int i = 0; i < j; ++i) {
//     min_test.push_back(rand() % 20000000);
//   }
//
//
//   { Timer timer("Serial Min " + std::to_string(j));
//   	int j = *std::min_element(min_test.begin(), min_test.end());
// 		std::cout << std::to_string(j) << std::endl;
//   }
//
//   { Timer timer("jb_parallel parallel_min " + std::to_string(j));
//      auto j = parallel_min(min_test.begin(), min_test.end());
// 		 std::cout << std::to_string(j) << std::endl;
//   }
//
// }
//
// 	for (int N = 1; N < pow(10,9); N *= 10){
//   std::vector<double> a;
//
//    // Normal loop
//   a = std::vector<double>(N,4);
//   { Timer timer("Serial Parallel Loop " + std::to_string(N));
//     for (auto it = a.begin(); it < a.end(); ++it)
//       *it = std::exp(std::sqrt(*it));
//   }
//
//     // Parallel loop
//   a = std::vector<double>(N,4);
//   { Timer timer("OMP Parallel Loop " + std::to_string(N));
// #pragma omp parallel for
//     for (auto it = a.begin(); it < a.end(); ++it) {
//       *it = std::exp(std::sqrt(*it));
//     }
//   }
//
//   // Wrapped in a function
//   a = std::vector<double>(N,4);
//   { Timer timer("JB Parallel_Transform " + std::to_string(N));
//     // Create a lambda function and call parallel_transform
//     auto func = [](double ai) { return std::exp(std::sqrt(ai)); };
//     jb_parallel::parallel_transform(a.begin(), a.end(), func);
//   }
// }


// for (int N = 1; N < pow(10,9); N *= 10){
//   std::vector<int> x;
//
//   for (int i = 0; i < N; ++i) {
//     x.push_back(rand() % 20000000);
//   }
//
// 	std::vector<int> y = x;
//
//   { Timer timer("Serial Sort " + std::to_string(N) );
//     std::sort(y.begin(), y.end());
//   }
//
//   { Timer timer("Parallel Sort " + std::to_string(N));
//     parallel_sort(x.begin(),x.end());
//   }
// }
//
//   // Check that parallel sort is sorting correctly.
//   for (size_t i = 0; i < x.size() - 1; ++i) {
//     assert(x[i] <= x[i + 1]);
//   }
//

for (int N = 1; N < pow(10,9); N *=10){

  std::vector<int> z(N, 2);
  int counter = 0;
  {Timer timer("Serial reduction " + std::to_string(N));
    auto f = [](int i){if (i % 2 == 0) return i; else return 0;};
    for (size_t i = 0; i < z.size(); ++i) {
      counter += f(z[i]);
    }
  }

  std::cout << counter << std::endl;


  counter = 0;
  {Timer timer("Parallel reduction " + std::to_string(N));
    auto f = [](int i){if (i % 2 == 0) return i; else return 0;};
    jb_parallel::parallel_reduction(z.begin(), z.end(), f, counter);
  }

  std::cout << counter << std::endl;
}


  return 0;
}