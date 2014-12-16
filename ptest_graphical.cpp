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
  std::default_random_engine generator(
  std::chrono::system_clock::now().time_since_epoch().count());
  std::uniform_int_distribution<int> distribution(1,12000000);
  auto gen = std::bind(distribution, generator);
  (void) gen;

// for (int N = 10000; N <= 163840000; N *=2){
//   std::vector<int> min_test;
//   for (int i = 0; i < N; ++i) {
//     min_test.push_back(rand() % 20000000);
//   }
//
//
//  	CS207::Clock clock;
//   clock.start();
//    int j = *std::min_element(min_test.begin(), min_test.end());
//  	double elapsed = clock.seconds();
//  	std::cout << N << " , " << elapsed;
//
//
//  	CS207::Clock clock2;
//   clock2.start();
//   auto k = parallel_min(min_test.begin(), min_test.end());
// 	elapsed = clock2.seconds();
// 	std::cout << " , " << elapsed << std::to_string(j) << std::to_string(k) << std::endl;
// }

//
// for (int N = 10000; N <= 163840000; N *=2){
//   std::vector<double> a;
//
//    // Normal loop
//   a = std::vector<double>(N,4);
//  	CS207::Clock clock;
//   clock.start();
//   for (auto it = a.begin(); it < a.end(); ++it)
//   	*it = std::exp(std::sqrt(*it));
// 	double elapsed = clock.seconds();
// 	std::cout << N << " , " << elapsed;
//
//     // Parallel loop
//   a = std::vector<double>(N,4);
//  	CS207::Clock clock2;
//   clock2.start();
// #pragma omp parallel for
//     for (auto it = a.begin(); it < a.end(); ++it) {
//       *it = std::exp(std::sqrt(*it));
//     }
// 		elapsed = clock2.seconds();
// 		std::cout << " , " << elapsed;
//
//
//   // Wrapped in a function
//   a = std::vector<double>(N,4);
//  	CS207::Clock clock3;
//   clock3.start();
//   auto func = [](double ai) { return std::exp(std::sqrt(ai)); };
//   jb_parallel::parallel_transform(a.begin(), a.end(), func);
// 	elapsed = clock3.seconds();
// 	std::cout << " , " << elapsed << std::endl;
//
//
// }
//
//
// for (int N = 10000; N <= 163840000; N *=2){
//   std::vector<int> x;
//
//   for (int i = 0; i < N; ++i) {
//     x.push_back(rand() % 20000000);
//   }
//
//  std::vector<int> y = x;
//
//  	CS207::Clock clock2;
//  	clock2.start();
//   std::sort(y.begin(), y.end());
// 	double elapsed = clock2.seconds();
// 	std::cout << N << " , " << elapsed;
// 	y.clear();
//
//  	CS207::Clock clock;
//  	clock.start();
// 	parallel_sort(x.begin(),x.end());
// 	elapsed = clock.seconds();
// 	std::cout << " , " << elapsed <<std::endl;
// }
//
//
//
//
for (int N = 10000; N <= 163840000; N *=2){
	CS207::Clock clock;
  std::vector<int> z(N, 20);
  int counter = 0;
	std::cout << N;
  auto f = [](int i){if (i % 6 == 0) return i; else return 0;};
	clock.start();
    for (size_t i = 0; i < z.size(); ++i)
      counter += f(z[i]);
	  double elapsed = clock.seconds();
	  std::cout << " , " << elapsed;

  int counter2 = 0;
	CS207::Clock clock2;
		clock2.start();
    jb_parallel::parallel_reduction(z.begin(), z.end(), f, counter2);
    elapsed = clock2.seconds();
    std::cout << " , " << elapsed;

  std::cout << " , " << counter << std::endl;

}


  return 0;
}
