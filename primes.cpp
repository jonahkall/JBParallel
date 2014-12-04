#include "CS207/Util.hpp"
#include <cmath>
#include <set>

/** Return true iff @a n is prime.
 * @pre @a n >= 0
 */
bool is_prime(int n)
{
  // Small speed up: if the number is even and not 2, it's not prime.
  if (n % 2 == 0 && n != 2)
    return false;
  // Small speed up: fill the primes container with some small primes.
  static std::set<int> primes_set = {2,3,5,7,11};
  std::cout << *(--primes_set.end()) << "\n";
  assert(n >= 0);
  // If is_prime hasn't been called yet with anything bigger than 121,
  // call it on all smaller input. This slows down a single call to 
  // is_prime but helps with memoization in the future.
  if (primes_set.size() == 5) {
    for (int i = 2; i < sqrt(n); ++i)
	if (is_prime(i)) {
          primes_set.insert(i);
        }
  }
  // Iterate over the set of primes <= sqrt(n). If any of these primes divide
  // n, it is not prime.
  for (auto it = primes_set.begin();
      it != primes_set.end(); ++it) {
    if (*it > sqrt(n)) {
      break;
    }
    if (n % *it == 0)
      return false;
  }
  // If the set of primes is too small, we may not have seen all the primes
  // up to sqrt(n), so check all numbers from max(primes_set) to sqrt(n).
  for (int i = *(--primes_set.end()); i < sqrt(n); ++i) {
    if (n % i == 0)
      return false;
  }
  // Cache the prime we found so that future calls to is_prime can use it.
  primes_set.insert(n);
  return true;
}

int main()
{
  while (!std::cin.eof()) {
    // How many primes to test? And should we print them?
    std::cerr << "Input Number: ";
    int n = 0;
    CS207::getline_parsed(std::cin, n);
    if (n <= 0)
      break;
    std::cout << is_prime(23493495);
    std::cerr << "Print Primes (y/n): ";
    char confirm = 'n';
    CS207::getline_parsed(std::cin, confirm);
    bool print_primes = (confirm == 'y' || confirm == 'Y');

    CS207::Clock timer;
    // Loop and count primes from 2 up to n
    int num_primes = 0;
    for (int i = 2; i <= n; ++i) {
      if (is_prime(i)) {
        ++num_primes;
        if (print_primes)
          std::cout << i << std::endl;
      }
    }

    double elapsed_time = timer.seconds();

    std::cout << "There are " << num_primes
              << " primes less than or equal to " << n << ".\n"
              << "Found in " << (1000 * elapsed_time) << " milliseconds.\n\n";
  }

  return 0;
}
