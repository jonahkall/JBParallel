/**
 * @file mtl_test.cpp
 * Test script for interfacing with MTL4 and it's linear solvers.
 */

// HW3: Need to install/include Boost and MTL in Makefile
#include <boost/numeric/mtl/mtl.hpp>
#include <boost/numeric/itl/itl.hpp>
#include <iostream>

using namespace itl;
using namespace mtl;

// This class provides an MTL-compatible identity matrix class (nicely,
// it doesn't have to store anything, because to multiply by a vector,
// it just spits back the vector itself)
class IdentityMatrix {
private:
  int dim_;
public:
  IdentityMatrix(int d, int) : dim_(d) {};

  int dim () const {
    return dim_;
  }

  /** Helper function to perform multiplication. Allows for delayed
  * evaluation of results.
  * Assign::apply(a, b) resolves to an assignment operation such as 
  * a += b, a -= b, or a = b.
  * @pre @a size(v) == size(w) */
  template <typename VectorIn , typename VectorOut , typename Assign>
      void mult(const VectorIn& v, VectorOut& w, Assign) const {
    assert(size(v) == size(w));
    for (auto i = 0; i < (int)size(w); ++i) {
      w[i] = Assign::apply(w[i], v[i]);
    }
  }

  template <typename Vector>
  mtl::vec::mat_cvec_multiplier<IdentityMatrix, Vector>
      operator*(const Vector& v) const {
    return mtl::vec::mat_cvec_multiplier <IdentityMatrix, Vector>(*this, v);
  }
};

/** Traits that MTL uses to determine properties of our IdentityMatrix. */
namespace mtl {
namespace ashape {

 /** Define IdentityMatrix to be a non-scalar type. */
template <>
struct ashape_aux <IdentityMatrix > {
   typedef nonscal type;
};
} // end namespace ashape

/** IdentityMatrix implements the Collection concept
 * with value_type and size_type */
template <>
struct Collection <IdentityMatrix > {
typedef double value_type;
typedef unsigned size_type;
};
} // end namespace mtl

inline std::size_t size(const IdentityMatrix& A) {
  return A.dim() * A.dim();
}

inline std::size_t num_rows(const IdentityMatrix& A) {
  return A.dim();
}

inline std::size_t num_cols(const IdentityMatrix& A) {
  return A.dim();
}

int main()
{
  const int size = 40, N = size * size;
  typedef IdentityMatrix matrix_type;
  matrix_type                   A(N, N);

  // Create an ILU(0) preconditioner
  //pc::ilu_0<matrix_type>        P(A);
  pc::identity<matrix_type> P(A);
  // Set b such that x == 1 is solution; start with x == 0
  dense_vector<double>          x(N, 1.0), b(N);
  b= A * x; x= 0;
  
  // Termination criterion: r < 1e-6 * b or N iterations
  noisy_iteration<double>       iter(b, 500, 1.e-6);
  
  // Solve Ax == b with left preconditioner P
  cg(A, x, b, P, iter);

  assert(x == b);

  return 0;
}
