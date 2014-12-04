/**
 * @file poisson.cpp
 * Test script for treating the Graph as a MTL Matrix
 * and solving a Poisson equation.
 *
 * @brief Reads in two files specified on the command line.
 * First file: 3D Points (one per line) defined by three doubles.
 * Second file: Eges (one per line) defined by 2 indices into the point list
 *              of the first file.
 *
 * Launches an SDLViewer to visualize the solution.
 */

#include <fstream>
#include <math.h>

#include "CS207/SDLViewer.hpp"
#include "CS207/Util.hpp"

#include <boost/numeric/mtl/mtl.hpp>
#include <boost/numeric/itl/itl.hpp>

#include "Graph.hpp"
#include "Point.hpp"
#include "BoundingBox.hpp"

using namespace itl;
using namespace mtl;

// Node value struct stores a tag if the node is on
// the boundary.
struct node_value_type {
  // true if it is on the boundary.
  bool tag;
};

// We don't need to store anything in edge values.
struct edge_value_type {
};

typedef Graph<node_value_type, edge_value_type> GraphType;
typedef int size_type;
typedef std::map<typename GraphType::node_type, unsigned> NodeMapType;

/** Remove all the nodes in graph @a g whose position is contained within
 * BoundingBox @a bb
 * @post For all i, 0 <= i < @a g.num_nodes(),
 *        not bb.contains(g.node(i).position())
 * Note: I did a visual check to see that my current method works!
 */
void remove_box(GraphType& g, const BoundingBox& bb) {
  // HW3: YOUR CODE HERE
  for (auto it = g.node_begin(); it != g.node_end(); ++it) {
      if (bb.contains((*it).position())) {
      for (auto e_it = (*it).edge_begin(); e_it != (*it).edge_end(); ++e_it) {
        (*e_it).node2().value().tag = true;
      }
      g.remove_node((*it));
      --it;
    }
  }
  return;
}

//This is the forcing function for the Poisson equation
double f(const Point& p) {
  return 5.0 * cos(norm_1(p));
}


// g returns the boundary conditions for the nodes:
// if g did not return -1, we are on a boundary.
double g(const Point& p) {
  if (norm_inf(p) == 1) {
    return 0.0;
  }
  else if (norm_inf(p-Point(0.6,0.6,0)) < 0.2) {
    return -0.2;
  }
  else if (norm_inf(p-Point(-0.6,-0.6,0)) < 0.2) {
    return -0.2;
  }
  else if (norm_inf(p-Point(-0.6,0.6,0)) < 0.2) {
    return -0.2;
  }
  else if (norm_inf(p-Point(0.6,-0.6,0)) < 0.2) {
    return -0.2;
  }
  else {
    BoundingBox bb = BoundingBox(Point(-0.6,-0.2,-1), Point(0.6,0.2,1));
    if (bb.contains(p)) {
      return 1.0;
    }
  }
  return -1.0;
}

// This class wraps a Graph and acts as a matrix which is MTL-compatible.
// (weak) RI: dim_ = g_->num_nodes()
//     g_ satisfies the graph concept and has boundary nodes tagged correctly.
struct GraphSymmetricMatrix {
  GraphType* g_;
  size_type dim_;

  size_type dim () const {
    return dim_;
  }

  /** Helper function to perform multiplication.
  * Assign::apply(a, b) resolves to an assignment operation
  * @pre @a size(v) == size(w) 
  * @pre all nodes on the boundaries are correctly tagged.
  */
  template <typename VectorIn, typename VectorOut, typename Assign>
      void mult(const VectorIn& v, VectorOut& w, Assign) const {
    assert(size(v) == size(w));
    for (auto it = g_->node_begin(); it != g_-> node_end(); ++it) {
      GraphType::Node n = *it;
      auto i = n.index();
      double val = 0.0;
      if (n.value().tag) {
        val += v[i];
      }
      else {
        val += -1.0*((double)n.degree()) * v[i];
        for (auto it2 = (*it).edge_begin(); it2 != (*it).edge_end(); ++it2) {
          if (!(*it2).node2().value().tag) {
            val += (double)v[(*it2).node2().index()];
          }
        }
      }
      Assign::apply(w[i], val);
    }
  }

  template <typename Vector>
  mtl::vec::mat_cvec_multiplier<GraphSymmetricMatrix, Vector>
      operator*(const Vector& v) const {
    return mtl::vec::mat_cvec_multiplier
        <GraphSymmetricMatrix, Vector>(*this, v);
  }

  // GraphSymmetricMatrix constructor which takes a pointer to a graph
  // and initializes the graph pointer to that pointer and the dimension
  // to the number of nodes in the graph.
  GraphSymmetricMatrix(GraphType* g) {
    g_ = g;
    dim_ = g->num_nodes();
  };
};

struct ColorFunctor {
  template <typename NODE>
    CS207::Color operator()(const NODE& n) {
      if (n.value().tag) {
        return CS207::Color::make_heat(0.8);
      }
      else {
        return CS207::Color::make_heat(0.3);
      }
    }
};

/** Traits that MTL uses to determine properties of our IdentityMatrix. */
namespace mtl {
namespace ashape {

 /** Define GraphSymmetricMatrix to be a non-scalar type. */
template <>
struct ashape_aux <GraphSymmetricMatrix> {
   typedef nonscal type;
};
} // end namespace ashape

/** GraphSymmetricMatrix implements the Collection concept
 * with value_type and size_type */
template <>
struct Collection <GraphSymmetricMatrix> {
typedef double value_type;
typedef unsigned size_type;
};
} // end namespace mtl

inline std::size_t size(const GraphSymmetricMatrix& A) {
  return A.dim() * A.dim();
}

inline std::size_t num_rows(const GraphSymmetricMatrix& A) {
  return A.dim();
}

inline std::size_t num_cols(const GraphSymmetricMatrix& A) {
  return A.dim();
}

// Functor which doesn't change x and y coordinates, but for a node n, 
// sets z value to be u[n.index()], to visualize convergence and solution
// of the poisson.
template<typename Vector>
struct NodePosition {
  Vector u_;
  NodePosition(Vector u) : u_(u) {};
  template <typename NODE>
    Point operator()(const NODE& node) {
      return Point(node.position().x, node.position().y, u_[node.index()]);
  }
};

// Functor which produces a color from a node and a vector with its min and
// max values specified.
template<typename Vector>
struct NodeColor {
  double max_;
  double min_;
  Vector u_;
  NodeColor(double max, double min, Vector u) : max_(max), min_(min), u_(u) {};
  template <typename NODE>
    CS207::Color operator()(const NODE& n) {
      // We don't want a division by zero error, so we check for that case.
      if (max_ == min_) {
        return CS207::Color::make_heat(0);
      }
      return CS207::Color::make_heat((u_[n.index()] - min_)/(max_ - min_));
    }
};

// visual_iteration is a class which inherits from and extends the iteration callback
// cyclic_iteration by overriding its "finished" methods which check for convergence
// to also visualize the current solution with SDLViewer.
template<class Vector, class Real>
class visual_iteration : public cyclic_iteration<Real> {

typedef cyclic_iteration<Real> super;

public:

  visual_iteration(CS207::SDLViewer* viewer, GraphType* g,
      dense_vector<Real>* u, const Vector& r0, int max_iter_, double tol_,
      double atol_ = Real(0)) : super(r0, max_iter_, tol_, atol_),
      viewer_(viewer), graph(g), u_(u) {
    node_map = viewer->empty_node_map(*graph);
    i_ = 0;
  }

  void draw () {
    auto node_map = viewer_->empty_node_map(*graph);
    viewer_->clear();
    // Update viewer with nodes' new positions
    Real min_val = (*u_)[0];
    Real max_val = (*u_)[0];
    for (int i = 1; i < (int)size(*u_); ++i){
      if((*u_)[i] < min_val) {
        min_val = (*u_)[i];
      }
      if ((*u_)[i] > max_val) {
        max_val = (*u_)[i];
      }
    }
    NodeColor<Vector> cf(max_val, min_val, *u_);
    NodePosition<Vector> np(*u_);
    viewer_->add_nodes(graph->node_begin(), graph->node_end(), cf, np, node_map);
    viewer_->add_edges(graph->edge_begin(), graph->edge_end(), node_map);
    viewer_->center_view();
    // Adjust this parameter to change how fast the solution loads. Sleeping for
    // some amount of time makes it easier to visualize convergence.
    CS207::sleep(0.24);
  }

  void print_resid() {
    //std::cout << "iteration " << i_ << ": resid " << super::resid() << "\n";
  }

  bool finished () {
    ++i_;
    draw();
    print_resid();
    return super::finished();
  }

  template <typename T>
  bool finished(const T& r) {
    ++i_;
    bool ret = super::finished(r);
    draw();
    print_resid();
    return ret;
  }

private:
  friend class cyclic_iteration<Real>;
  CS207::SDLViewer* viewer_;
  NodeMapType node_map;
  GraphType* graph;
  Vector* u_;
  int i_;
};

int main(int argc, char** argv)
{
  // Check arguments
  if (argc < 2) {
    std::cerr << "Usage: " << argv[0] << " NODES_FILE EDGES_FILE\n";
    exit(1);
  }

  // Define an empty Graph
  GraphType graph;

  // Create a nodes_file from the first input argument
  std::ifstream nodes_file(argv[1]);
  // Interpret each line of the nodes_file as a 3D Point and add to the Graph
  std::vector<typename GraphType::node_type> node_vec;
  Point p;
  while (CS207::getline_parsed(nodes_file, p))
    node_vec.push_back(graph.add_node(2*p - Point(1,1,0)));

  // Create a tets_file from the second input argument
  std::ifstream tets_file(argv[2]);
  // Interpret each line of the tets_file as four ints which refer to nodes
  std::array<int,4> t;
  while (CS207::getline_parsed(tets_file, t)) {
    graph.add_edge(node_vec[t[0]], node_vec[t[1]]);
    graph.add_edge(node_vec[t[0]], node_vec[t[2]]);
    graph.add_edge(node_vec[t[1]], node_vec[t[3]]);
    graph.add_edge(node_vec[t[2]], node_vec[t[3]]);
  }

  for (auto it = graph.node_begin(); it != graph.node_end(); ++it) {
    (*it).value().tag = false;
  }

  // Get the edge length, should be the same for each edge
  double h = graph.edge(0).length();
  // Make holes in our Graph
  remove_box(graph, BoundingBox(Point(-0.8+h,-0.8+h,-1), Point(-0.4-h,-0.4-h,1)));
  remove_box(graph, BoundingBox(Point( 0.4+h,-0.8+h,-1), Point( 0.8-h,-0.4-h,1)));
  remove_box(graph, BoundingBox(Point(-0.8+h, 0.4+h,-1), Point(-0.4-h, 0.8-h,1)));
  remove_box(graph, BoundingBox(Point( 0.4+h, 0.4+h,-1), Point( 0.8-h, 0.8-h,1)));
  remove_box(graph, BoundingBox(Point(-0.6+h,-0.2+h,-1), Point( 0.6-h, 0.2-h,1)));

  // Make sure nodes are default not on the boundary.
  for (auto it = graph.node_begin(); it != graph.node_end(); ++it) {
    if (g((*it).position()) != -1) {
      (*it).value().tag = true;
    }
  }

  int N = graph.num_nodes();
  typedef GraphSymmetricMatrix matrix_type;
  dense_vector<double> b(N, 0.0);

  // Initialize b as specified on the pset.
  for (auto it = graph.node_begin(); it != graph.node_end(); ++it) {  
    if (g((*it).position()) != -1) {
      b[(*it).index()] = g((*it).position());
    }
    else {
      double adj_sum = pow(h,2) * f((*it).position());
      for (auto inc_it = (*it).edge_begin(); inc_it != (*it).edge_end(); ++inc_it) {
        if (g((*inc_it).node2().position()) != -1) {
          adj_sum -= g((*inc_it).node2().position());
        }
      }
      b[(*it).index()] = adj_sum;
    }
  }

  matrix_type gsm = GraphSymmetricMatrix(&graph);

  dense_vector<double> u(N, 0.0);
  u = 0;

  CS207::SDLViewer viewer;
  viewer.launch();
  visual_iteration<dense_vector<double>, double> iter(&viewer, &graph, &u, b, 500, 1.e-10);

  // Solve Au == b
  cg(gsm, u, b, iter);

  return 0;
}
