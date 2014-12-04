/**
 * @file subgraph.cpp
 * Test script for viewing a subgraph from our Graph
 *
 * @brief Reads in two files specified on the command line.
 * First file: 3D Points (one per line) defined by three doubles
 * Second file: Tetrahedra (one per line) defined by 4 indices into the point
 * list
 */

#include <fstream>
#include <iterator>

#include "CS207/SDLViewer.hpp"
#include "CS207/Util.hpp"

#include "Graph.hpp"

/** An iterator that skips over elements of another iterator based on whether
 * those elements satisfy a predicate.
 *
 * Given an iterator range [@a first, @a last) and a predicate @a pred,
 * this iterator models a filtered range such that all i with
 * @a first <= i < @a last and @a pred(*i) appear in order of the original range.
 * RI: (it_ == end) || p_(*it)
 */
template <typename Pred, typename It>
class filter_iterator
    : private equality_comparable<filter_iterator<Pred,It>> {
 public:
  // Get all of the iterator traits and make them our own
  typedef typename std::iterator_traits<It>::value_type        value_type;
  typedef typename std::iterator_traits<It>::pointer           pointer;
  typedef typename std::iterator_traits<It>::reference         reference;
  typedef typename std::iterator_traits<It>::difference_type   difference_type;
  typedef typename std::input_iterator_tag                     iterator_category;

  typedef filter_iterator<Pred,It> self_type;

  // Constructor for filter_iterator: makes sure that iterator starts somewhere
  // that the predicate is satisfied.
  filter_iterator(const Pred& p, const It& first, const It& last)
      : p_(p), it_(first), end_(last) {
    if (first == last)
      return;
    while (!p_(*it_)) {
      ++it_;
      if (first == last)
        return;
    }
  }

  // Dereferences by dereferencing It it_.
  // Iterator must be to a valid value (RI).
  value_type operator*() const {
    return *it_;
  }

  // Increment operator.  If possible, increments 
  self_type& operator++() {
    if (it_ == end_)
      return *this;
    ++it_;
    if (it_ == end_) {
      return *this;
    }
    while (it_ != end_){
      if (!p_(*it_)) {
        ++it_;
      }
      else {
        return *this;
      }
    }
    return *this;
  }

  /** Compare two filter iterators for equality by checking their current 
   * It it_ and end_ values are the same.
   * @pre: fi.p_ is the same as p_, otherwise behavior is undefined.
   */
  bool operator==(const self_type& fi) const {
    return (fi.it_ == it_ && fi.end_ == end_);
  }

 private:
  Pred p_;
  It it_;
  It end_;
};

/** Helper function for constructing filter_iterators.
 *
 * Usage:
 * // Construct an iterator that filters odd values out and keeps even values.
 * std::vector<int> a = ...;
 * auto it = make_filtered(a.begin(), a.end(), [](int k) {return k % 2 == 0;});
 */
template <typename Pred, typename Iter>
filter_iterator<Pred,Iter> make_filtered(const Iter& it, const Iter& end,
                                         const Pred& p) {
  return filter_iterator<Pred,Iter>(p, it, end);
}

// This predicate makes the graph more sparse by getting rid of about 2/3
// of the nodes, namely any nodes which have indices not 1 mod 3.  This
// maintains local connectivity in strips of 2, leading to an interesting
// visual pattern.
struct SparsePredicate {
  template <typename NODE>
  bool operator()(const NODE& n) {
    return (n.index() % 3 == 1);
  }
};


int main(int argc, char** argv)
{
  // Check arguments
  if (argc < 3) {
    std::cerr << "Usage: " << argv[0] << " NODES_FILE TETS_FILE\n";
    exit(1);
  }

  // Construct a Graph
  typedef Graph<int, int> GraphType;
  GraphType graph;
  std::vector<GraphType::node_type> nodes;

  // Create a nodes_file from the first input argument
  std::ifstream nodes_file(argv[1]);
  // Interpret each line of the nodes_file as a 3D Point and add to the Graph
  Point p;
  while (CS207::getline_parsed(nodes_file, p))
    nodes.push_back(graph.add_node(p));

  // Create a tets_file from the second input argument
  std::ifstream tets_file(argv[2]);
  // Interpret each line of the tets_file as four ints which refer to nodes
  std::array<int,4> t;
  while (CS207::getline_parsed(tets_file, t))
    for (unsigned i = 1; i < t.size(); ++i)
      for (unsigned j = 0; j < i; ++j)
        graph.add_edge(nodes[t[i]], nodes[t[j]]);

  // Print out the stats
  std::cout << graph.num_nodes() << " " << graph.num_edges() << std::endl;

  // Launch the SDLViewer
  CS207::SDLViewer viewer;
  auto node_map = viewer.empty_node_map(graph);
  SparsePredicate sparse_pred;

  // Use make_filtered helper to make iterators for passing to add_nodes and
  // launch the viewer.
  auto input_it_begin = make_filtered(graph.node_begin(),
      graph.node_end(), sparse_pred);
  auto input_it_end = make_filtered(graph.node_end(), graph.node_end(),
      sparse_pred);
  viewer.add_nodes(input_it_begin, input_it_end, node_map);
  viewer.add_edges(graph.edge_begin(), graph.edge_end(), node_map);
  viewer.launch();

  return 0;
}
