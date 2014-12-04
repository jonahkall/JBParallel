/**
 * @file shortest_path.cpp
 * Test script for using our templated Graph to determine shortest paths.
 *
 * @brief Reads in two files specified on the command line.
 * First file: 3D Points (one per line) defined by three doubles
 * Second file: Tetrahedra (one per line) defined by 4 indices into the point
 * list
 */

#include <vector>
#include <fstream>
#include <queue>

#include "CS207/SDLViewer.hpp"
#include "CS207/Util.hpp"
#include "CS207/Color.hpp"

#include "Graph.hpp"

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

  // Constructor
  filter_iterator(const Pred& p, const It& first, const It& last)
      : p_(p), it_(first), end_(last) {
    // HW1 #4: YOUR CODE HERE
    if (first == last)
      return;
    while (!p_(*it_)) {
      ++it_;
      if (first == last)
        return;
    }
  }

  // HW1 #4: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:

  // ADD SPEC
  value_type operator*() const {
    return *it_;
  }

  // ADD SPEC
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

  // ADD SPEC -- notes that we assume preds are the same, else undefined behavior
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

// This predicate makes the graph more sparse by getting rid of about 2/3 of the nodes.
struct SparsePredicate {
  template <typename NODE>
  bool operator()(const NODE& n) {
    return (n.index() % 3 == 1);
  }
};

/** Comparator that compares the distance from a given point p.
 */
struct MyComparator {
   Point p_;
   MyComparator(const Point& p) : p_(p) {
   };

   template <typename NODE>
   bool operator()(const NODE& node1, const NODE& node2) const {
    // Compute the Euclidean distance between node1 and p_ and node2 and p_.
    auto d_x_1 = (node1.position().x - p_.x);
    auto d_x_2 = (node2.position().x - p_.x);
    auto d_y_1 = (node1.position().y - p_.y);
    auto d_y_2 = (node2.position().y - p_.y);
    auto d_z_1 = (node1.position().z - p_.z);
    auto d_z_2 = (node2.position().z - p_.z);
    auto d_1 = d_x_1*d_x_1 + d_y_1*d_y_1 + d_z_1*d_z_1;
    auto d_2 = d_x_2*d_x_2 + d_y_2*d_y_2 + d_z_2*d_z_2;
    if (d_1 < d_2)
      return true;
    return false;
  }
};

struct ColorFunctor {
  Point::value_type max_;
  ColorFunctor(Point::value_type max) : max_(max) {};
  template <typename NODE>
    CS207::Color operator()(const NODE& n) {
      return CS207::Color::make_heat(n.value()/max_);
    }
};

/** Calculate shortest path lengths in @a g from the nearest node to @a point.
 * @param[in,out] g Input graph
 * @param[in] point Point to find the nearest node to.
 * @post Graph has modified node values indicating the minimum path length
 *           to the nearest node to @a point
 * @post Graph nodes that are unreachable to the nearest node to @a point have
 *           the value() -1.
 * @return The maximum path length found.
 *
 * Finds the nearest node to @a point and treats that as the root node for a
 * breadth first search.
 * This sets node's value() to the length of the shortest path to
 * the root node. The root's value() is 0. Nodes unreachable from
 * the root have value() -1.
 */
int shortest_path_lengths(Graph<int, int>& g, const Point& point) {
  MyComparator comp(point);
  typedef Graph<int, int>::Node Nd;
  Nd n = *std::min_element(g.node_begin(), g.node_end(), comp);
  std::vector<bool> visited(g.num_nodes());
  std::fill(visited.begin(), visited.end(), false);
  n.value() = 0;
  typedef Graph<int, int>::size_type St;

  // Execute a breadth first search using a queue of node ids.
  std::queue<St> q;
  q.push(n.index());
  auto max = 0;
  while (q.size() != 0) {
    St tmp = q.front();
    q.pop();
    visited[tmp] = true;
    for (auto ii = g.node(tmp).edge_begin(); ii != g.node(tmp).edge_end(); ++ii) {
        if (visited[(*ii).node2().index()] == false) {
          q.push((*ii).node2().index());
          visited[(*ii).node2().index()] = true;
          (*ii).node2().value() = g.node(tmp).value() + 1;
          if ((*ii).node2().value() > max)
            max = (*ii).node2().value();
        }
    }
  }
  for (size_t i = 0; i < visited.size(); ++i) {
    if (visited[i] == false) {
      g.node(i).value() = -1;
    }
  }
  return max;
}



int main(int argc, char** argv)
{
  // Check arguments
  if (argc < 3) {
    std::cerr << "Usage: " << argv[0] << " NODES_FILE TETS_FILE\n";
    exit(1);
  }

  // Construct a Graph
  typedef Graph<int,int> GraphType;
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

  auto x = shortest_path_lengths(graph, Point(-1,0,1));
  ColorFunctor cf(x);

  // Launch the SDLViewer
  CS207::SDLViewer viewer;
  auto node_map = viewer.empty_node_map(graph);
  SparsePredicate sparse_pred;
  filter_iterator<SparsePredicate, GraphType::NodeIterator> input_it_begin(
      sparse_pred, graph.node_begin(), graph.node_end());
  filter_iterator<SparsePredicate, GraphType::NodeIterator> input_it_end(
      sparse_pred, graph.node_end(), graph.node_end());
  viewer.add_nodes(input_it_begin, input_it_end, cf, node_map);
  viewer.add_edges(graph.edge_begin(), graph.edge_end(), node_map);
  viewer.launch();
  return 0;
}

