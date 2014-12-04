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

// Simply converts the distance in the graph to a color value via
// make_heat.  Path length values are normalized by the maximum path length.
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

  // Assign shortest path lengths and declare a color functor with the max
  // path length as its normalization constant.
  auto x = shortest_path_lengths(graph, Point(-1,0,1));
  ColorFunctor cf(x);

  // Launch the SDLViewer
  CS207::SDLViewer viewer;
  auto node_map = viewer.empty_node_map(graph);
  viewer.add_nodes(graph.node_begin(), graph.node_end(), cf, node_map);
  viewer.add_edges(graph.edge_begin(), graph.edge_end(), node_map);
  viewer.launch();
  return 0;
}
