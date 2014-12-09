#include <vector>
#include <fstream>
#include <queue>
#include <iostream>
#include <cstdlib>

#include "CS207/SDLViewer.hpp"
#include "CS207/Util.hpp"
#include "CS207/Color.hpp"
#include "Point.hpp"

#include "Graph.hpp"
#include "jb_parallel.hpp"

struct node_value_type {
  int cluster_assignment;
  bool changed;
};

typedef Graph<node_value_type, int> GraphType;

template<typename NodeType>
struct find_centers {
  std::vector<Point>* centers_;
  std::vector<int>* center_counts_;
  void operator () (NodeType n) {
    int ca = n.value().cluster_assignment;
    ++((*center_counts_)[ca]);
    (*centers_)[ca] += n.position();
  }
  find_centers(std::vector<Point>* centers, std::vector<int>* center_counts)
      : centers_(centers), center_counts_(center_counts){};
};

template<typename NodeType>
struct assign_centers {
  std::vector<Point> centers_;
  int k_;
  void operator () (NodeType n) {
    int min_cluster_index = -1;
    float min_error_seen = 1000000;
    for (int i = 0; i < k_; ++i) {
      auto curr_error = norm(n.position() - centers_[i]);
      if (curr_error < min_error_seen) {
        min_cluster_index = i;
        min_error_seen = curr_error;
      }
    }
    if (n.value().cluster_assignment != min_cluster_index) {
      n.value().cluster_assignment = min_cluster_index;
      n.value().changed = true;
      return;
    }
    else {
      n.value().changed = false;
    }
  }
  assign_centers(std::vector<Point> centers, int k) : centers_(centers), k_(k){};
};

template<typename GraphType>
void kmeans_f(GraphType& g, int k) {
  srand((unsigned)time(0));
  for (auto it = g.node_begin(); it != g.node_end(); ++it) {
    (*it).value().cluster_assignment = rand() % k;
  }
  std::vector<Point> centers(k);
  std::vector<int> center_counts(k);
  for (int i = 0; i < k; ++i) {
    centers[i] = Point(0,0,0);
    center_counts[i] = 0;
  }


  // While none of the cluster assignments change
  while (1) {
    //find_centers<typename GraphType::Node> fc(&centers, &center_counts);
    // Can't be parallelized because of data dependencies.
    for (auto it = g.node_begin(); it != g.node_end(); ++it) {
      int ca = (*it).value().cluster_assignment;
      ++center_counts[ca];
      centers[ca] += (*it).position();
      (*it).value().changed = true;
    }

    for (int i = 0; i < k; ++i) {
      centers[i] = centers[i]/static_cast<float>(center_counts[i]);
    }

    // loop over all data and assign things to their closest center
    bool some_assignment_changed = false;
    assign_centers<typename GraphType::Node> ac(centers, k);
    jb_parallel::for_each(g.node_begin(), g.node_end(), ac);
    for (auto it = g.node_begin(); it != g.node_end(); ++it) {
      if ((*it).value().changed) {
        some_assignment_changed = true;
        break;
      }
    }
    if (!some_assignment_changed) {
      break;
    }
  }
}

// Simply converts the distance in the graph to a color value via
// make_heat.  Path length values are normalized by the maximum path length.
struct ClusterColorFunctor {
  float k_;
  ClusterColorFunctor(float k) : k_(k) {};
  template <typename NODE>
    CS207::Color operator()(const NODE& n) {
      //std::cout << n.value().cluster_assignment << " out of " << k_ << std::endl;
      return CS207::Color::make_heat(
          static_cast<float>(n.value().cluster_assignment)/k_);
    }
};

int main(int argc, char** argv)
{
  // Check arguments
  if (argc < 2) {
    std::cerr << "Usage: " << argv[0] << " NODES_FILE\n";
    exit(1);
  }

  // Construct a Graph

  typedef Graph<node_value_type, int> GraphType;
  GraphType graph;
  std::vector<GraphType::node_type> nodes;

  // Create a nodes_file from the first input argument
  std::ifstream nodes_file(argv[1]);
  // Interpret each line of the nodes_file as a 3D Point and add to the Graph
  Point p;
  while (CS207::getline_parsed(nodes_file, p))
    nodes.push_back(graph.add_node(p));

  // Print out the stats
  std::cout << graph.num_nodes() << " " << graph.num_edges() << std::endl;

  // Assign shortest path lengths and declare a color functor with the max
  // path length as its normalization constant.

  int num_clusters = 4;
  { jb_parallel::Timer timer("KMeans Clustering");
    kmeans_f(graph, num_clusters);
  }
  ClusterColorFunctor cf(num_clusters);

  // Launch the SDLViewer
  CS207::SDLViewer viewer;
  auto node_map = viewer.empty_node_map(graph);
  viewer.add_nodes(graph.node_begin(), graph.node_end(), cf, node_map);
  viewer.launch();
  return 0;
}
