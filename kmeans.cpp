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

// Node value type for graph representation of clustering
// problem.  This allows us to quickly access which cluster a
// node 
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

// This functor takes in a node and the vector of cluster
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


/** BRIEF DESCRIPTION
 * @param NAME DESCRIPTION
 * @return DESCRIPTION OF RETURN VALUE
 * @pre PRECONDITION
 * @post POSTCONDITION
 *
 * LONGER DESCRIPTION, RUNTIME, EXAMPLES, ETC 
 */
template<typename GraphType>
void kmeans_f(GraphType& g, int k, bool opt) {
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


  // While some cluster assignment has changed, continue iterating.
  while (1) {
    // Can't easily be parallelized because of data dependencies.
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
    if (!opt) {
      std::for_each(g.node_begin(), g.node_end(), ac);
    }
    else {
      jb_parallel::for_each(g.node_begin(), g.node_end(), ac);
    }
    // This loop can also be made parallel by using parallel reducer,
    // but this did not lend speedup.
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

  // Construct a Graph
  typedef Graph<node_value_type, int> GraphType;
  GraphType graph;
  std::vector<GraphType::node_type> nodes;

  // Create a nodes_file from the first input argument
  // Specify a 1 at the command line for parallelization on,
  // and a 0 for parallelization off.
  std::string PATH_TO_NODES;
  bool opt = true;
  if (argc == 2) {
    PATH_TO_NODES = "large_clustering_problem2.nodes";
    if (std::stoi(argv[1]) == 0) {
      opt = false;
    }
  }
  else if (argc == 3){
    PATH_TO_NODES = argv[1];
    if (std::stoi(argv[2]) == 0) {
      opt = false;
    }
  }
  else {
    std::cout << "Usage: ./kmeans NODES_FILE 0/1" << std::endl;
    exit(1);
  }

  std::ifstream nodes_file(PATH_TO_NODES);

  // Interpret each line of the nodes_file as a 3D Point and add to the Graph
  Point p;
  while (CS207::getline_parsed(nodes_file, p))
    nodes.push_back(graph.add_node(p));

  // Print out the stats
  std::cout << graph.num_nodes() << " " << graph.num_edges() << std::endl;

  // Define the number of clusters, and begin clustering.
  int num_clusters = 4;
  { jb_parallel::Timer timer("KMeans Clustering");
    kmeans_f(graph, num_clusters, opt);
  }

  // Initialize functor to view the clusters.
  ClusterColorFunctor cf(num_clusters);

  // Launch the SDLViewer to visualize clustering
  // Note that for simple, dense collections of nodes like most of those
  // we deal with in this class, this will simply produce a voronoi
  // tessellation.
  CS207::SDLViewer viewer;
  auto node_map = viewer.empty_node_map(graph);
  viewer.add_nodes(graph.node_begin(), graph.node_end(), cf, node_map);
  viewer.launch();
  return 0;
}
