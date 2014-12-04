#include <fstream>
#include <iterator>

#include "CS207/SDLViewer.hpp"
#include "CS207/Util.hpp"

#include "Graph.hpp"

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
  //graph.add_node(Point(1,1,1), 3);
  //printf("val is %d\n", graph.node(0).value());
  //graph.add_node(Point(1,2,1));
  //printf("default val is: %d\n", graph.node(1).value());
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
  while (CS207::getline_parsed(tets_file, t)) {
    for (unsigned i = 1; i < t.size(); ++i) {
      for (unsigned j = 0; j < i; ++j) {
        graph.add_edge(nodes[t[i]], nodes[t[j]]);
      }
    }
  }

  // Print out the stats
  std::cout << "\n\n" << graph.num_nodes() << " " << graph.num_edges() << std::endl;
  // Launch the SDLViewer
  CS207::SDLViewer viewer;
  auto node_map = viewer.empty_node_map(graph);
  viewer.add_nodes(graph.node_begin(), graph.node_end(), node_map);
  std::cout << "about to add edges\n";
  fflush(stdout);
  viewer.add_edges(graph.edge_begin(), graph.edge_end(), node_map);
  std::cout << "\n\n\n\nfinished adding edges\n\n\n";
  fflush(stdout);
  viewer.launch();

  /*for (auto ni = graph.node_begin(); ni != graph.node_end(); ++ni) {
    std::cout << "entering loop yo\n";
    fflush(stdout);
    auto node = *ni;
    std::cout << node.index() << "  ";
    std::cout << node.value() << std::endl;
  }*/
  //for (auto ii = graph.node(2).edge_begin(); ii != graph.node(2).edge_end(); ++ii) {
  //  std::cout << (*ii).node1().index() << " " << (*ii).node2().index() << "\n";
  //}

  // Set the viewer
  //viewer.draw_graph(graph);
  viewer.center_view();

  return 0;
}