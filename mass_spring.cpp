/**
 * @file mass_spring.cpp
 * Implementation of mass-spring system using Graph
 *
 * @brief Reads in two files specified on the command line.
 * First file: 3D Points (one per line) defined by three doubles
 * Second file: Tetrahedra (one per line) defined by 4 indices into the point
 * list
 */

#include <fstream>

#include "CS207/SDLViewer.hpp"
#include "CS207/Util.hpp"
#include "CS207/Color.hpp"

#include "Graph.hpp"
#include "Point.hpp"


// Gravity in meters/sec^2
static constexpr double grav = 9.81;

/** Custom structure of data to store with Nodes */
struct NodeData {
  Point velocity;  //< Node velocity
  double mass;     //< Node mass
};

/** Custom struct to represent data for each edge:
 */
struct EdgeData {
  double K; //< Edge spring constant
  double L; //< Edge rest length
};

typedef Graph<NodeData, EdgeData> GraphType;
typedef typename GraphType::node_type Node;
typedef typename GraphType::edge_type Edge;

/** Change a graph's nodes according to a step of the symplectic Euler
 *    method with the given node force.
 * @param[in,out] g      Graph
 * @param[in]     t      The current time (useful for time-dependent forces)
 * @param[in]     dt     The time step
 * @param[in]     force  Function object defining the force per node
 * @return the next time step (usually @a t + @a dt)
 *
 * @tparam G::node_value_type supports ???????? YOU CHOOSE
 * @tparam F is a function object called as @a force(n, @a t),
 *           where n is a node of the graph and @a t is the current time.
 *           @a force must return a Point representing the force vector on Node
 *           at time @a t.
 */
template <typename G, typename F, typename Constraint>
double symp_euler_step(G& g, double t, double dt, F force, Constraint c) {
  // Compute the {n+1} node positions
  for (auto it = g.node_begin(); it != g.node_end(); ++it) {
    auto n = *it;

    // Update the position of the node according to its velocity
    // x^{n+1} = x^{n} + v^{n} * dt
    n.position() += n.value().velocity * dt;
  }

  // Apply the constraint.
  c(g,t);

  // Compute the {n+1} node velocities
  for (auto it = g.node_begin(); it != g.node_end(); ++it) {
    auto n = *it;
    // v^{n+1} = v^{n} + F(x^{n+1},t) * dt / m
    n.value().velocity += force(n, t) * (dt / n.value().mass);
  }

  return t + dt;
}

// Function object which represents the downward force of gravity,
// -mg in the z direction.
struct GravityForce {
  Point operator()(Node n, double t) {
    (void) t;
    return Point(0,0, -n.value().mass * grav);
  }
};

// Function object which represents the spring force between nodes,
// using a sum of all the hooke's law forces on a node.
struct MassSpringForce {
  Point operator ()(Node n, double t) {
    (void) t;
    Point force_sum = Point(0,0,0);
    for (auto it = n.edge_begin(); it != n.edge_end(); ++it) {
      Point diff = n.position() - (*it).node2().position();
      force_sum += -(*it).value().K*(diff/norm(diff))*
          (norm(diff) - (*it).value().L);
    }
    return force_sum;
  }
};

// Damping force -cv_i which counters the current velocity of
// a given node.
struct DampingForce {
  Point operator () (Node n, double t) {
    (void) t;
    return -(1/n.parent_size()) * n.value().velocity;
  }
};

// Utility for testing, alwyas returns 0 force, not used anywhere currently.
struct ZeroForce {
  Point operator () (Node n, double t) {
    (void) n; (void) t;
    return Point(0,0,0);
  }
};

// Functor which represents a generalized combination of forces.
// Note: must satisfy the force concept.
template<typename force_type_1, typename force_type_2>
struct GeneralForce {
  Point operator () (Node n, double t) {
    (void) t;
    return f1_(n,t) + f2_(n,t);
  }
  force_type_1 f1_;
  force_type_2 f2_;
  GeneralForce(force_type_1 f1, force_type_2 f2) : f1_(f1), f2_(f2) {};
};

// Function which takes 2 function objects representing forces and returns
// a GeneralForce which represents their combination.
template<typename force_type_1, typename force_type_2>
GeneralForce<force_type_1, force_type_2> make_combined_force
    (force_type_1 f1, force_type_2 f2) {
  return GeneralForce<force_type_1, force_type_2>(f1, f2);
}

// Function which takes 3 function objects representing forces and returns a
// General Force representing their combination (note the template arguments:
// as talked about in class, GeneralForces can be nested arbitrarily, allowing 
// easy combination of lots of forces).
template<typename force_type_1, typename force_type_2, typename force_type_3>
GeneralForce<GeneralForce<force_type_1, force_type_2>,force_type_3>
    make_combined_force (force_type_1 f1, force_type_2 f2, force_type_3 f3) {
  return make_combined_force(make_combined_force(f1, f2), f3);
}

/** Force Functor which combines MassSpringForce, DampingForce, GravityForce 
 * for this problem set.
 */
struct Problem1Force {
  /** Return the force being applied to @a n at time @a t.
   *
   */
  Point operator () (Node n, double t) {
    if (n.position() == Point(0,0,0) || n.position() == Point(1,0,0)) {
      return Point(0,0,0);
    }
    return make_combined_force<GravityForce, MassSpringForce, DampingForce>
        (GravityForce(), MassSpringForce(),DampingForce())(n, t);
  }
};

// Functor which represents a combination of two constraints.  First applies
// constraint 1, then applies constraint 2.
template <class Constraint1, class Constraint2>
struct CombinedConstraint {
  void operator() (GraphType& g, double t) {
    (void) t;
    //linear search
    C1_(g,t);
    C2_(g,t);
  }
  Constraint1 C1_;
  Constraint2 C2_;
  CombinedConstraint(Constraint1 C1, Constraint2 C2) : C1_(C1), C2_(C2) {};
};

// Function object which represents a z-plane which the nodes cannnot
// pass through.
struct PlaneConstraint {
  void operator() (GraphType& g, double) {
    for (auto it = g.node_begin(); it != g.node_end(); ++it) {
      // If the constraint fails, fix the node.
      if (dot((*it).position(), Point(0,0,1)) < -0.75) {
        (*it).position().z = -0.75;
        (*it).value().velocity.z = 0;
      }
    }
  }
};

// Function object which represents an impassable sphere.
struct SphereConstraint {
  void operator() (GraphType& g, double) {
    Point c = Point(0.5, 0.5, -0.5);
    for (auto it = g.node_begin(); it != g.node_end(); ++it) {
      // If the constraint fails, fix the node.
      if (norm((*it).position() - c) < 0.15) {
        Point R_i = ((*it).position()-c)/norm((*it).position()-c);
        (*it).value().velocity =
            (*it).value().velocity - dot((*it).value().velocity, R_i)*R_i;
        (*it).position() =
            (((*it).position() - c)/norm((*it).position()-c)) * 0.15 + c;
      }
    }
  }
};

// Function object which represents a sphere as defined in the HW2 spec, and
// nodes which pass through will be deleted.
struct SphereDelConstraint {
  void operator() (GraphType& g, double) {
    auto end_it = g.node_end();
    for (auto it = g.node_begin(); it != end_it; ++it) {
      // If the constraint fails, delete the node.
      if (norm((*it).position() - Point(0.5,0.5,-0.5)) < 0.15) {
        int res = g.remove_node(*it);
        // Use decrement operator to leverage the fact that all indexes
        // after current index are decremented (the -1 check should be
        // unnecessary but makes sure that we actually deleted a node before
        // we decrement the iterator).
        if (res != -1) {
          --it;
        }
      }
      end_it = g.node_end();
    }
  }
};

int main(int argc, char** argv) {
  // Check arguments
  if (argc < 3) {
    std::cerr << "Usage: " << argv[0] << " NODES_FILE TETS_FILE\n";
    exit(1);
  }

  // Construct a graph
  GraphType graph;

  // Create a nodes_file from the first input argument
  std::ifstream nodes_file(argv[1]);
  // Interpret each line of the nodes_file as a 3D Point and add to the Graph
  std::vector<Node> nodes;
  Point p;
  while (CS207::getline_parsed(nodes_file, p))
    nodes.push_back(graph.add_node(p));

  // Create a tets_file from the second input argument
  std::ifstream tets_file(argv[2]);
  // Interpret each line of the tets_file as four ints which refer to nodes
  std::array<int,4> t;
  while (CS207::getline_parsed(tets_file, t)) {
    for (unsigned i = 1; i < t.size(); ++i) {
      graph.add_edge(nodes[t[0]], nodes[t[1]]);
      graph.add_edge(nodes[t[0]], nodes[t[2]]);

      // Diagonal edges: include as of HW2 #2
      graph.add_edge(nodes[t[0]], nodes[t[3]]);
      graph.add_edge(nodes[t[1]], nodes[t[2]]);

      graph.add_edge(nodes[t[1]], nodes[t[3]]);
      graph.add_edge(nodes[t[2]], nodes[t[3]]);
    }
  }

  // Set initial conditions for nodes,
  // construct Forces/Constraints
  auto mass = 1/static_cast<float>(graph.num_nodes());

  // Set initial velocities to zero and masses to 1/num_nodes().
  for (auto it = graph.node_begin(); it != graph.node_end(); ++it) {
    (*it).value().velocity = Point(0,0,0);
    (*it).value().mass = mass;
  }

  // Set rest lengths of edges equal to initial distance between the two nodes,
  // set spring constant to 100.
  for (auto ei = graph.edge_begin(); ei != graph.edge_end(); ++ei) {
    (*ei).value().L = norm((*ei).node1().position() - (*ei).node2().position());
    (*ei).value().K = 100.0;
  }

  // Print out the stats
  std::cout << graph.num_nodes() << " " << graph.num_edges() << std::endl;

  // Launch the SDLViewer
  CS207::SDLViewer viewer;
  auto node_map = viewer.empty_node_map(graph);
  viewer.launch();

  viewer.add_nodes(graph.node_begin(), graph.node_end(), node_map);
  viewer.add_edges(graph.edge_begin(), graph.edge_end(), node_map);
  viewer.center_view();

  // Begin the mass-spring simulation
  // Decrease timestep if you don't want grid3 to crash.
  double dt = 0.001;
  double t_start = 0.0;
  double t_end   = 5.0;

  // Combines the plane constraint and the sphere constraint with deletion.
  auto C = CombinedConstraint<PlaneConstraint, SphereDelConstraint>
      (PlaneConstraint(),SphereDelConstraint());

  for (double t = t_start; t < t_end; t += dt) {
    symp_euler_step(graph, t, dt, Problem1Force(), C);

    viewer.clear();
    node_map.clear();

    // Update viewer with nodes' new positions
    viewer.add_nodes(graph.node_begin(), graph.node_end(), node_map);
    viewer.add_edges(graph.edge_begin(), graph.edge_end(), node_map);

    viewer.set_label(t);

    // These lines slow down the animation for small graphs, like grid0_*.
    // Feel free to remove them or tweak the constants.
    if (graph.size() < 100)
      CS207::sleep(0.001);
  }

  return 0;
}
