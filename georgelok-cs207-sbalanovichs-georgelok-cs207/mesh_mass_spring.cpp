/**
 * @file mesh_mass_spring.cpp
 * Implementation of mass-spring system using Mesh
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

#include "Mesh.hpp"
#include "Point.hpp"
#include "../jb_parallel.hpp"
#include <omp.h>

// Gravity in meters/sec^2
static constexpr double grav = 9.81;
static constexpr double K = 100.;
static constexpr double C = 30.;

/** Custom structure of data to store with Nodes */
struct NodeData {
  Point velocity;  //< Node velocity
  double mass;     //< Node mass
};

typedef Mesh<NodeData, double, int> MeshType;
typedef typename MeshType::node_type Node;
typedef typename MeshType::edge_type Edge;
typedef typename MeshType::triangle_type Triangle;

template <typename F>
struct vmod {
  F force;
  double dt;
  double t;
  void operator () (Node n) {
   // n.value().velocity += force(n, t) * (dt / n.value().mass);
  }
  vmod(F forcei, double dti, double ti) :
      force(forcei), dt(dti), t(ti) {};
};

template <typename F>
struct pmod {
  double dt;
  void operator () (Node n) {
    n.position() += n.value().value_.velocity * dt;
  }
  pmod(double dti) : dt(dti) {};
};

// template<typename Tri>
// struct msv_mod {
//   volume
//   void operator () (Tri t) {
//     double volume = 0.;
//     Point sn = t.surface_normal();
//     int val = t.value();
//     if ((sn.z < 0 && val > 0) || (sn.z > 0 && val < 0))
//       sn *= -1;
//     // Increment the volume
//     volume += (sn.z * val * t.area()) * 
//               ((t.node(0).position().z + 
//                 t.node(1).position().z + 
//                 t.node(2).position().z) / 3.);
//   }
// };



/** Change a mesh's nodes according to a step of the symplectic Euler
 *    method with the given node force.
 * @param[in,out] m      Mesh
 * @param[in]     t      The current time (useful for time-dependent forces)
 * @param[in]     dt     The time step
 * @param[in]     force  Function object defining the force per node
 * @return the next time step (usually @a t + @a dt)
 *
 * @tparam G::node_value_type supports any struct with velocity and mass
 * @tparam F is a function object called as @a force(n, @a t),
 *           where n is a node of the mesh and @a t is the current time.
 *           @a force must return a Point representing the force vector on Node
 *           at time @a t.
 */
template <typename M, typename F, typename C>
double symp_euler_step(M& m, double t, double dt, F force, C constraints) {
  // Compute the {n+1} node positions
  /*for (auto it = m.node_begin(); it != m.node_end(); ++it) {
    auto n = *it;

    // Update the position of the node according to its velocity
    // x^{n+1} = x^{n} + v^{n} * dt
    n.position() += n.value().value_.velocity * dt;
  }*/
  pmod<F> pm(dt);
  std::for_each(m.node_begin(), m.node_end(), pm);

  // Enfore constraints after position update but before force calculation.
  constraints(m, t);

  // Compute the {n+1} node velocities
  for (auto it = m.node_begin(); it != m.node_end(); ++it) {
    auto n = *it;
    // v^{n+1} = v^{n} + F(x^{n+1},t) * dt / m
    n.value().value_.velocity += force(n, t) * (dt / n.value().value_.mass);
  }

  return t + dt;
}

/** Set the direction flags on each triangle in order to 
 * keep track of the direction in which they are facing
 * @param[in] mesh, the current Mesh
 *
 * @pre: @a mesh is a valid mesh 
 * @post: for i in range(mesh), i.value() = -1 if i faces away from center and 1 otherwise
 *
 * Complexity: O(n)
 */
template <typename M>
void set_tri_directions(M& mesh) {
  // Set flags for surface normals
  // Approximate center of sphere
  Point sphere_center = Point(0,0,0);
  for (auto it = mesh.node_begin(); it != mesh.node_end(); ++it) {
    sphere_center += (*it).position();
  }
  sphere_center /= mesh.num_nodes();

  for(auto it = mesh.triangle_begin(); it != mesh.triangle_end(); ++it) {
    // Approximate center of triangle
    auto tri = (*it);
    Point tri_center = Point(0,0,0);
    for(unsigned i = 0; i < 3; ++i) {
      tri_center  += tri.node(i).position();
    }  
    tri_center /= 3.0;

    // Determine direction of surface normal.
    Point d_vector = tri_center - sphere_center;
    if(d_vector.z * tri.surface_normal().z > 0) {
      tri.value() = 1;
    } else {
      tri.value() = -1; 
    }
  }
}

/** A boolean functor that compares two nodes
 *  based on their velocity values and returns true
 *  if the first node's velocity is smaller than the
 *  second one's */
struct VelComparator {
  template <typename Node>
  bool operator()(const Node& n1, const Node& n2) const {
    return (norm(n1.value().value_.velocity) < norm(n2.value().value_.velocity));
  }
};

/** A NodeColor functor that colors nodes based off
 *  of node velocities */
struct VelHeatMap {
  public:
    template <typename Node>
    // Returns a heat value for a node based off velocity(n)
    CS207::Color operator()(const Node& n) const {
      return CS207::Color::make_heat(norm(n.value().value_.velocity)/max_vel_);
    }
    /* Constructor */
    VelHeatMap(double max_vel) : max_vel_(max_vel) {}
  private:
    const double max_vel_;
};

struct volume_inc {
  double operator () (Triangle t) {
    double volume = 0.;
    Point sn = t.surface_normal();
    int val = t.value();
    if ((sn.z < 0 && val > 0) || (sn.z > 0 && val < 0))
      sn *= -1;
    // Increment the volume
    volume += (sn.z * val * t.area()) * 
              ((t.node(0).position().z + 
                t.node(1).position().z + 
                t.node(2).position().z) / 3.);
    return volume;
  }
};

/** Calculate the volume of the mesh shape
 * @param[in] m, the Mesh
 * 
 * @return the volume of the mesh
 *
 * Complexity: O(n)
 */
template <typename M>
double mesh_shape_volume(M* m) {
  double volume = 0.;
  volume_inc vi;
  jb_parallel::parallel_reduction(m->triangle_begin(), m->triangle_end(), vi, volume);
  return volume;
}

/** Air Pressure force functor for the mesh.
 * Return the force being applied to @a n at time @a t.
 */
struct AirPressureForce {
  // Constructor. Establishes the current sphere's volume
  AirPressureForce(double volume) : volume_(volume) {}

  // Return the air pressure force being applied to @a n at time @a t.
  Point operator()(Node n, double t) {
    (void) t;

    Point n_i = Point(0,0,0);

    for(unsigned i = 0; i < n.value().neighbor_triangles_.size(); ++i) {
      auto tri = n.value().neighbor_triangle(i);
      n_i += (C) * tri.area() *  tri.value() * tri.surface_normal();
    }
    return n_i;
  }
  private:
    double volume_;
};

/** Gravity force functor for the mesh.
 * Return the force being applied to @a n at time @a t.
 */
struct GravityForce {
  Point operator()(Node n, double t) {
    (void) t;
    return Point(0,0, -grav * n.value().value_.mass);
  }
};

/** Mass spring force functor for the mesh.
 * Return the force being applied to @a n at time @a t.
 */
struct MassSpringForce {
  Point operator()(Node n, double t) {
    (void) t;

    // initialize force
    Point force = Point(0,0,0);

    // calculate spring force
    for(auto it = n.edge_begin(); it != n.edge_end(); ++it) {
      MeshType::edge_type temp_edge = *it;
      MeshType::node_type temp_node = temp_edge.node2();
      Point diff = n.position() - temp_node.position();
      force += (-K) * (diff) / norm(diff) * (norm(diff) - temp_edge.value().value_);      
    }
    return force;
  }
};

/** Wind force functor for the mesh.
 * Return the force being applied to @a n at time @a t.
 */
struct WindForce {
  // Constructor.  Establishes constants for interaction and air velocty
  WindForce(const double& c, const double& w) : c_(c), w_(w) {}; 

  // Return the wind force being applied to @a n at time @a t.
  Point operator()(Node n, double t) {
    (void) t;

    // initialize, calculate, and normalize the n_i surface normal
    Point n_i = Point(0,0,0);   
    for(unsigned i = 0; i < n.value().neighbor_triangles_.size(); ++i) {
      n_i += n.value().neighbor_triangle(i).surface_normal();
    }
    // n_i = n_i /norm(n_i);
    // return the force based on the wind formula
    return c_ * ((w_ - n.value().value_.velocity) * n_i) * n_i;
  }
private:
    double c_;
    double w_;
};

/** Damping force functor for the mesh.
 * Return the force being applied to @a n at time @a t.
 */
struct DampingForce {
  double c_;

  // Constructor. Establishes constant for damping.
  DampingForce(const double& c) : c_(c) {};  

  // Return the damping force being applied to @a n at time @a t.
  Point operator()(Node n, double t) {
    (void) t;
    return -c_ * n.value().value_.velocity;
  }
};



/** MetaForce construction as shown in class
 * Constructs a new force given two forces.
 * @param[in] force_1      First force
 * @param[in] force_2      Second force
 * 
 * @return a function object such output() = force_1_() + force_2_()
 *
 * @tparam force_1, force_2 are function objects called as @a force_#(n, @a t),
 *           where n is a node of the mesh and @a t is the current time.
 *           @a force must return a Point representing the force vector on Node
 *           at time @a t.
 */
template<typename F1, typename F2>
struct MetaForce {
  Point operator()(Node n, double t) {
    return force_1_(n, t) + force_2_(n, t);
  }

  MetaForce(F1 force_1, F2 force_2) 
    : force_1_(force_1), force_2_(force_2) {
  }

  private:
    F1 force_1_;
    F2 force_2_;
};

/** Combines two forces to create a new force
 *
 * @pre  Inputs can be called as @a fun(n, @a t),
 *        where n is a node and t is time.
 * Returns a MetaForce
 */
template<typename F1, typename F2>
MetaForce<F1, F2> make_combined_force(F1 force_1, F2 force_2) {
  return MetaForce<F1, F2>(force_1, force_2);
}

/** Combines three forces to create a new force
 *
 * @pre  Inputs can be called as @a fun(n, @a t),
 *        where n is a node and t is time.
 * Returns a MetaForce
 */
template<typename F1, typename F2, typename F3>
MetaForce<F1, MetaForce<F2, F3>> make_combined_force(F1 force_1, F2 force_2, F3 force_3) {
  MetaForce<F2, F3> f2_f3_combo = MetaForce<F2, F3>(force_2, force_3);
  return MetaForce<F1, MetaForce<F2, F3>>(force_1, f2_f3_combo);
}

/** Checks nodes in the graph to see if they are at a certain point and fixes them
 * @post n.position() == n.position() is true for all nodes at (0,0,0) and (0,0,1)
 * O(n) time
 * */
template <typename G>
struct FixedConstraint {
  void operator()(G& g, double t) {
    (void) t;
    for(auto it = g.node_begin(); it != g.node_end(); ++it) {
      if((*it).position() == Point(0,0,0) || (*it).position() == Point(1,0,0) ) {
        (*it).value().value_.velocity = Point(0,0,0);
      }
    }
  }
};

/** Constrains points to be above z = -0.75
 * O(N) Time
 */
template <typename M>
struct PlaneConstraint {

  void operator()(M& m, double t) {
    (void) t;
    for(auto it = m.node_begin(); it != m.node_end(); ++it) {
      if (dot((*it).position(), Point(0,0,1)) < -2) {
        (*it).position().z = -2;
        (*it).value().value_.velocity.z = 0;
      }
    }
  }
};

/** Constrains points to be outside sphere of radius r with center c
 * O(N) Time
 */
template <typename M>
struct SphereConstraint {
  SphereConstraint(Point c, double r) : c_(c), r_(r) {}

  void operator()(M& m, double t) {
    (void) t;
    for(auto it = m.node_begin(); it != m.node_end(); ++it) {
      if ( norm((*it).position() - c_) < r_ ) {
        Point Ri = ((*it).position() - c_) / norm((*it).position() - c_);
        (*it).position() = Ri * 0.15 + c_;
        (*it).value().value_.velocity -= dot((*it).value().value_.velocity, Ri)*Ri;
      }
    }
  }
  private:
    Point c_;
    double r_;
};

/** Constrains points to radius of r 
 * O(N) Time
 */
template <typename M>
struct BallConstraint {

  // Constructor. Establishes the ball radius constant
  BallConstraint(double r) : r_(r) {}

  void operator()(M& m, double t) {
    (void) t;
    Point sphere_center = Point(0,0,0);
    // Calculate center of sphere
    for(auto it = m.node_begin(); it != m.node_end(); ++it) {
      sphere_center += (*it).position();
    }
    sphere_center /= m.num_nodes();

    // Constrain nodes to a particular radius.
    for(auto it = m.node_begin(); it != m.node_end(); ++it) {
      Point vector = (*it).position() - sphere_center;
      if(norm(vector) > r_) {
        vector = vector / norm(vector) * r_;
        (*it).position() = vector + sphere_center;
      }

    }
  }
  private:
    double r_;
};

/** MetaConstraint in the same vein as the MetaForce
 * Constructs a new constraint given two constraints.
 * @param[in] constraint_1      First constraint
 * @param[in] constraint_2      Second constraint
 * 
 * @return a function object such output() = constraint_1_() + constraint_2_()
 *
 * @tparam constraint_1, constraint_2 are function objects called as @a pconstraint_#(g, @a t),
 *           where g is a mesh and @a t is the current time.
 *            no return value
 */
template<typename M, typename C1, typename C2>
struct MetaConstraint {
  void operator()(M& m, double t) {
    constraint_1_(m, t);
    constraint_2_(m, t);
  }

  MetaConstraint(C1 constraint_1, C2 constraint_2) 
    : constraint_1_(constraint_1), constraint_2_(constraint_2) {
  }

  C1 constraint_1_;
  C2 constraint_2_;
};

/** Combines two constraints to create a new constraint
 *
 * @pre  Inputs can be called as @a fun(g, @a t),
 *        where g is a node and t is time.
 * Returns a MetaConstraint
 */
template<typename C1, typename C2, typename M = MeshType>
MetaConstraint<M, C1, C2> make_combined_constraint(C1 constraint_1, C2 constraint_2) {
  return MetaConstraint<M, C1, C2>(constraint_1, constraint_2);
}

/** Combines three constraints to create a new constraint
 *
 * @pre  Inputs can be called as @a fun(g, @a t),
 *        where g is a node and t is time.
 * Returns a MetaConstraint
 */
template<typename C1, typename C2, typename C3, typename M = MeshType>
MetaConstraint<M, C1, MetaConstraint<M, C2, C3>> make_combined_constraint(C1 constraint_1, C2 constraint_2, C3 constraint_3) {
  MetaConstraint<M, C2, C3> c2_c3_combo = MetaConstraint<M, C2, C3>(constraint_2, constraint_3);
  return MetaConstraint<M, C1, MetaConstraint<M, C2, C3>>(constraint_1, c2_c3_combo);
}

int main(int argc, char** argv) {
  // Check arguments
  if (argc < 3) {
    std::cerr << "Usage: " << argv[0] << " NODES_FILE TETS_FILE\n";
    exit(1);
  }

  // Construct a mesh
  MeshType mesh;
  std::vector<typename MeshType::node_type> mesh_node;

  // Read all Points and add them to the Mesh
  std::ifstream nodes_file(argv[1]);
  Point p;
  while (CS207::getline_parsed(nodes_file, p)) {
    mesh_node.push_back(mesh.add_node(p));
  }

  // Read all mesh triangles and add them to the Mesh
  std::ifstream tris_file(argv[2]);
  std::array<int,3> t;
  while (CS207::getline_parsed(tris_file, t)) {
    mesh.add_triangle(mesh_node[t[0]], mesh_node[t[1]], mesh_node[t[2]]);
  }
  
  // Zero intial velocities and have constant density.
  for(MeshType::node_iterator it = mesh.node_begin(); it != mesh.node_end(); ++it) {
    MeshType::node_type temp_node = *it;
    temp_node.value().value_.velocity = Point(0,0,0);
    temp_node.value().value_.mass = 1 / double(mesh.num_nodes());
  }

  // Set initial lengths to initial distances.
  for(MeshType::edge_iterator it = mesh.edge_begin(); it != mesh.edge_end(); ++it) {
    MeshType::edge_type temp_edge = *it;
      temp_edge.value().value_ = temp_edge.length();
  }

  // Print out the stats
  //std::cout << mesh.num_nodes() << " " << mesh.num_edges() << std::endl;

  // Launch the SDLViewer
  //CS207::SDLViewer viewer;
  //auto node_map = viewer.empty_node_map(mesh);
  //viewer.launch();

  //viewer.add_nodes(mesh.node_begin(), mesh.node_end(), node_map);
  //viewer.add_edges(mesh.edge_begin(), mesh.edge_end(), node_map);

  //viewer.center_view();


  set_tri_directions(mesh);

  //std::cout << "Starting\n";

  // Begin the mass-spring simulation
  double dt = 0.0001;
  double t_start = 0.0;
  double t_end   = 0.001;
  double max_vel = -UINT_MAX;
  {jb_parallel::Timer timer("");
  for (double t = t_start; t < t_end; t += dt) {
    auto force = make_combined_force(
        make_combined_force(MassSpringForce(), 
                            GravityForce()),
        AirPressureForce(mesh_shape_volume(&mesh)),
        // WindForce(3., 0.7),
        DampingForce(1/(mesh.num_nodes() * 10))   
    );
    auto constraints = make_combined_constraint(
        FixedConstraint<MeshType>(),
        BallConstraint<MeshType>(1.9),
        PlaneConstraint<MeshType>()
        );
    symp_euler_step(mesh, t, dt, force, constraints);
    // Clear the viewer's nodes and edges.
    //viewer.clear();
    //node_map.clear();

    // Update the viewer with new node positions and color
    // auto max_node = *std::max_element(mesh.node_begin(), 
    //                                               mesh.node_end(), 
    //                                               VelComparator());
    // max_vel = std::max(max_vel, norm(max_node.value().value_.velocity));

    // Update viewer with nodes' new positions and new edges.
    //viewer.add_nodes(mesh.node_begin(), mesh.node_end(), VelHeatMap(max_vel), node_map);
    //viewer.add_edges(mesh.edge_begin(), mesh.edge_end(), node_map);

    //viewer.set_label(t);
  }
  }

  return 0;
}
