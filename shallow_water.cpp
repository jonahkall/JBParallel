/**
 * @file shallow_water.cpp
 * Implementation of a shallow water system using Mesh
 *
 * @brief Reads in two files specified on the command line.
 * First file: 3D point list (one per line) defined by three doubles
 * Second file: Triangles (one per line) defined by 3 indices into the point list
 */

#include <fstream>
#include <cmath>

#include "CS207/SDLViewer.hpp"
#include "CS207/Util.hpp"
#include "CS207/Color.hpp"

#include "Point.hpp"
#include "Mesh.hpp"
#include "jb_parallel.hpp"

// Standard gravity (average gravity at Earth's surface) in meters/sec^2
static constexpr double grav = 9.80665;
typedef unsigned uid_type;

/** Water column characteristics */
struct QVar {
  double h;   // Height of column
  double hx;  // Height times average x velocity of column
  double hy;  // Height times average y velocity of column

  /** Default constructor.
   *
   * A default water column is 1 unit high with no velocity. */
  QVar()
      : h(1), hx(0), hy(0) {
  }
  /** Construct the given water column. */
  QVar(double h_, double hx_, double hy_)
      : h(h_), hx(hx_), hy(hy_) {
  }

  // Basic arithmetic operation operators on QVars.
  QVar& operator*=(double b) {
    h *= b;
    hx *= b;
    hy *= b;
    return *this;
  }

  QVar& operator/=(double b) {
    h /= b;
    hx /= b;
    hy /= b;
    return *this;
  }

  QVar& operator-=(double b) {
    h -= b;
    hx -= b;
    hy -= b;
    return *this;
  }

  QVar& operator+=(double b) {
    h += b;
    hx += b;
    hy += b;
    return *this;
  }

  QVar& operator+=(const QVar& b) {
    h += b.h;
    hx += b.hx;
    hy += b.hy;
    return *this;
  }

  QVar& operator-=(const QVar& b) {
    h -= b.h;
    hx -= b.hx;
    hy -= b.hy;
    return *this;
  }

};

QVar operator*(QVar a, double b) {
  return a *= b;
}

QVar operator+(QVar a, QVar b) {
  return a += b;
}

QVar operator*(double b, QVar a) {
  return a *= b;
}

QVar operator/(QVar a, double b) {
  return a/=b;
}

// HW4B: Placeholder for Mesh Type!
// Define NodeData, EdgeData, TriData, etc
// or redefine for your particular Mesh

struct node_value {
  QVar q;
};

struct tri_value {
	QVar q_k_bar_;
};

struct edge_value {
  QVar flux_;
  Point norm_vec_;
};


typedef Mesh<node_value,edge_value,tri_value> MeshType;


/** Function object for calculating shallow-water flux.
 *          |n
 *   T_k    |---> n = (nx,ny)   T_m
 *   QBar_k |                   QBar_m
 *          |
 * @param[in] nx, ny Defines the 2D outward normal vector n = (@a nx, @a ny)
 *            from triangle T_k to triangle T_m. The length of n is equal to the
 *            the length of the edge, |n| = |e|.
 * @param[in] dt The time step taken by the simulation. Used to compute the
 *               Lax-Wendroff dissipation term.
 * @param[in] qk The values of the conserved variables on the left of the edge.
 * @param[in] qm The values of the conserved variables on the right of the edge.
 * @return The flux of the conserved values across the edge e
 */
struct EdgeFluxCalculator {
  QVar operator()(double nx, double ny, double dt,
                  const QVar& qk, const QVar& qm) {
    // Normalize the (nx,ny) vector
    double n_length = std::sqrt(nx*nx + ny*ny);
    nx /= n_length;
    ny /= n_length;

    // The velocities normal to the edge
    double wm = (qm.hx*nx + qm.hy*ny) / qm.h;
    double wk = (qk.hx*nx + qk.hy*ny) / qk.h;

    // Lax-Wendroff local dissipation coefficient
    double vm = sqrt(grav*qm.h) + sqrt(qm.hx*qm.hx + qm.hy*qm.hy) / qm.h;
    double vk = sqrt(grav*qk.h) + sqrt(qk.hx*qk.hx + qk.hy*qk.hy) / qk.h;
    double a  = dt * std::max(vm*vm, vk*vk);

    // Helper values
    double scale = 0.5 * n_length;
    double gh2   = 0.5 * grav * (qm.h*qm.h + qk.h*qk.h);

    // Simple flux with dissipation for stability
    return QVar(scale * (wm*qm.h  + wk*qk.h)           - a * (qm.h  - qk.h),
                scale * (wm*qm.hx + wk*qk.hx + gh2*nx) - a * (qm.hx - qk.hx),
                scale * (wm*qm.hy + wk*qk.hy + gh2*ny) - a * (qm.hy - qk.hy));
  }
};

/** Node position function object for use in the SDLViewer. */
struct NodePosition {
  template <typename NODE>
  Point operator()(const NODE& n) {
    return Point(n.position().x, n.position().y, n.value().q.h);
  }
};

template<typename Tri>
struct FluxUpdater{
  double dt_;
  void operator () (Tri t) {
    QVar flux_sum = QVar(0,0,0);

    // set the fluxes correctly
    if (t == t.edge1().triangle1()) {
      flux_sum += t.edge1().value().flux_;
    }
    else {
      flux_sum -= t.edge1().value().flux_;
    }

    if (t == t.edge2().triangle1()) {
      flux_sum += t.edge2().value().flux_;
    }
    else {
      flux_sum -= t.edge2().value().flux_;
    }

    if (t == t.edge3().triangle1()) {
      flux_sum += t.edge3().value().flux_;
    }
    else {
      flux_sum -= t.edge3().value().flux_;
    }

    t.value().q_k_bar_ -= (dt_/t.area()) * flux_sum;
  }
  FluxUpdater(double dt) : dt_(dt) {};
};
/** Integrate a hyperbolic conservation law defined over the mesh m
 * with flux functor f by dt in time.
 */
template <typename MESH, typename FLUX>
double hyperbolic_step(MESH& m, FLUX& f, double t, double dt) {
  // HW4B: YOUR CODE HERE
  // Step the finite volume model in time by dt.

  // Pseudocode:
  // Compute all fluxes. (before updating any triangle Q_bars)
  // For each triangle, update Q_bar using the fluxes as in Equation 8.
  //  NOTE: Much like symp_euler_step, this may require TWO for-loops
  // maybe make static, see what happens
  //std::vector<QVar> fluxes_(m.num_triangles());

  // iterate over triangles, iterate over edges, calculate new QVar based on 
  // fluxes

  for (auto it = m.edge_begin(); it != m.edge_end(); ++it) {
    // iterate over the triangles adjacent to that edge, at most two.
		auto t1 = (*it).triangle1();
    auto t2 = (*it).triangle2();
    // Ghost triangle
    QVar tmp = t1.value().q_k_bar_;
    Point norm_vec_tmp = (*it).value().norm_vec_;
    if (t2.index() == -1) {
      (*it).value().flux_ = f(norm_vec_tmp.x, norm_vec_tmp.y, dt,
                                tmp, QVar(tmp.h, 0, 0));  
    }
    else {
      (*it).value().flux_ = f(norm_vec_tmp.x, norm_vec_tmp.y, dt, tmp,
                                t2.value().q_k_bar_);
    }
  }

  FluxUpdater<typename MESH::Triangle> fu(dt);
  std::for_each(m.triangle_begin(), m.triangle_end(), fu);

  return t + dt;
}

/** Convert the triangle-averaged values to node-averaged values for viewing. */
template <typename MESH>
void post_process(MESH& m) {
  // HW4B: Post-processing step
  // Translate the triangle-averaged values to node-averaged values
  // Implement Equation 9 from your pseudocode here
  QVar sum(0,0,0);
  double area_sum = 0.0;;
  for (auto it = m.node_begin(); it != m.node_end(); ++it) {
    area_sum = 0.0;
    sum = QVar(0,0,0);
    for (auto it2 = (*it).triangle_begin(); it2 != (*it).triangle_end(); ++it2) {
      sum += (*it2).area() * (*it2).value().q_k_bar_;
      area_sum += (*it2).area();
    }
    sum /= area_sum;
    (*it).value().q = sum;
  }
}

// Simple function double -> double which returns 1 if
// input is negative, otherwise 0.
double H(double x) {
  if (x < 0) {
    return 1;
  }
  else {
    return 0;
  }
}

// The DamBreak functor simulates a dam break and a temporary waterfall as the two
// sides of the dam equalize.  Provides an operator on a position which returns
// a QVar.
struct DamBreak {
  QVar operator()(Point& position) {
    double new_pos = 1.0 + ((0.75)*H(position.x));
    QVar tmp;
    tmp.h = new_pos;
    tmp.hy = 0;
    tmp.hx= 0;
    return tmp;
  }
};

// The Pebble functor simulates the first set of initial conditions specified in
// HW4b, simulating a pebble drop which causes a depression. Provides an operator
// on a position which returns a QVar.
struct Pebble {
  QVar operator() (Point& p) {
    QVar tmp;
    tmp.h = 1 - 0.75 * exp(-80 * ((p.x - 0.75)*(p.x - 0.75) + p.y * p.y));
    tmp.hy = 0;
    tmp.hx = 0;
    return tmp;
  }
};

// The Wave functor simulates one set of initial conditions specified in HW4b, simulating
// the release of a large column of water causing a sharp, fast wave. Provides an operator
// on a position which returns a QVar.
struct Wave {
  QVar operator() (Point& p) {
    QVar tmp;
    tmp.h = 1 + 0.75 * H(-80 * ((p.x - 0.75)*(p.x - 0.75) + p.y * p.y - 0.0225));
    tmp.hy = 0;
    tmp.hx = 0;
    return tmp;
  }
};

// For fun...would use max values of q.h if we wanted to do a data agnostic
// color functor.
struct ColorFunctor {
  template <typename NODE>
  CS207::Color operator()(const NODE& n) {
    if (n.value().q.h/3.0 >= 0 && n.value().q.h/3.0 <= 1){
      return CS207::Color::make_heat(n.value().q.h/3.0);
    }
    else {
       return CS207::Color::make_heat(0.99);
    }
  }
};

int main(int argc, char* argv[])
{
  // Check arguments
  if (argc < 3) {
    std::cerr << "Usage: shallow_water NODES_FILE TRIS_FILE [INIT_CONDITION_TYPE]\n";
    exit(1);
  }

  int type;
  if (argc == 4) {
    type = std::stoi(argv[3]);
  }
  else {
    type = 0;
  }

  MeshType mesh;
  // HW4B: Need node_type before this can be used!
  std::vector<typename MeshType::node_type> mesh_node;

  // Read all Points and add them to the Mesh
  std::ifstream nodes_file(argv[1]);
  Point p;
  while (CS207::getline_parsed(nodes_file, p)) {
    // HW4B: Need to implement add_node before this can be used!
    mesh_node.push_back(mesh.add_node(p));
  }

  // Read all mesh triangles and add them to the Mesh
  std::ifstream tris_file(argv[2]);
  std::array<int,3> t;
  while (CS207::getline_parsed(tris_file, t)) {
    // HW4B: Need to implement add_triangle before this can be used!
    mesh.add_triangle(mesh_node[t[0]], mesh_node[t[1]], mesh_node[t[2]]);
  }

  // Print out the stats
  std::cout << mesh.num_nodes() << " "
            << mesh.num_edges() << " "
            << mesh.num_triangles() << std::endl;

  // HW4B Initialization
  // Set the initial conditions
  // Perform any needed precomputation

  /* Set the initial conditions */
  double max_h = 0.0;
  for (auto it = mesh.node_begin(); it != mesh.node_end(); ++it) { 
    auto n = *it;
    if (type == 0) {
      auto init_cond = DamBreak();
      n.value().q = init_cond(n.position()); 
    }
    if (type == 1) {
      auto init_cond = Pebble();
      n.value().q = init_cond(n.position()); 
    }
    if (type == 2) {
      auto init_cond = Wave();
      n.value().q = init_cond(n.position()); 
    }
    max_h = std::max(max_h, n.value().q.h); 
  }

  // Set the initial values of the triangles to the average of their nodes
  for (auto it = mesh.triangle_begin(); it != mesh.triangle_end(); ++it) {
    auto t = *it; 
      t.value().q_k_bar_.h = (t.node(0).value().q.h + t.node(1).value().q.h + t.node(2).value().q.h) / 3.0; 
  }

  for (auto it = mesh.edge_begin(); it != mesh.edge_end(); ++it) {
    // iterate over the triangles adjacent to that edge, at most two.
    auto t1 = (*it).triangle1();
    Point norm_vec = (*it).normal(t1);
    norm_vec = norm_vec/norm(norm_vec);
    norm_vec *= (*it).length();
    (*it).value().norm_vec_= norm_vec;
  }

  // Launch the SDLViewer
  CS207::SDLViewer viewer;
  viewer.launch();

  auto node_map = viewer.empty_node_map(mesh);
  viewer.add_nodes(mesh.node_begin(), mesh.node_end(),
                   CS207::DefaultColor(), NodePosition(), node_map);
  viewer.add_edges(mesh.edge_begin(), mesh.edge_end(), node_map);
  viewer.center_view();

  double min_edge_length = 100000.0;
  for (auto it = mesh.edge_begin(); it != mesh.edge_end(); ++it) {
    min_edge_length = std::min(min_edge_length, (*it).length());
  }

  // HW4B: Timestep
  // CFL stability condition requires dt <= dx / max|velocity|
  // For the shallow water equations with u = v = 0 initial conditions
  //   we can compute the minimum edge length and maximum original water height
  //   to set the time-step
  // Compute the minimum edge length and maximum water height for computing dt
  double dt = 0.25 * min_edge_length / (sqrt(grav * max_h));

  double t_start = 0;
  double t_end = 0.2;

  // Preconstruct a Flux functor
  EdgeFluxCalculator f;
  ColorFunctor cf;

  // Begin the time stepping
  { jb_parallel::Timer timer("Hyperbolic step");
  for (double t = t_start; t < t_end; t += dt) {

    // Step forward in time with forward Euler
    hyperbolic_step(mesh, f, t, dt);
    // Update node values with triangle-averaged values
    post_process(mesh);

    // Update the viewer with new node positions
    viewer.add_nodes(mesh.node_begin(), mesh.node_end(),
                     cf, NodePosition(), node_map);
    viewer.set_label(t);

    // These lines slow down the animation for small meshes.
		if (mesh.num_nodes() < 100)
    	CS207::sleep(0.03);
  }
  }

  return 0;
}
