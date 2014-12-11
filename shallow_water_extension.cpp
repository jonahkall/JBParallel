/**
 * @file shallow_water.cpp
 * Implementation of a shallow water system using Mesh
 * CS207 HW4
 * Editors: Wenshuai Ye && Yuhao Zhu
 * @brief Reads in two files specified on the command line.
 * First file: 3D point list (one per line) defined by three doubles
 * Second file: Triangles (one per line) defined by 3 indices into the point list
 */

#include <fstream>
#include <cmath>
#include "CS207/SDLViewer.hpp"
#include "CS207/Util.hpp"
#include "CS207/Color.hpp"
#include <iostream>
#include "Point.hpp"
#include "Mesh2.hpp"

// Standard gravity (average gravity at Earth's surface) in meters/sec^2
static constexpr double grav = 9.80665;

typedef double value_type;
typedef unsigned size_type;

/** Water column characteristics */
struct QVar {
  double h;   // Height of column
  double hx;  // Height times average x velocity of column
  double hy;  // Height times average y velocity of column

  /** Default constructor.
   *
   * A default water column is 1 unit high with no velocity. */
  QVar(): h(1), hx(0), hy(0) {
  }
  /** Construct the given water column. */
  QVar(double h_, double hx_, double hy_): h(h_), hx(hx_), hy(hy_) {
  }

  /** Add scalar @a b to this QVar */
  QVar& operator+=(value_type b) {
    h  += b;
    hx += b;
    hy += b;
    return *this;
  }
  /** Subtract scalar @a b from this QVar */
  QVar& operator-=(value_type b) {
    h  -= b;
    hx -= b;
    hy -= b;
    return *this;
  }
  /** Scale this QVar up by scalar @a b */
  QVar& operator*=(value_type b) {
    h  *= b;
    hx *= b;
    hy *= b;
    return *this;
  }
  /** Scale this QVar down by scalar @a b */
  QVar& operator/=(value_type b) {
    h  /= b;
    hx /= b;
    hy /= b;
    return *this;
  }
  /** Add QVar @a b to this QVar */
  QVar& operator+=(const QVar& b) {
    h  += b.h;
    hx += b.hx;
    hy += b.hy;
    return *this;
  }
  /** Subtract QVar @a b from this QVar */
  QVar& operator-=(const QVar& b) {
    h  -= b.h;
    hx -= b.hx;
    hy -= b.hy;
    return *this;
  }
  QVar operator*(value_type b) const {
    return {h*b, hx*b, hy*b};
  }
};


struct Force {
  double Fk;
  double Fm;
  double Ak;
  double Am;
  double rho;
};


struct TriValue{
  double area;
  QVar Q;
};


// HW4B: Placeholder for Mesh Type!
// Define Nodevalue, Edgevalue, Trivalue. Edgevalue is only for placeholder
// or redefine for your particular Mesh
typedef Mesh<QVar, int, TriValue> MeshType;
typedef MeshType::Node Node;
typedef MeshType::Edge Edge;
typedef MeshType::Triangle Triangle;

struct ExternalObj1 {
    double t_;
    ExternalObj1(double t): t_(t) {}
    double operator()(Point p) {
	if ((p.x-3*t_)*(p.x-3*t_)/5 + (p.y)*(p.y)/2 < 0.01)
		return 1000;
	
	return 0.;
    }
    bool isBoat(Point p) {
	return ((p.x-3*t_)*(p.x-3*t_)/5 + (p.y)*(p.y)/2 < 0.01);
    }
};

struct ExternalObj2 {
    double t_;
    ExternalObj2(double t): t_(t) {}
    double operator()(Point p) {
	double ship_center_x = .5 * cos(10*t_) + .2;
	double ship_center_y = .5 * sin(10*t_) + .2;
	if ((p.x-ship_center_x)*(p.x-ship_center_x) + (p.y-ship_center_y)*(p.y-ship_center_y) < 0.01)
		return 500.;
	return 0.;
    }

    bool isBoat(Point p) {
	double ship_center_x = .5 * cos(10*t_) + .2;
	double ship_center_y = .5 * sin(10*t_) + .2;
	return ((p.x-ship_center_x)*(p.x-ship_center_x) + (p.y-ship_center_y)*(p.y-ship_center_y) < 0.01);
    }
};

template <typename F1, typename F2>
struct combine_obj2 {
  F1 f1_;
  F2 f2_;
  combine_obj2(const F1& f1, const F2& f2): f1_(f1), f2_(f2) {
  }
  double operator()(Point p) {
    return f1_(p) + f2_(p);
  }
};  

// Create the combine two force //
template <typename F1, typename F2>
combine_obj2<F1, F2> make_combined_obj(F1 f1, F2 f2){
  return combine_obj2<F1, F2>(f1, f2);
}


/* Oil Color for Shallow Water*/
struct OilColor {
   CS207::Color operator()(const Node& node) {
	float value = 0.;
	for (auto it = node.triangle_begin();it != node.triangle_end();++it){
		if ((*it).value().Q.h > value)
			value = ((*it).value().Q.h-1)*100;
	}
  
   float proportion = (float) exp(value)/(1+exp(value));

   //return CS207::Color::make_heat(proportion); 
   return CS207::Color::make_heat(.5); 
   }
};



/* Color for Ship and water*/
template <typename B1, typename B2>
struct ShipColor {
    B1 b1_;
    B2 b2_;
    CS207::Color operator()(const Node& node) {
	Point p = node.position();       
	if (b1_.isBoat(p)) 
		return CS207::Color::make_heat(.99);
	else if (b2_.isBoat(p))
		return CS207::Color::make_heat(.66);
	
	return CS207::Color::make_heat(.2);
    }

    ShipColor(B1 b1, B2 b2) : b1_(b1), b2_(b2) {}

};

template <typename B1, typename B2>
ShipColor<B1, B2> make_obj_colors(B1 b1, B2 b2){
   return ShipColor<B1, B2>(b1, b2);
}

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
                  const QVar& qk, const QVar& qm, 
		  Force& F) {
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
    //double gh2   = 0.5 * grav * (qm.h*qm.h + qk.h*qk.h) + obj.F * (qm.h+qk.h)/(obj.ro*obj.A);
    double gh2   = 0.5 * grav * (qm.h*qm.h + qk.h*qk.h);
    gh2 += (F.Fk * qk.h / F.Ak + F.Fm * qm.h / F.Am) / F.rho; 

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
    // HW4B: You may change this to plot something other than the
    // positions of the nodes
    return Point(n.position().x, n.position().y, n.value().h);
  }
};

/* Model the surface*/
struct Bathymetry{
   double operator()(Point p){
      return 1;
   }
   double partialx(Point p){
      return 0.;
   }
   double partialy(Point p){
      return 0.;
   }
};



/** Integrate a hyperbolic conservation law defined over the mesh m
 * with flux functor f by dt in time.
 */
template <typename MESH, typename FLUX, typename OBJ>
double hyperbolic_step(MESH& m, FLUX& f, OBJ ext_obj, Bathymetry b, double t, double dt, double rho = 1.0) {
  // Step the finite volume model in time by dt.
  // Pseudocode:
  // Compute all fluxes. (before updating any triangle Q_bars)
  // For each triangle, update Q_bar using the fluxes as in Equation 8.
  // NOTE: Much like symp_euler_step, this may require TWO for-loops

  for (auto it=m.triangle_begin(); it!=m.triangle_end(); ++it) {
    QVar temp_sum = QVar(0.0,0.0,0.0);
    // calculate the center point of the triangle
    Point tri_center(0,0,0);
    Force F;
    F.rho = rho;
    for (size_type i = 1; i <= 3; ++i) {
        tri_center = tri_center + (*it).triangle_node(i).position();
    }
    tri_center = tri_center / 3.;
    QVar surface = QVar(0, -grav*(*it).value().Q.h*b.partialx(tri_center),-grav*(*it).value().Q.h*b.partialy(tri_center));
    F.Fk = ext_obj(tri_center);
    //double Ak = (*it).area();
    F.Ak = 1.;

    for (size_type i=1; i<=3; ++i) {
      
      QVar temp;
      if ((*it).triangle_edge(i).degree() == 1) {
	F.Fm = 0.;
        F.Am = 1.;
        temp = f((*it).normal_vector(i).x, (*it).normal_vector(i).y, dt, (*it).value().Q, QVar((*it).value().Q.h,0.0,0.0), F);
      }
      else {
        Point adj_tri_center(0,0,0);
        for (size_type j = 1; j <= 3; ++j) {
            adj_tri_center = adj_tri_center + (*it).triAdj_triangle(i).triangle_node(j).position();
        }
        adj_tri_center = adj_tri_center / 3.;
        F.Fm = ext_obj(adj_tri_center);
        //double Am = (*it).triAdj_triangle(i).area();
        F.Am = 1.;

        temp = f((*it).normal_vector(i).x, (*it).normal_vector(i).y, dt, (*it).value().Q, (*it).triAdj_triangle(i).value().Q,F);
      }
      temp *= dt;
      temp /= (*it).value().area;
      temp_sum += temp;
    }

    (*it).value().Q -= temp_sum;
    (*it).value().Q += surface*dt;
  }

  return t + dt;
};



/** Convert the triangle-averaged values to node-averaged values for viewing. */
template <typename MESH>
void post_process(MESH& m) {
  // Translate the triangle-averaged values to node-averaged values
  // Implement Equation 9 from your pseudocode here
  for (auto it=m.node_begin(); it!=m.node_end(); ++it) {
    double total_area = 0.0;
    (*it).value() = QVar(0.0, 0.0, 0.0);
    for (auto iit=(*it).triangle_begin(); iit!=(*it).triangle_end(); ++iit) {
      total_area += (*iit).value().area;
      QVar temp = (*iit).value().Q;
      temp *=  (*iit).value().area;
      (*it).value() += temp;
    }
    (*it).value() /= total_area;
  }
};



int main(int argc, char* argv[])
{
  // Check arguments
  if (argc < 3) {
    std::cerr << "Usage: shallow_water NODES_FILE TRIS_FILE\n";
    exit(1);
  }

  MeshType mesh;
  // HW4B: Need node_type before this can be used!
#if 1
  std::vector<typename MeshType::node_type> mesh_node;
#endif

  // Read all Points and add them to the Mesh
  std::ifstream nodes_file(argv[1]);
  Point p;
  while (CS207::getline_parsed(nodes_file, p)) {
    // HW4B: Need to implement add_node before this can be used!
#if 1
    mesh_node.push_back(mesh.add_node(p));
#endif
  }

  // Read all mesh triangles and add them to the Mesh
  std::ifstream tris_file(argv[2]);
  std::array<int,3> t;
  while (CS207::getline_parsed(tris_file, t)) {
    // HW4B: Need to implement add_triangle before this can be used!
#if 1
    mesh.add_triangle(mesh_node[t[0]], mesh_node[t[1]], mesh_node[t[2]]);
#endif
  }

  // Print out the stats
  std::cout << mesh.num_nodes() << " "
            << mesh.num_edges() << " "
            << mesh.num_triangles() << std::endl;

  // HW4B Initialization
  // Set the initial conditions
  // Perform any needed precomputation

// Case 1: Pebble being thrown in and causing a depression
  for (auto it=mesh.node_begin(); it!=mesh.node_end(); ++it) {
    (*it).value() = QVar(0.0,0.0,0.0);
    (*it).value().h = 1.0;
    //(*it).value().h = 1.0-0.75*exp(-80.0*(pow(((*it).position().x-0.75),2)+pow((*it).position().y,2)));
  }

#if 0
// Case 2: Releases a large column of water and causes a very sharp, fast wave
  for (auto it=mesh.node_begin(); it!=mesh.node_end(); ++it) {
    (*it).value() = QVar(0.0,0.0,0.0);
    double temp = ((pow(((*it).position().x-0.75),2)+pow((*it).position().y,2))-pow(0.15,2));
    if (temp < 0)  (*it).value().h = 1.0+0.75;
    else (*it).value().h = 1.0;
  }


// Case 3: A dam break and a temporary waterfall as the two sides attempt to equalize
  for (auto it=mesh.node_begin(); it!=mesh.node_end(); ++it) {
    (*it).value() = QVar(0.0,0.0,0.0);
    if ((*it).position().x < 0) (*it).value().h = 1.0+0.75;
    else (*it).value().h = 1.0;
  }
#endif

  // Set triangle values
  for (auto it=mesh.triangle_begin(); it!=mesh.triangle_end(); ++it) {
    (*it).value().Q = QVar(0.0,0.0,0.0);
    (*it).value().Q += (*it).triangle_node(1).value();
    (*it).value().Q += (*it).triangle_node(2).value();
    (*it).value().Q += (*it).triangle_node(3).value();
    (*it).value().Q /= 3.0;
  }

  // Perform any needed precomputation
  for (auto it=mesh.triangle_begin(); it!=mesh.triangle_end(); ++it) {
    (*it).value().area = (*it).area();
  }


  // Launch the SDLViewer
  CS207::SDLViewer viewer;
  viewer.launch();

  // HW4B: Need to define Mesh::node_type and node/edge iterator
  // before these can be used!
#if 1
  auto node_map = viewer.empty_node_map(mesh);
  viewer.add_nodes(mesh.node_begin(), mesh.node_end(),
                   CS207::DefaultColor(), NodePosition(), node_map);
  viewer.add_edges(mesh.edge_begin(), mesh.edge_end(), node_map);
#endif
  viewer.center_view();

  // HW4B: Timestep
  // CFL stability condition requires dt <= dx / max|velocity|
  // For the shallow water equations with u = v = 0 initial conditions
  //   we can compute the minimum edge length and maximum original water height
  //   to set the time-step
  // Compute the minimum edge length and maximum water height for computing dt
 double min_edge_length = mesh.edge(0).length();
  for (auto it=mesh.edge_begin(); it!=mesh.edge_end(); ++it)
    if ((*it).length() < min_edge_length) min_edge_length = (*it).length();

  double max_height = 0.0;
  for (auto it=mesh.node_begin(); it!=mesh.node_end(); ++it)
    if ((*it).value().h > max_height) max_height = (*it).value().h;

#if 1
  double dt = 0.25 * min_edge_length / (sqrt(grav * max_height));
#else
  // Placeholder!! Delete me when min_edge_length and max_height can be computed!
  //double dt = 0.1;
#endif
  double t_start = 0;
  double t_end = 10;

  // Preconstruct a Flux functor
  EdgeFluxCalculator f;
  //ExternalObj F;
  Bathymetry b;
  // Begin the time stepping
  for (double t = t_start; t < t_end; t += dt) {
    // Step forward in time with forward Euler
    hyperbolic_step(mesh, f, make_combined_obj(ExternalObj1(t),ExternalObj2(t)), b, t, dt, 1000.0);

    // Update node values with triangle-averaged values
    post_process(mesh);

    // Update the viewer with new node positions
    // HW4B: Need to define node_iterators before these can be used!
#if 1
    viewer.add_nodes(mesh.node_begin(), mesh.node_end(),
                     make_obj_colors(ExternalObj1(t),ExternalObj2(t)), NodePosition(), node_map);

    /*viewer.add_nodes(mesh.node_begin(), mesh.node_end(),
                     OilColor(), NodePosition(), node_map);*/
#endif
    viewer.set_label(t);

    // These lines slow down the animation for small meshes.
    // Feel free to remove them or tweak the constants.
    if (mesh.num_nodes() < 100)
      CS207::sleep(0.05);
  }

  return 0;
}
