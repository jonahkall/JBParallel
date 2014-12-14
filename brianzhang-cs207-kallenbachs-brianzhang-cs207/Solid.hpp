// Look at Interactables.hpp instead

#ifndef CS207_SOLID_HPP
#define CS207_SOLID_HPP


/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <vector>
#include <cassert>
#include <omp.h>


#include "CS207/Util.hpp"
#include "Point.hpp"
#include "Graph.hpp"
#include "../jb_parallel.hpp"

typedef Graph<double, double> GraphType;
using namespace jb_parallel;

class Interactables {
public:
	//using GraphType = Graph<Point, double>;

	template<typename Node>
	struct position_mod {
		double dt_;
		void operator ()(Node n) {
			n.position() += n.value() * dt_;
		}
		position_mod(double dt) : dt_(dt) {};
	};

	template<typename Node>
	struct fixed_mod {
		Graph<Point,double>* g_;
		void operator ()(Node n) {
			if (g_->has_node(n))
	    	n.value() = Point(0,0,0);
	  }
	  fixed_mod(Graph<Point, double>* g) : g_(g) {};
	};

	struct vel_mod {
		double& grav;
		double& t;
		double& dt;
		double& K_;
		double& mass_;
		double& alpha_;
		void operator () (Graph<Point,double>::Node n) {
			// Mass-spring force
			Point f = Point(0,0,0);
			for (auto it2 = n.edge_begin(); it2 != n.edge_end(); ++it2) {
				auto e = *it2;
				auto n2 = e.node2();
				double L = e.value();
				f += -K_ * (n.position() - n2.position()) / e.length() * (e.length() - L);
			}

			// Gravity force
			f += Point(0,0,-mass_*grav);

			// Damping force
			f -= alpha_*n.value();

			n.value() += f * (dt / mass_);
		}
		vel_mod(double gravi, double ti, double dti, double K_i, double mass_i, double alpha_i) :
				grav(gravi), t(ti), dt(dti), K_(K_i), mass_(mass_i), alpha_(alpha_i) {};
	};

	template<typename Ball, typename Sheet, typename Node>
	struct interact_mod {
		Ball b;
		Sheet s;
		double grav;
		double t;
		double dt;
		void operator () (Node n) {
			if (norm(n.position() - b.c_) < b.r_) {
		    Point R = n.position() - b.c_;
		    assert(norm(R) > 0);
		    R /= norm(R);
		    // R is a unit vector pointing from the center of the sphere towards our point
		    // f represents the total force that the point in the sheet experiences
		    Point f = Point(0,0,0);

		    // Mass-spring force
		    for (auto it2 = n.edge_begin(); it2 != n.edge_end(); ++it2) {
		      auto e = *it2;
		      auto n2 = e.node2();
		      double L = e.value();
		      f += -s.K_ * (n.position() - n2.position()) / e.length() * (e.length() - L);
		    }

		    Point f2 = R*dot(f, R);

		    n.position() = b.c_ + R*b.r_;
		    // n.value() += f * (dt / mass_);
		    b.v_ += f2 * (dt / b.m_);	// reaction force to the force on the point
		    n.value() -= R*dot(n.value(), R);
		    n.value() += R*dot(b.v_, R);
			}
	}
	interact_mod(Ball bi, Sheet si, double gravi, double ti, double dti) :
		b(bi), s(si), grav(gravi), t(ti), dt(dti) {};
};
	class Ball {
	
	using GraphType = Graph<double, double>;

	public:
		void draw(CS207::SDLViewer& viewer) {

		  auto node_map = viewer.empty_node_map(g_);
		  viewer.add_nodes(g_.node_begin(), g_.node_end(), node_map);
		  viewer.add_edges(g_.edge_begin(), g_.edge_end(), node_map);
			return;
		}

		double step(double grav, double t, double dt) {
			// update position
			translate(v_ * dt);
			
			// update velocity
		  v_ += Point(0,0,-grav*dt);

			return t + dt;
		}

		Ball(Point c = Point(0,0,0), double r = 1.0, Point v = Point(0,0,0), double m=1.0)
		{
		  std::string NODE_PATH = "data/sphere1082.nodes";
		  std::string TRIS_PATH = "data/sphere1082.tris";
 			
		  GraphType ball;
		  Point p;
		  std::vector<typename GraphType::node_type> nodes;

		  std::ifstream nodes_file(NODE_PATH);
		  std::ifstream tris_file(TRIS_PATH);

		  while (CS207::getline_parsed(nodes_file, p)) {
		    nodes.push_back(ball.add_node(p*r+c));
		  }

		  // Interpret each line of the tris_file as three ints which refer to nodes
		  std::array<int,3> t;
		  while (CS207::getline_parsed(tris_file, t))
		    for (unsigned i = 1; i < t.size(); ++i)
		      for (unsigned j = 0; j < i; ++j)
		        ball.add_edge(nodes[t[i]], nodes[t[j]]);

		  g_ = ball;
		  c_ = c;
		  r_ = r;
		  v_ = v;
		  m_ = m;
		}

		~Ball()
		{
			g_.clear();
		}

	private:
		void translate(Point vec)
		{
			c_ += vec;
		  auto lam = [&](GraphType::Node n) {n.position()+=vec;};
		  jb_parallel::for_each(g_.node_begin(), g_.node_end(), lam);
		}

		Point c_;		// center
		double r_;	// radius
		Point v_;		// velocity
		double m_;	// mass
		GraphType g_;	// graph used to draw the ball, gets translated with the ball

		friend class Interactables;
	};

	class Floor {

	using GraphType = Graph<double, double>;

	public:
		void draw(CS207::SDLViewer& viewer) {
			auto node_map = viewer.empty_node_map(g_);
		  viewer.add_nodes(g_.node_begin(), g_.node_end(), node_map);
		  viewer.add_edges(g_.edge_begin(), g_.edge_end(), node_map);
			return;
		}

		// Constructor
		Floor(double z, double ratio )
		{
			std::string NODE_PATH = "data/grid2.nodes";
		  std::string TETS_PATH = "data/grid2.tets";
				
		  GraphType graph;
		  Point p;
		  std::vector<typename GraphType::node_type> nodes;

		  std::ifstream nodes_file(NODE_PATH);
		  std::ifstream tets_file(TETS_PATH);

		  while (CS207::getline_parsed(nodes_file, p)) {
		    nodes.push_back(graph.add_node((p - Point(0.5,0.5,0))*5 + Point(0,0,z) ));
		  }

		  // Interpret each line of the tets_file as three ints which refer to nodes
		  std::array<int,4> t;
		  while (CS207::getline_parsed(tets_file, t))
		    for (unsigned i = 1; i < t.size(); ++i) {

		      graph.add_edge(nodes[t[0]], nodes[t[1]]);
		      graph.add_edge(nodes[t[0]], nodes[t[2]]);
		      graph.add_edge(nodes[t[1]], nodes[t[3]]);
		      graph.add_edge(nodes[t[2]], nodes[t[3]]);
		 
		    }

		  g_ = graph;
		  z_ = z;
		  ratio_ = ratio;

		}

	private:
		double z_;
		double ratio_;
		GraphType g_;

		friend class Interactables;
		
	};

	class Sheet {

	using GraphType = Graph<Point, double>;
	typedef GraphType::Node Node;
	typedef GraphType::Edge Edge;
	// Each Node stores a velocity; each Edge stores a rest length

	public:
		void draw(CS207::SDLViewer& viewer) {

		  auto node_map = viewer.empty_node_map(g_);
		  viewer.add_nodes(g_.node_begin(), g_.node_end(), node_map);
		  viewer.add_edges(g_.edge_begin(), g_.edge_end(), node_map);
			return;
		}

		double step(double grav, double t, double dt) {
			// update position
		  position_mod<GraphType::Node> pm(dt);
		  for_each(g_.node_begin(), g_.node_end(), pm);

		  /*for (auto it = g_.node_begin(); it != g_.node_end(); ++it) {
		    auto n = *it;
		    n.position() += n.value() * dt;
		  }*/

			// update velocity
			vel_mod vm(grav,t,dt,K_,mass_,alpha_);
			for_each(g_.node_begin(), g_.node_end(), vm);
			// Forces = gravity, mass-spring, damping, node constraints
		  /*for (auto it = g_.node_begin(); it != g_.node_end(); ++it) {
		    auto n = *it;

		    // Mass-spring force
		    Point f = Point(0,0,0);
		    for (GraphType::incident_iterator it2 = n.edge_begin(); it2 != n.edge_end(); ++it2) {
		      GraphType::Edge e = *it2;
		      GraphType::Node n2 = e.node2();
		      double L = e.value();
		      f += -K_ * (n.position() - n2.position()) / e.length() * (e.length() - L);
		    }

		    // Gravity force
		    f += Point(0,0,-mass_*grav);

		    // Damping force
		    f -= alpha_*n.value();

		    n.value() += f * (dt / mass_);
		  }*/

		  // apply fixed node constraint
		  fixed_mod<GraphType::Node> fm(&g_);
		 	jb_parallel::for_each(fixed_nodes_.begin(), fixed_nodes_.end(), fm);
	// #pragma omp parallel for
	// 	  for (auto it = fixed_nodes_.begin(); it != fixed_nodes_.end(); ++it) {
	// 	  	auto n = *it;
	//       if (g_.has_node(n))
	//         n.value() = Point(0,0,0);
	// 	  }

			return t + dt;
		}

		//constructor
		Sheet(GraphType& g, double K, double alpha, double mass, 
								std::vector<GraphType::Node> fixed_nodes) :
			g_(g), K_(K), alpha_(alpha), mass_(mass), fixed_nodes_(fixed_nodes)
				{};

		~Sheet()
		{
			g_.clear();
		}
	private:
		GraphType& g_; // sheet coordinates are stored as a mass-spring graph
		double K_; // spring const
		double alpha_; // damping const
		double mass_; // mass of each node

		std::vector<GraphType::Node> fixed_nodes_;

		friend class Interactables;
		
	};

	using GraphType = Graph<Point, double>;

	void interact(Ball& b, Sheet& s, double grav, double t, double dt) {
		interact_mod<Ball, Sheet, GraphType::Node> im(b,s,grav,t,dt);
		for_each(s.g_.node_begin(), s.g_.node_end(), im);
  }

	void interact(Ball& b, Floor& f) {
		double diff = b.c_.z - f.z_;

		if (diff  < b.r_)
		{
			b.v_.z *= -f.ratio_;
			b.translate(Point(0,0, b.r_ - diff));
		}
	}
};

#endif