#pragma once


#include <fstream>
#include <algorithm>

#include "CS207/SDLViewer.hpp"
#include "CS207/Util.hpp"
#include "CS207/Color.hpp"

#include "Graph.hpp"
#include "Point.hpp"

using namespace std;

/** @file Mesh.hpp
 * @brief A Mesh is composed of nodes, edges, and triangles such that:
 *  -- All triangles have three nodes and three edges.
 *  -- All edges belong to at least one triangle and at most two triangles.
 */

/** @class Mesh
 * @brief A template for 3D triangular meshes.
 *
 * Users can add triangles and retrieve nodes, edges, and triangles.
 */
template <typename N, typename E, typename T>
class Mesh {

public:
  // Public typedefs.
	typedef long size_type;
  typedef N node_value_type;
  typedef E edge_value_type;
  typedef T triangle_value_type;

private:

  /* Internal structs store templated value type (exposed to the user) 
   * and extra information allowing us to write iterators.
   **/
  struct internal_node {
    node_value_type nvt;
    vector<size_type> adj_tris_to_node;
  };

  struct internal_edge {
    edge_value_type evt;
    vector<size_type> adj_tris_to_edge;
  };

  struct triangle_info {
    triangle_value_type tvt_;
    // the three possible triangles that could be adjacent
    vector<size_type> adj_tris_to_tri;
    vector<size_type> nodes_;
  };

public:
  typedef Graph<internal_node, internal_edge> GraphType;

  typedef typename GraphType::Node DepNode;
  typedef typename GraphType::Edge DepEdge;
  typedef typename GraphType::EdgeIterator DepEdgeIter;

  // Forward declarations.
  class MeshNode;
  typedef MeshNode node_type;

  class MeshEdge;
  typedef MeshEdge edge_type;
	
	class Triangle;
	class nodes_triangle_iterator;


  typedef long uid_type;

 private:

  GraphType g_;
  vector<triangle_info> tris_;

 public:

  //////////////////
  /// MESH NODES ///
  //////////////////

  /** @class Mesh::MeshNode
   * @brief Class representing the mesh's nodes.
   *
   * MeshNode objects are used to access information about the Mesh's nodes.
   * Most operations simply forward to the underlying graph.  We simply store
   * the index of the node in the underlying graph corresponding to the MeshNode,
   * allowing us to exploit the same proxy pattern as graph.
   */
  class MeshNode : private totally_ordered<MeshNode> {

		private:
			Mesh* m_;
			size_type i_;
	    MeshNode(const Mesh* m, size_type i) :
	      m_(const_cast<Mesh*>(m)), i_(i)  {
	    }

      friend class Mesh;

			DepNode node(){
				return m_->g_.node(i_);
			}

		public:
			// All functions use GraphType::Node functions
			// All run-times dependent on graph run-times
			
      // Construct an invalid node.
      MeshNode() {
        m_ = NULL;
        i_ = -1;
      }

			const Point& position() const {
				return  m_->g_.node(i_).position();
			}
			
      /* Return the position of this MeshNode */
			Point& position() {
				return m_->g_.node(i_).position();
			}
			
      /* Return the value of a MeshNode (uses the value element of
       * internal_node) 
       **/
			node_value_type& value(){
				return  m_->g_.node(i_).value().nvt;
			}

      const node_value_type& value() const {
        return  m_->g_.node(i_).value().nvt;
      }
			
			bool operator==(const MeshNode& x) const {
				return  m_->g_.node(i_) == x.node();
			}

			bool operator<(const MeshNode& x) const {
				return  m_->g_.node(i_) < x.m_->g_.node(x.i_);
			}

			size_type index() const{
				return i_;
			}
			
      // Begin iterator for incident triangles to a MeshNode.
			nodes_triangle_iterator triangle_begin(){
				return nodes_triangle_iterator(m_, m_->g_.node(i_).value().adj_tris_to_node.begin());
			}

      // End iterator for incident triangles to a MeshNode
			nodes_triangle_iterator triangle_end(){
				return nodes_triangle_iterator(m_, m_->g_.node(i_).value().adj_tris_to_node.end());
			}
  };


  ////////////////
  // MESH EDGES //
  ////////////////

  /** @class Mesh::MeshEdge
   * @brief Class representing the mesh's edges.
   *
   * Edges are order-insensitive pairs of nodes. Two Edges with the same nodes
   * are considered equal if they connect the same nodes, in either order.
   */
  class MeshEdge : private totally_ordered<MeshEdge>{
    private:
      Mesh* m_;
			size_type index1_;
			size_type index2_;			
			MeshEdge(const Mesh* m, size_type index1, size_type index2) : m_(const_cast<Mesh*>(m)),
			index1_(index1), index2_(index2) {}
      friend class Mesh;

    public: 

      double length() {
        auto e = m_->g_.add_edge(MeshNode(m_, index1_).node(), MeshNode(m_, index2_).node());
        return e.length();
      }

      /** Return a node of this Edge */
      MeshNode node1() const {
        return MeshNode(m_, index1_);
      }

      /** Return the other node of this Edge */
      MeshNode node2() const {
        return MeshNode(m_, index2_);
      }

      bool operator==(const MeshEdge& x) const {
        if (m_ != x.m_){
          return false;
        }
        else {
          return (std::min(x.index1_, x.index2_) ==
              std::min(index1_, index2_) && std::max(x.index1_, index1_)
                  == std::max(index1_, index2_));
        }
        return true;
      }

      bool operator<(const MeshEdge& x) const {
        if (g_ < x.g_) {
          return true;
      }

      size_type x_1 = std::min(x.index1_, x.index2_);
      size_type x_2 = std::max(x.index1_, x.index2_);

      // Compare the first node ids.
      if (std::min(index1_, index2_) != x_1) {
        return std::min(index1_, index2_) < x_1;
      }
      // If those are equal, compare the second node ids.
      else if (std::max(index1_, index2_) != x_2) {
        return std::max(index1_, index2_) < x_2;
      }
      // Both nodes have equal ids.
      else {
        return false;
      }
    }

    edge_value_type& value() {
      // can't construct by uids, so have to use add_edge here, which
      // is hack-y but only way to make it work.
      // also realistically this function is never going to get called....lol.
      auto e = m_->g_.add_edge(MeshNode(m_, index1_).node(), MeshNode(m_, index2_).node());
      return e.value().evt;
    }	
		
		/* Returns the first triangle associated with an edge
		* comlexity is dependant on complexity of Graph::Edge.value(),
		* which is currently O(d)
		*/
		Triangle triangle1(){
      auto e = m_->g_.add_edge(MeshNode(m_, index1_).node(), MeshNode(m_, index2_).node());
			auto index = e.value().adj_tris_to_edge[0];
			auto ti = m_->tris_[index];
			return Triangle(m_, index, ti.nodes_[0], ti.nodes_[1], ti.nodes_[2]);
		}
		
		/* Returns the second triangle associated with an edge
		* if it exists, or else an invalid triangle given by the
		* default triangle constructor. 
		* Comlexity is dependant on complexity of Graph::Edge.value(),
		* which is currently O(d)
		*/
		Triangle triangle2(){
      auto e = m_->g_.add_edge(MeshNode(m_, index1_).node(), MeshNode(m_, index2_).node());
			if (e.value().adj_tris_to_edge.size() == 2){
				auto index = e.value().adj_tris_to_edge[1];
				auto ti = m_->tris_[index];
				return Triangle(m_, index, ti.nodes_[0], ti.nodes_[1], ti.nodes_[2]);
			}
			else{
				return Triangle();
			}
		}

		/** Returns the outward normal vector from the edge with
		* respect to a triangle. 
		*/
    const Point normal(Triangle& t) const {
			// sets out_index to the node in Triangle t that is not 
			// in the edge
      int out_index = 0;
      for (int i = 0; i < 3; ++i) {
        if (t.node(i).index() != index1_ && t.node(i).index() != index2_) {
          out_index = i;
          break;
        }
      }

     
      // Caculates normal vector by computing the cross product between our edge 
			// and an adjacent edge in the triangle and then crossing that resultant
			// vector with our edge once again
			Point normal_vector = 
        cross(cross((m_->g_.node(index1_).position() - m_->g_.node(index2_).position()),
        (t.node(out_index).position() - m_->g_.node(index2_).position())),
        (m_->g_.node(index1_).position() - m_->g_.node(index2_).position()));

      return -1*normal_vector;
    }

  };

  //////////////////////
 /// MESH TRIANGLES ///
//////////////////////

 /** @class Mesh::MeshTriangle
  * @brief Class representing the mesh's triangles.
  *
  * Triangles are order-insensitive triplets of nodes. Two triangles with the same nodes
  * are considered equal if they connect the same nodes, in either order.
	* Since they are a 2-d geometric shape, we can compute their area
  */
  class Triangle : private totally_ordered<Triangle>{
	private:
    Mesh* m_;
		uid_type uid_;

    // indices of the three nodes belonging to the triangle
    size_type tid_1_;
    size_type tid_2_;
    size_type tid_3_;

    friend class Mesh;

    /* Triangle Private Constructor ensures numeric ordering 
		* of node id's so that same triangle is constructed 
		* regardless of order of user input
		*/
    Triangle(const Mesh* m, uid_type uid, size_type tid1, size_type tid2, size_type tid3) 
        : m_(const_cast<Mesh*>(m)), uid_(uid) {
      tid_1_ = std::min(std::min(tid1,tid2),tid3);
      tid_3_ = std::max(std::max(tid1,tid2),tid3);
      if (tid1 != tid_1_ && tid1 != tid_3_) {
        tid_2_ = tid1;
      }
      if (tid2 != tid_1_ && tid2 != tid_3_) {
        tid_2_ = tid2;
      }
      if (tid3 != tid_1_ && tid3 != tid_3_) {
        tid_2_ = tid3;
      }
    };

	public:
		// Returns index of triangle
    size_type index () {
      return uid_;
    }

    bool operator ==(const Triangle& t2) const {
      return (t2.m_ == m_ && uid_ == t2.uid_);
    }

		// Construct an invalid triangle
		Triangle() :
			m_(nullptr),uid_(-1),tid_1_(-1),tid_2_(-1), tid_3_(-1) {};
		
		/* Returns the ith node of a triangle
	 	 *  @[pre] 0<= i < 3	
		 */
    MeshNode node(int i) {
      return MeshNode(m_, m_->tris_[uid_].nodes_[i]);
    }
		/** Returns the edges of a triangle. O(1) */
		MeshEdge edge1(){
			return MeshEdge(m_, m_->tris_[uid_].nodes_[0], m_->tris_[uid_].nodes_[1]);
		}
		MeshEdge edge2(){
			return MeshEdge(m_, m_->tris_[uid_].nodes_[1], m_->tris_[uid_].nodes_[2]);
		}
		MeshEdge edge3(){
			return MeshEdge(m_, m_->tris_[uid_].nodes_[0], m_->tris_[uid_].nodes_[2]);
		}
		
    /** Computes area of a triangle. O(1) */
    double area () {
      Point p1 = this->node(0).position();
      Point p2 = this->node(1).position();
      Point p3 = this->node(2).position();

      Point v1 = p2-p1;
      Point v2 = p3-p1;

      return norm(cross(v1,v2))/2.0;
    }

		/** Returns value of triangle */
    triangle_value_type& value () {
      return m_->tris_[uid_].tvt_;
    }

  };

  /** Return the number of nodes in the mesh. */
  size_type num_nodes() const {
    return g_.num_nodes();
  }

  /** Return the number of edges in the mesh. */
  size_type num_edges() const {
    return g_.num_edges();
  }

  /** Return the number of triangles in the mesh. */
  size_type num_triangles() const {
    return tris_.size();
  }

  MeshNode add_node(const Point& p) {
    auto n = g_.add_node(p);
    return MeshNode(this, n.index());
  }

 	 /**
	 * Creates a new triangle, adds it to tris_, and updates all of the relevant
   * information in node, edge, and triangle info structs 
   * so that all iterators are valid. 
	 * @[pre] @a m1_, @a m_2, @a m_3 are valid MeshNodes and are not colinear
	 * @[pre] there does not already exist two other triangles that share an edge
	 * between any two of @a m1_, @a m2_, @a m3_
	 * @[post] allows to iterate over adjacent triangles to a triangle,
	 * MeshEdge, or MeshNode
	 * @returns newly constructed triangle
	 * Complexity currently limited by edge functions, which is O(d)
	 */
  const Triangle add_triangle(MeshNode& m1_, MeshNode& m2_, MeshNode& m3_) {
		// computes the index of the triangle to be added
    size_type new_uid = tris_.size();
		// adds edges to graph
    auto e1_ = g_.add_edge(m1_.node(), m2_.node());
    auto e2_ = g_.add_edge(m1_.node(), m3_.node());
    auto e3_ = g_.add_edge(m2_.node(), m3_.node());
    triangle_info ti;
		// updates adjacent triangles to nodes to include new triangle
    m1_.node().value().adj_tris_to_node.push_back(new_uid);
    m2_.node().value().adj_tris_to_node.push_back(new_uid);
    m3_.node().value().adj_tris_to_node.push_back(new_uid);
		// adds nodes to internal triangle
		ti.nodes_.push_back(m1_.index());
		ti.nodes_.push_back(m2_.index());
		ti.nodes_.push_back(m3_.index());
		// if edge already has two triangles
    if (e1_.value().adj_tris_to_edge.size() >= 2 ||
        e2_.value().adj_tris_to_edge.size() >= 2 ||
        e3_.value().adj_tris_to_edge.size() >= 2) {
      // the user is trying to add a third triangle to an edge, don't let them.
      assert(false);
    }

		// if an edge already is part of a triangle
		// add our new triangle to the adjacent triangle list of the old triangle
		// and vice-versa
    if (e1_.value().adj_tris_to_edge.size() == 1) {
      ti.adj_tris_to_tri.push_back(e1_.value().adj_tris_to_edge[0]);
      this->tris_[e1_.value().adj_tris_to_edge[0]].adj_tris_to_tri.push_back(new_uid);
    }
    if (e2_.value().adj_tris_to_edge.size() == 1) {
      ti.adj_tris_to_tri.push_back(e2_.value().adj_tris_to_edge[0]);
      this->tris_[e2_.value().adj_tris_to_edge[0]].adj_tris_to_tri.push_back(new_uid);
    }
    if (e3_.value().adj_tris_to_edge.size() == 1) {
      ti.adj_tris_to_tri.push_back(e3_.value().adj_tris_to_edge[0]);
      this->tris_[e3_.value().adj_tris_to_edge[0]].adj_tris_to_tri.push_back(new_uid);
    }
		// add new triangle's index to list of adjacent triangles
		// per edge
    e1_.value().adj_tris_to_edge.push_back(new_uid);
    e2_.value().adj_tris_to_edge.push_back(new_uid);    
    e3_.value().adj_tris_to_edge.push_back(new_uid);
		// push_back triangle info to triangle_info list    		
    tris_.push_back(ti);
		// return triangle
    return Triangle(this, new_uid, m1_.index(), m2_.index(), m3_.index());
  }

 ///////////////
 // Iterators //
 ///////////////

 /** @class Mesh::mesh_node_iterator
  * @brief Iterator class for MeshNodes. A two-way iterator. 
  * Really just a wrapper for an index: since the nodes are stored in a vector,
  * iterating is easy.
  */
  class mesh_node_iterator : private totally_ordered<mesh_node_iterator> {
   public:
    // These type definitions help us use STL's iterator_traits.
    /** Element type. */
    typedef MeshNode value_type;
    /** Type of pointers to elements. */
    typedef MeshNode* pointer;
    /** Type of references to elements. */
    typedef MeshNode& reference;
    /** Iterator category. */
    typedef std::input_iterator_tag iterator_category;
    /** Difference between iterators */
    typedef std::ptrdiff_t difference_type;

    /** Construct an invalid mesh_node_iterator. */
    mesh_node_iterator() {
      m_ = NULL;
      node_index_ = -1;
    }

    // Returns the current index into the node_info vector.
    size_type index () const {
      return node_index_;
    }

    // Dereference.  Returns a MeshNode with the current index on the fly by
    // forwarding to graph.
    MeshNode operator*() const {
      return MeshNode(m_, node_index_);
    }

    // Increment the mesh_node_iterator.  Simply increments the current node index.
    // Returns the current node iterator after increment.
    mesh_node_iterator& operator++() {
      ++node_index_;
      return *this;
    }

    mesh_node_iterator& operator--() {
      assert(node_index_ > 0);
      --node_index_;
      return *this;
    }

    // Checks if two node iterators are equal by checking that they have the
    // same graph pointer and the same node index.
    bool operator==(const mesh_node_iterator& mni) const {
      return (mni.index() == node_index_ && mni.m_ == m_);
    }

   private:
    friend class Mesh;
    size_type node_index_;
    Mesh* m_;
    mesh_node_iterator(const Mesh* m, size_type index = 0) :
        node_index_(index), m_(const_cast<Mesh*>(m)){}
  };

  // Returns a node_iterator indexed to the first node in the graph.
  mesh_node_iterator node_begin() const {
    return mesh_node_iterator(this, 0);
  }

  // Returns a node_iterator indexed to the last node in the graph.
  mesh_node_iterator node_end() const {
    return mesh_node_iterator(this, g_.num_nodes());
  }

 /** @class Mesh::mesh_edge_iterator
  * @brief Iterator class for MeshEdges. A forward iterator. 
  * Implemented by simply wrapping an edge_iterator in the underlying graph.
  */
  class mesh_edge_iterator : private totally_ordered<mesh_edge_iterator> {
    public:
      // These type definitions help us use STL's iterator_traits.
      /** Element type. */
      typedef MeshEdge value_type;
      /** Type of pointers to elements. */
      typedef MeshEdge* pointer;
      /** Type of references to elements. */
      typedef MeshEdge& reference;
      /** Iterator category. */
      typedef std::input_iterator_tag iterator_category;
      /** Difference between iterators */
      typedef std::ptrdiff_t difference_type;

      /** Construct an invalid EdgeIterator. */
      mesh_edge_iterator() {
        dei_ = DepEdgeIter();
      }

      MeshEdge operator*() const {
        auto e = *dei_;
        return MeshEdge(m_, e.node1().index(), e.node2().index());
      }

      mesh_edge_iterator& operator++() {
        ++dei_;
        return *this;
      }

      bool operator==(const mesh_edge_iterator& mei) const {
        return (dei_ == mei.dei_ && m_ == mei.m_);
      }
    
    private:
      friend class Mesh;
      DepEdgeIter dei_;
      Mesh* m_;

      // Private constructor.
      mesh_edge_iterator(DepEdgeIter dei, const Mesh* m)
          : dei_(dei), m_(const_cast<Mesh*>(m)) {};
  };

  // Returns an iterator to the first MeshEdge in the graph.
  mesh_edge_iterator edge_begin() const {
    return mesh_edge_iterator(g_.edge_begin(), this);
  };

  // Returns a mesh_edge_iterator to the end of the mesh's edges in
  // adjacency list format.
  mesh_edge_iterator edge_end() const {
    return mesh_edge_iterator(g_.edge_end(), this);
  };
	

 /** @class Mesh::mesh_triangle_iterator
  * @brief Iterator class for Triangles. A forward iterator. */
  class mesh_triangle_iterator : private totally_ordered<mesh_triangle_iterator> {
	public:
 
    // These type definitions help us use STL's iterator_traits.
    /** Element type. */
    typedef Triangle value_type;
    /** Type of pointers to elements. */
    typedef Triangle* pointer;
    /** Type of references to elements. */
    typedef Triangle& reference;
    /** Iterator category. */
    typedef std::input_iterator_tag iterator_category;
    /** Difference between iterators */
    typedef std::ptrdiff_t difference_type;

    /** Construct an invalid EdgeIterator. */
    mesh_triangle_iterator() {
			m_ = nullptr;
			index_ = -1;
    }

    // Dereference overload, returns the triangle from our current index.
    Triangle operator*() const {
      auto t = m_->tris_[index_];
      return Triangle(m_, index_, t.nodes_[0], t.nodes_[1], t.nodes_[2]);
    }

    /** Operator++ overload, simply increments the index we are wrapping
     *  around.
     */
    mesh_triangle_iterator& operator++() {
      ++index_;
      return *this;
    }

    bool operator==(const mesh_triangle_iterator& mti) const {
      return (index_ == mti.index_ && m_ == mti.m_);
    }

    mesh_triangle_iterator& operator[](int n){
      assert(n < m_->tris_.size() && n >= 0);
      index_ = n;
      return *this;
    }
    
    int operator-(mesh_triangle_iterator a){
      return index_ - a.index_;
    }

    mesh_triangle_iterator& operator+=(int n){
      index_ += n;
      return *this;
    }
    
    mesh_triangle_iterator& operator-=(int n){
      index_ -= n;
      return *this;
    }
  
  private:
    friend class Mesh;
    size_type index_;
    Mesh* m_;

    // Private constructor
    mesh_triangle_iterator(size_type index, const Mesh* m) :
        index_(index), m_(const_cast<Mesh*>(m)) {};
	};

  friend mesh_triangle_iterator operator+(mesh_triangle_iterator a, int n){
    return a+=n;
  }
  
  friend mesh_triangle_iterator operator+(int n, mesh_triangle_iterator a){
    return a+=n;
  }
  
  friend mesh_triangle_iterator operator-(mesh_triangle_iterator a, int n){
    return a-=n;
  }
	
  // Returns an iterator to the first triangle in the graph.
  mesh_triangle_iterator triangle_begin() const {
    return mesh_triangle_iterator(0, this);
  };

  // Returns a mesh_triangle_iterator to the end of the graph's triangles
  mesh_triangle_iterator triangle_end() const {
    return mesh_triangle_iterator(tris_.size(), this);
  };
	
 /** @class Mesh::nodes_triangle_iterator
  * @brief Iterator class for incident Triangles of MeshNode.
	* A forward iterator. Uses underlying vector iterator of tris_.
  * Nodes provide triangle_begin and triangle_end member functions.
	*/
  class nodes_triangle_iterator : private totally_ordered<nodes_triangle_iterator> {
	public:
 
    // These type definitions help us use STL's iterator_traits.
    /** Element type. */
    typedef Triangle value_type;
    /** Type of pointers to elements. */
    typedef Triangle* pointer;
    /** Type of references to elements. */
    typedef Triangle& reference;
    /** Iterator category. */
    typedef std::input_iterator_tag iterator_category;
    /** Difference between iterators */
    typedef std::ptrdiff_t difference_type;

    /** Construct an invalid EdgeIterator. */
    nodes_triangle_iterator() {
			m_ = nullptr;
			it_ = nullptr;
    }

    Triangle operator*() const {
      auto t = m_->tris_[*it_];
      return Triangle(m_, *it_, t.nodes_[0], t.nodes_[1], t.nodes_[2]);
    }

    nodes_triangle_iterator& operator++() {
      ++it_;
      return *this;
    }

    // == operator overload for nodes' triangle iterator
    bool operator==(const nodes_triangle_iterator& mti) const {
      return (it_ == mti.it_ && m_ == mti.m_);
    }
  
  private:
    friend class Mesh;
    vector<size_type>::iterator it_;
    Mesh* m_;

    // Private constructor
    nodes_triangle_iterator(const Mesh* m, vector<size_type>::iterator it)
        : it_(it), m_(const_cast<Mesh*>(m)) {};
	};
};

