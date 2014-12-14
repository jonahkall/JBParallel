#include "Graph.hpp"
#include "CS207/Util.hpp"
#include "Point.hpp"

#include <climits>
#include <unordered_map>

#pragma once
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

 private:
  // Internal structs
  struct internal_node;
  struct internal_edge;
  struct internal_triangle;

 public:
  // Typedef everything
  typedef unsigned size_type;

  typedef Point point_type;

  typedef N node_value_type;
  typedef E edge_value_type;
  typedef T triangle_value_type;
  typedef Mesh<N, E, T> mesh_type;

  typedef Graph<internal_node, internal_edge> NodeGraph;
  typedef Graph<internal_triangle> TriangleGraph;

  typedef typename NodeGraph::node_type node_type;
  typedef typename NodeGraph::edge_type edge_type;
  typedef typename NodeGraph::node_iterator node_iterator;
  typedef typename NodeGraph::edge_iterator edge_iterator;

  class Triangle;
  typedef Triangle triangle_type;

  class TriangleIterator;
  typedef TriangleIterator triangle_iterator;

  class NeighborIterator;
  typedef NeighborIterator neighbor_iterator;

  typedef typename TriangleGraph::node_type t_node_type;
  typedef typename TriangleGraph::node_iterator t_node_iterator;

  /** Construct an empty mesh. 
   * Complexity: O(1).
   */
  Mesh() :
    n_graph_(NodeGraph()),
    t_graph_(TriangleGraph()),
    internal_triangles_(std::unordered_map<size_type, internal_triangle>())
  {}

  /** Remove all nodes, edges, and triangles from this mesh.
   * @post num_nodes() == 0 && num_edges() == 0 && num_triangles() == 0
   *
   * Invalidates all outstanding Node, Edge, and Triangle objects.
   * Complexity: O(1).
   */
  void clear() {
    n_graph_.clear();
    t_graph_.clear();
    internal_triangles_.clear();
  }

  /** Return the number of triangles in the mesh.
   * Complexity: O(1).
   */
  size_type size() {
    return num_triangles();
  }

  ////////////////////
  // NODE OPERATORS //
  ////////////////////

  /** Return the number of nodes in the mesh.
   * Complexity: O(1).
   */
  size_type num_nodes() const {
    return n_graph_.num_nodes();
  }

  /** Determine if this Node belongs to this Mesh
   * @return True if @a n is currently a Node of this Mesh
   *
   * Complexity: O(1).
   */
  bool has_node(const node_type& n) const {
    return n_graph_.has_node(n);
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  node_type node(size_type i) const {
    return n_graph_.node(i);
  }

   /** Add a node to the mesh, returning the added node.
   * @param[in] position The new node's position
   * @post new size() == old size() + 1
   * @post result_node.index() == old size()
   *
   * Complexity: O(1) 
   */
  node_type add_node(const Point& position, const node_value_type& value = node_value_type()) {
    internal_node new_node = internal_node(this);
    new_node.value_= value;
    return n_graph_.add_node(position, new_node);
  }

  /** Returns a node_iterator pointing to the first node
   * of the mesh.
   *
   * If the mesh is empty, explicitly returns node_end();
   * Complexity: O(1).
   */
  node_iterator node_begin() const {
    return n_graph_.node_begin();
  }

  /** Returns a node_iterator with a valid mesh
   * and invalid position num_nodes() 
   * Complexity: O(1).
   */
  node_iterator node_end() const {
    return n_graph_.node_end();
  }

  /** Returns a vector of indicies of Triangles adjacent to the node
   * of the mesh.
   *
   * Complexity: O(d) : d is the degree of the node
   */
  std::vector<size_type> adj_tri(const node_type& n) const {
    return n.value().neighbor_triangles_;
  }

  // TODO: make these two a part of the Node class instead of the Mesh class
  neighbor_iterator triangle_begin(const node_type& n) const {
    // assert_invariants();
    if(n.value().neighbor_triangles_.size())
      return NeighborIterator(this, n.index(), 0);
    else
      return triangle_end();
  }

  /** Returns an incident_iterator pointing to an invalid edge
   * of the adjacency list with index_ = graph_->adj_list_[uid_].size()
   */
  neighbor_iterator triangle_end(const node_type& n) const {
    // assert_invariants();
    return NeighborIterator(this, n.index(), n.neighbor_triangles_.size());
  }

  ////////////////////
  // EDGE OPERATORS //
  ////////////////////

  /** Return the total number of edges in the mesh.
   *
   * Complexity: O(1).
   */
  size_type num_edges() const {
    return n_graph_.num_edges();
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this mesh
   * @return true if, for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: O(d), where d = max(degree(@a a), degree(@a b))
   */
  bool has_edge(const node_type& a, const node_type& b) const {
    return n_graph_.has_edge(a,b);
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: O(num_edges())
   * Depreciated, use an iterator instead.
   */
  edge_type edge(size_type i) const {
    return n_graph_.edge(i);
  }

  /** Returns an edge_iterator pointing to the first edge
   * of the adjacency list
   *
   * If the mesh is empty, explicitly returns edge_end();
   * Complexity: O(1). 
   */
  edge_iterator edge_begin() const {
    return n_graph_.edge_begin();
  }

  /** Returns an edge_iterator pointing to the end
   * of the adjacency list with uid1_ of num_edges_
  * Complexity: O(1).
   */
  edge_iterator edge_end() const {
    return n_graph_.edge_end();
  }

  ////////////////
  // TRIANGLES  //
  ////////////////

  /** @class Mesh::Triangle
   * @brief Class representing the mesh's triangles.
   *
   * Triangles are order-insensitive triples of nodes where any pair of
   * nodes are a valid edge in the mesh.  
   * Two Triangles with the same nodes are considered equal if they 
   * connect the same node in either order.
   */
  class Triangle {
   public:
    /** Construct an invalid Triangle. */
    Triangle() {
      tr_id_ = UINT_MAX;
    }

    /** Return this triangle's index, a number in the range [0, mesh_size).
     *
     * @pre This Triangle is valid and is an element of a Mesh 
     */
    const size_type& index() const {
      return tr_id_;
    }

    /** Returns the value of the triangle.
     *
     * @pre The triangle is valid.
     */
    triangle_value_type& value() {
      return mesh_->internal_triangles_.at(index()).value_;
    }


    /** Returns the value of the triangle.
     *
     * @pre The triangle is valid.
     * @post The value returned is read-only.
     */
    const triangle_value_type& value() const {
      return mesh_->internal_triangles_.at(index()).value_;
    }

    /** Returns a node in the triangle
     * 
     * Complexity: O(1)
     */
    node_type node(size_type i) const {
      return mesh_->n_graph_.node(mesh_->internal_triangles_.at(index()).primal_nodes_[i]);
    }

    /** Returns an edge in the triangle.
     * edge(n) returns the edge oposite node(n)
     * Complexity: O(1)
     */
    edge_type edge(const size_type i) const {
      return mesh_->internal_triangles_.at(index()).primal_edges_[i];
    }

    /** Returns an outward normal vector corresponding to edge(i) proportional to edge(i)
     * Complexity: O(1)
     */    
    point_type normal(size_type i) const {
      Point a = node((i+1)%3).position();
      Point b = node((i+2)%3).position();

      // Get the surface normal and normalize
      Point normal = surface_normal(i);
      // Calculate vector BA and get its cross product with previous normal
      Point ba = b - a;
      return cross(ba, normal); 
    }

    /** Returns an outward surface unit normal of the triangle
     * Complexity: O(1)
     */    
    point_type surface_normal(size_type i = 0) const {
      // Point C = node(i), Points A and B are the other nodes in the triangle
      Point a = node((i+1)%3).position();
      Point b = node((i+2)%3).position();
      Point c = node(i).position();
      // Calculate the vectors AC and BC
      Point ac = a - c;
      Point bc = b - c;
      // Get the cross product of AC and BC
      Point normal = cross(ac, bc);
      return normal / norm(normal);
    }

    /** Test whether this triangle and @a x are equal.
     * Equal triangles are from the same mesh and have the same nodes.
     *
     * If either triangle is invalid, false is returned.
     * Complexity: O(1)
     */
    bool operator==(const Triangle& x) const {
      if(mesh_ == NULL || x.mesh_ == NULL) return false;
      return (mesh_ == x.mesh_ && tr_id_ == x.tr_id_);
    }

    /** Returns true if this triangle < @a x
     * Complexity: O(1)
     */
    bool operator<(const Triangle& x) const {
      if(mesh_ != NULL && mesh_ == x.mesh_)
        return (tr_id_ < x.tr_id_);
      else
        return (mesh_ < x.mesh_);
    }

     /** Returns the area of the triangle.
     * Utilizes Heron's formula
     *
     * Complexity: O(1)
     */
    double area() const {
      // Get the edge lengths
      double l0 = edge(0).length();
      double l1 = edge(1).length();
      double l2 = edge(2).length();
      // Calculate s = perimeter/2
      double s = (l0 + l1 + l2)/2.;
      // Heron's formula
      return sqrt(s * (s - l0) * (s - l1) * (s - l2));
    }
    

    /** Returns the neighbor triangle corresponding to edge(i).
     *
     * Complexity: O(1)
     */
    Triangle neighbor(int i) const {
      size_type index = edge(i).value().tri1_;
      if(index == tr_id_) {
        index = edge(i).value().tri2_;
        if(index != UINT_MAX) 
          return Triangle(mesh_, index);
        else 
          return Triangle();
      } else 
        return Triangle(mesh_, index);
    }
    
   private:
    // Allow Mesh to access Triangle's private member data and functions.
    friend class Mesh;

    // Pointer back to Mesh container
    // RI: mesh_ is valid.
    mesh_type* mesh_;

    // This triangle's index
    size_type tr_id_;

    void assert_invariants() const {
      assert(mesh_ != NULL);
    }

    // Private Constructor
    Triangle(const mesh_type* mesh, size_type tr_id)
        : mesh_(const_cast<mesh_type*>(mesh)), tr_id_(tr_id) {}
  };

  class TriangleIterator : private totally_ordered<TriangleIterator> {
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

    /** Construct an invalid TriangleIterator. */
    TriangleIterator() {}

    /** Returns the triangle at the current position
     * @pre The iterator must be valid
     */
    value_type operator*() const {
      assert(mesh_ != NULL);
      return mesh_->triangle((*it_).value().index_); 
    }

    /** Sets this iterator to the next position
     * and returns this iterator
     *
     * @pre The iterator must be valid
     * @post index_ is a valid triangle or mesh_->size() if we're at the end.
     */
    TriangleIterator& operator++() {
      assert(mesh_ != NULL);
      ++it_;
      return *this;
    }

    TriangleIterator& operator[](int n){
      it_ = it_[n];
      return *this;
    }

    TriangleIterator& operator+=(int n){
      it_ += n;
      return *this;
    }
    
    TriangleIterator& operator-=(int n){
      it_ -= n;
      return *this;
    }

    //TriangleIterator& operator-(TriangleIterator& a) {
      //it_ = it_ - a.it_;
    //  return *this;
    //}

    int operator-(TriangleIterator& a) {
      return it_ - a.it_;
    }

    /** Checks to see if two iterators are equal
     *
     * Returns true if both iterators have a valid mesh and
     * mesh_ and index_ are equal.
     * Otherwise, returns false.
     */
    bool operator==(const TriangleIterator& x) const {
      return (mesh_ != NULL && mesh_ == x.mesh_ && it_ == x.it_);
    }

   private:
    friend class Mesh;

    // Pointer back to Mesh container
    mesh_type* mesh_; 

    // The iterator's position in the graph's nodes
    t_node_iterator it_;   

    TriangleIterator(const mesh_type* mesh, t_node_iterator it)
        : mesh_(const_cast<mesh_type*>(mesh)), it_(it) {
    }

  };

  friend TriangleIterator operator+(TriangleIterator a, int n){
    return a+=n;
  }
  
  friend TriangleIterator operator+(int n, TriangleIterator a){
    return a+=n;
  }
  
  friend TriangleIterator operator-(TriangleIterator a, int n){
    return a-=n;
  }

  /** @class Graph::NeighborIterator
   * @brief Iterator class for triangles incident to a node. A forward iterator. */
  class NeighborIterator : private totally_ordered<NeighborIterator> {
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

    /** Construct an invalid NeighborIterator. */
    NeighborIterator() {}

    /** Returns the incident edge at the current position
     * @pre The iterator must be valid
     */
     Triangle operator*() const {
      assert(mesh_ != NULL);
      // assert(mesh_->uid2index_[home_uid_] < mesh_->num_nodes());
      // assert(away_index_ < mesh_->adj_list_[mesh_->uid2index_[home_uid_]].size());

      return mesh_->triangle(mesh_->node(n_index_).neighbor_triangles_[tri_index_]);
    }

    /** Sets this iterator to the next position
     * and returns this iterator
     *
     * @pre The iterator must be valid
     * @post mesh->node(mesh_->adj_list_[uid1_][index_]) is a valid node 
     * or index_ = mesh_->adj_list_[uid1_].size() if we're at the end.
     */
    NeighborIterator& operator++() {
      assert(mesh_ != NULL);
      // assert(mesh_->uid2index_[home_uid_] < mesh_->num_nodes());
      // assert(away_index_ < mesh_->adj_list_[mesh_->uid2index_[home_uid_]].size());

      ++tri_index_;
      return *this;
    }

    /** Checks to see if two iterators are equal
     *
     * Returns true if mesh_, uid1_, and index_ are equal.
     * Otherwise, returns false.
     */
    bool operator==(const NeighborIterator& x) const {
      if(mesh_ == NULL || x.mesh_ == NULL) return false;
      return (mesh_ == x.mesh_ && n_index_ == x.n_index_ && tri_index_ == x.tri_index_);
    }

   private:
    friend class Mesh;

    // Pointer back to Mesh container
    mesh_type* mesh_; 

    // The index of the source node.
    size_type n_index_;

    // Index of the triangle in the neighbor array
    size_type tri_index_;

    NeighborIterator(const mesh_type* mesh, size_type n_index, size_type tri_index)
        : mesh_(const_cast<mesh_type*>(mesh)), n_index_(n_index), tri_index_(tri_index) {
    }
  };

  /** Return the triangle with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: O(1)
   */
  triangle_type triangle(size_type i) const {
    assert(i < t_graph_.size());
    return Triangle(this, t_graph_.node(i).value().index_);
  }

  /** Adds a triangle to the mesh if it doesn't exist.
   * Return the new triangle or the existing triangle.
   * Can invalid edge and triangle indicies.
   * Complexity: O(d), where d is the max(degree(@a n1), degree(@a n2), degree(@a n3))
   */
  triangle_type add_triangle(const node_type& n1, const node_type& n2, const node_type& n3) {
    
    // Assert that nodes exist
    assert(has_node(n1) && has_node(n2) && has_node(n3));
    size_type index = triangle(n1,n2,n3);
    // Return the triangle if it already exists in the mesh.
    if(index != UINT_MAX)
      return Triangle(this, index);

    // Check to see whether it is possible to add the triangle
    // Return an invalid triangle otherwise.
    // (i.e. each edge has a free slot for a triangle)
    if(has_edge(n1,n2)) {
      edge_type e = n_graph_.add_edge(n1, n2);
      if(e.value().tri1_ != UINT_MAX && e.value().tri2_ != UINT_MAX)
        return Triangle();
    } 
    if(has_edge(n2,n3)) {
      edge_type e = n_graph_.add_edge(n2, n3);
      if(e.value().tri1_ != UINT_MAX && e.value().tri2_ != UINT_MAX)
        return Triangle();
    } 
    if(has_edge(n1,n3)) {
      edge_type e = n_graph_.add_edge(n1, n3);
      if(e.value().tri1_ != UINT_MAX && e.value().tri2_ != UINT_MAX)
        return Triangle();
    } 
    // Add the triangle node to the dual graph
    t_node_type new_tri = t_graph_.add_node(Point(0,0,0));

    // Set the primal nodes and edges of added triangle
    new_tri.value().primal_nodes_[0] = n1.index();
    new_tri.value().primal_nodes_[1] = n2.index();
    new_tri.value().primal_nodes_[2] = n3.index();
    new_tri.value().primal_edges_[0] = n_graph_.add_edge(n3, n2);
    new_tri.value().primal_edges_[1] = n_graph_.add_edge(n3, n1);
    new_tri.value().primal_edges_[2] = n_graph_.add_edge(n2, n1);

    // Test whether triangle already exists by comparing e0.tri1 with e1.tri1
    // Also return invalid triangle if making one is Impossibru
    for (size_type i = 0; i < 3; ++i) {
      // Set added triangle as a neighbor to each primal node
      n_graph_.node(new_tri.value().primal_nodes_[i]).value().neighbor_triangles_.push_back(new_tri.index());
      // Set added triangle as the first or second neighbor to each primal edge
      if (new_tri.value().primal_edges_[i].value().tri1_ != UINT_MAX)
        new_tri.value().primal_edges_[i].value().tri2_ = new_tri.index();
      else
        new_tri.value().primal_edges_[i].value().tri1_ = new_tri.index();
    }
    // Set the value and index of the node and push onto internal trianglse vector
    new_tri.value().value_ = triangle_value_type();
    new_tri.value().index_ = new_tri.index();
    internal_triangles_[new_tri.index()] = new_tri.value();
    return Triangle(this, new_tri.index());

  }

  /** Grabs the index of the triangle with the three nodes
   *  @return index of triangle defined by input nodes if triangle exists, else UINT_MAX
   * Complexity: O(d), where d is the max(degree(@a n1), degree(@a n2), degree(@a n3))
   */
  size_type triangle(const node_type& n1, const node_type& n2, const node_type& n3) {
    if(!has_node(n1) || !has_node(n2) || !has_node(n3)){
      return UINT_MAX;
    }
    // Using edges to check for a common triangle instead of using 
    // the adjacent triangle vectors of each node is less pretty but has
    // a faster runtime. (O(d) vs O(d log d))
    size_type index1;
    size_type index2;
    if(has_edge(n1, n2)) {
      edge_type edge = n_graph_.add_edge(n1, n2);
      index1 = edge.value().tri1_;
      index2 = edge.value().tri2_;
      if(has_edge(n1, n3)) {
        size_type final_index;
        edge = n_graph_.add_edge(n1, n3);
        if(index1 == edge.value().tri1_ || index1 == edge.value().tri2_) {
          final_index = index1;
        } else if(index2 == edge.value().tri1_ || index2 == edge.value().tri2_) {
          final_index = index2;         
        } else
          return UINT_MAX;

        if(has_edge(n2, n3)) {
          edge = n_graph_.add_edge(n2,n3);
          if(final_index == edge.value().tri1_ || final_index == edge.value().tri2_) {
            return final_index;
          } else
            return UINT_MAX;
        } else
          return UINT_MAX;
      } else 
        return UINT_MAX;
    } else 
      return UINT_MAX;
  }

  /** Test whether a triangle is in the mesh.
   *
   * Complexity: O(1)
   */
  bool has_triangle(const triangle_type& t) {
    return t_graph_.has_node(t.index_);
  }

  /** Return the number of triangles in the mesh.
   * Complexity: O(1).
   */
  size_type num_triangles() const {
    return t_graph_.num_nodes();
  }

  /** Returns an triangle_iterator pointing to the first triangle
   * of the Mesh
   *
   * Complexity: O(1).
   */
  triangle_iterator triangle_begin() const {
    return triangle_iterator(this, t_graph_.node_begin());
  };

  /** Returns an invalid triangle_iterator.
   * Complexity: O(1).
   */
  triangle_iterator triangle_end() const {
    return triangle_iterator(this, t_graph_.node_end());
  }

 private:

  // Define the primal and dual graphs
  NodeGraph n_graph_;
  TriangleGraph t_graph_;

  // Struct for proxy representation of triangles.
  struct internal_triangle {
    std::array<size_type, 3> primal_nodes_;
    std::array<edge_type, 3> primal_edges_;
    std::array<size_type, 3> neighbor_triangles_;
    triangle_value_type value_;
    size_type index_;
  };
  
  // Struct for proxy representation of internal nodes
  struct internal_node {

    // Default Constructor
    internal_node(const mesh_type* mesh)
      : mesh_(const_cast<mesh_type*>(mesh)) {
      value_ = node_value_type();
      neighbor_triangles_ = std::vector<size_type>();
    }

    triangle_type neighbor_triangle(size_type i) {
      assert(i < neighbor_triangles_.size());
      return mesh_->triangle(neighbor_triangles_[i]);
    }
    mesh_type* mesh_;
    std::vector<size_type> neighbor_triangles_;
    node_value_type value_;
  }; 
  
  // Struct for proxy representation of internal edges
  struct internal_edge {
    // Default Constructor
    internal_edge() {
      tri1_ = UINT_MAX;
      tri2_ = UINT_MAX;
      value_ = edge_value_type();
    }
    size_type tri1_;
    size_type tri2_;
    edge_value_type value_;
  };

  /**
   *  AF(T) = T = {n_1, n_2, ...n_(num_triangles()-1)}
   *  RI:  internal_triangles_.size() == num_triangles();
   */
  std::unordered_map<size_type, internal_triangle> internal_triangles_;
};
