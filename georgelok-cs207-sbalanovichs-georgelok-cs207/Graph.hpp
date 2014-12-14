#ifndef CS207_GRAPH_HPP
#define CS207_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <vector>
#include <cassert>

#include "CS207/Util.hpp"
#include "Point.hpp"
#include <limits>

/** Alias for max int size. **/
unsigned int MAX_UNSIGNED_INT_SIZE = std::numeric_limits<unsigned int>::max();

/** @class Graph
 * @brief A template for 3D undirected graphs.
 *
 * Users can add and retrieve nodes and edges. Edges are unique (there is at
 * most one edge between any pair of distinct nodes).
 */
template <typename V = int, typename E = double>
class Graph {
 private:
  /** Type of points in this graph **/
  typedef Point point_type;

  /** Type of node_pair_uid **/
  typedef unsigned long long node_pair_uid_type;

  /** Predeclare intenal structs **/
  struct internal_node;
  struct internal_edge;

 public:

  /////////////////////////////
  // PUBLIC TYPE DEFINITIONS //
  /////////////////////////////

  /** Synonym for V **/
  typedef V node_value_type;

  /** Synonym for E **/
  typedef E edge_value_type;

  /** Type of this graph. */
  typedef Graph graph_type;

  /** Predeclaration of Node type. */
  class Node;
  /** Synonym for Node (following STL conventions). */
  typedef Node node_type;

  /** Predeclaration of Edge type. */
  class Edge;
  /** Synonym for Edge (following STL conventions). */
  typedef Edge edge_type;

  /** Type of indexes and sizes.
      Return type of Graph::Node::index(), Graph::num_nodes(),
      Graph::num_edges(), and argument type of Graph::node(size_type) */
  typedef unsigned size_type;

  /** Type of node iterators, which iterate over all graph nodes. */
  class NodeIterator;
  /** Synonym for NodeIterator */
  typedef NodeIterator node_iterator;

  /** Type of edge iterators, which iterate over all graph edges. */
  class EdgeIterator;
  /** Synonym for EdgeIterator */
  typedef EdgeIterator edge_iterator;

  /** Type of incident iterators, which iterate incident edges to a node. */
  class IncidentIterator;
  /** Synonym for IncidentIterator */
  typedef IncidentIterator incident_iterator;

  ////////////////////////////////
  // CONSTRUCTOR AND DESTRUCTOR //
  ////////////////////////////////

  /** Construct an empty graph. */
  Graph() : 
      num_edges_(0), 
      first_edge_index_(0), 
      internal_nodes_({}), 
      uid2index_({}), 
      empty_node_uids_({}), 
      adj_list_({}), 
      internal_edges_({}), 
      empty_edge_indicies_({}) {}

  /** Default destructor */
  ~Graph() = default;

  /////////////
  // General //
  /////////////

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    return internal_nodes_.size();
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   * Complexity: O(1).
   */
  void clear() {
    num_edges_ = 0;
    first_edge_index_ = 0;
    internal_nodes_.clear();
    uid2index_.clear();
    empty_node_uids_.clear();
    adj_list_.clear();
    internal_edges_.clear();
    empty_edge_indicies_.clear();
  }

  /////////////////
  // GRAPH NODES //
  /////////////////

  /** @class Graph::Node
   * @brief Class representing the graph's nodes.
   *
   * Node objects are used to access information about the Graph's nodes.
   */
  class Node : private totally_ordered<Node> {
   public:
    /** Construct an invalid node.
     *
     * Valid nodes are obtained from the Graph class, but it
     * is occasionally useful to declare an @i invalid node, and assign a
     * valid node to it later. For example:
     *
     * @code
     * Graph::node_type x;
     * if (...should pick the first node...)
     *   x = graph.node(0);
     * else
     *   x = some other node using a complicated calculation
     * do_something(x);
     * @endcode
     */
    Node() {}

    /** Return this node's position.
     *
     * @pre This Node is valid and is an element of a Graph 
     * @post The value returned is read only
     */    
    const Point& position() const {
      assert_invariants();
      return graph_->internal_nodes_[index()].position_;
    }

    /** Return this node's position.
     *
     * @pre This Node is valid and is an element of a Graph 
     */
    Point& position() {
      assert_invariants();
      return graph_->internal_nodes_[index()].position_;
    }

    /** Return this node's index, a number in the range [0, graph_size).
     *
     * @pre This Node is valid and is an element of a Graph 
     */
    size_type index() const {
      assert_invariants();
      return graph_->uid2index_[uid_];
    }

    /** Valid nodes are equal if they have the same graph and same uid.
     *
     * If either node is invalid, false is returned.
     */
    bool operator==(const Node& x) const {
      if(graph_ == NULL || x.graph_ == NULL) return false;
      return (graph_ == x.graph_ && uid_ == x.uid_);
    }

    /** Returns the node with first the smallest graph and then
     * the smallest uid
     *
     * @pre this is a valid node and x is a valid node. 
     */
    bool operator<(const Node& x) const {
      if(graph_ != NULL && graph_ == x.graph_)
        return (uid_ < x.uid_);
      else
        return (graph_ < x.graph_);
    }

    /** Returns the value of the node.
     *
     * @pre The node is valid.
     */
    node_value_type& value() {
      assert_invariants();
      return graph_->internal_nodes_[index()].value_;
    }

    /** Returns the value of the node.
     *
     * @pre The node is valid.
     * @post The value returned is read-only.
     */
    const node_value_type& value() const {
      assert_invariants();
      return graph_->internal_nodes_[index()].value_;
    }

    /** Returns the degree of the node.
     *
     * @pre The node is valid.
     */
    size_type degree() const {
      assert_invariants();
      return graph_->adj_list_[index()].size();
    }

    /** Returns an incident_iterator pointing to the first node
     * of the adjacency list of a given home_node
     *
     * If the list is empty, returns edge_end()
     */
    incident_iterator edge_begin() const {
      assert_invariants();
      if(graph_->adj_list_[index()].size() )
        return IncidentIterator(graph_, uid_, 0);
      else
        return edge_end();
    }

    /** Returns an incident_iterator pointing to an invalid edge
     * of the adjacency list with index_ = graph_->adj_list_[uid_].size()
     */
    incident_iterator edge_end() const {
      assert_invariants();
      return IncidentIterator(graph_, uid_, graph_->adj_list_[index()].size());
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;

    // Pointer back to Graph container
    // RI: graph_ is valid.
    graph_type* graph_;

    // This node's UID
    // RI: uid2index_[uid_] < graph_->num_nodes()
    size_type uid_;

    void assert_invariants() const {
      assert(graph_ != NULL);
      assert(uid_ < graph_->uid2index_.size());
      assert(graph_->uid2index_[uid_] < graph_->num_nodes());   
    }

    // Private Constructor
    Node(const graph_type* graph, size_type uid)
        : graph_(const_cast<graph_type*>(graph)), uid_(uid) {}
  };

  /** Synonym for size(). */
  size_type num_nodes() const {
    return size();
  }

  /** Add a node to the graph, returning the added node.
   * @param[in] position The new node's position
   * @post new size() == old size() + 1
   * @post result_node.index() == old size()
   *
   * Complexity: O(1) 
   */
  Node add_node(const Point& position, const node_value_type& value = node_value_type()) {
  	// Get uid for node
  	size_type uid;

  	if (!empty_node_uids_.empty()) {
  		uid = empty_node_uids_.back();
  		empty_node_uids_.pop_back();
  		uid2index_[uid] = num_nodes();
  	}
  	else {
  		uid = uid2index_.size();
  		uid2index_.push_back(num_nodes());
  	}

    // Create and add proxy node
    internal_node new_node = {uid, uid2index_[uid], position, value};
    internal_nodes_.push_back(new_node);

    // Add entry in adjacency list
    std::vector<size_type> new_list;
    adj_list_.push_back(new_list);

    // Update first_edge_home_node if needed
    if(num_edges_ == 0) first_edge_index_ = num_nodes();

    // Return our node.
    return Node(this, uid);
  }

  /** Remove a node from the graph.
   * @param[in] n Node to be removed
   * @return 1 if old has_node(n), 0 otherwise
   *
   * @post new size() == old size() - result.
   *
   * Can invalidate outstanding iterators. 
   * If old has_node(@a n), then @a n becomes invalid, as do any
   * other Node objects equal to @a n. All other Node objects remain valid.
   *
   * Also, for all nodes n2 such that has_edge(n, n2) is true, all such edges
   * are also invalidated.
   *
   * Complexity: The larger of O(N) and 
   *              O(n.degree() *d) : d is the maximum degree of the graph.
   */
  size_type remove_node(const Node& n) {
    if(!has_node(n))
      return 0;

    // Get index of original node
    size_t index = n.index();

    // Update number of edges
    num_edges_ -= adj_list_[index].size();

    // Remove all edges associated with this node.
    for(auto it = n.edge_begin(); it != n.edge_end(); ++it) {
      size_type node2_index = (*it).node2().index();
      for(auto adj_it = adj_list_[node2_index].begin(); adj_it != adj_list_[node2_index].end(); ++adj_it) {
        if(internal_edges_[(*adj_it)].uid1_ == n.uid_ || internal_edges_[(*adj_it)].uid2_ == n.uid_) {
          empty_edge_indicies_.push_back((*adj_it));
          internal_edges_[(*adj_it)].uid1_ = 0;
          internal_edges_[(*adj_it)].uid2_ = 0; 
          adj_list_[node2_index].erase(adj_it);
          break;
        }
      }
    }
    adj_list_.erase(adj_list_.begin() + index);

    // Erase item from nodes vector.
    auto it = internal_nodes_.erase(internal_nodes_.begin() + index);

    // Update all metadata
    for(; it != internal_nodes_.end(); ++it) {
      (*it).index_ -= 1;
      uid2index_[(*it).uid_] -= 1;
    }

    // Add uid to empty nodes
    empty_node_uids_.push_back(n.uid_);
    return 1;
  }

  /** Remove the node pointed to by n_it from the graph
   * @param[in] n_it iterator to Node to be removed
   * @return node_iterator equivalent to old ++n_it
   * The returned node iterator is guaranteed to be either a valid
   * iterator or new node_end()
   *
   * @post new size() == old size() - -1
   *
   * Invalidates old node_end() and the value pointed to by other iterators
   * may change. 
   *
   * If old has_node(@a *n_it), then @a *n_it becomes invalid, as do any
   * other Node objects equal to @a *n_it. All other valid Node objects remain valid.
   * 
   * Also, for all nodes n2 such that has_edge(*n_it, n2) is true, all such edges
   * are also invalidated.
   *
   * Complexity: The larger of O(N) and 
   *              O(n.degree() *d) : d is the maximum degree of the graph.
   */ 
  node_iterator remove_node(node_iterator n_it) {
    remove_node(*n_it);
    return n_it;
  }  

  /** Determine if this Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    return (n.index() < num_nodes()) && (node(n.index()) == n);
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    assert(i < num_nodes());
    return Node(this, internal_nodes_[i].uid_);
  }

  /////////////////
  // GRAPH EDGES //
  /////////////////

  /** @class Graph::Edge
   * @brief Class representing the graph's edges.
   *
   * Edges are order-insensitive pairs of nodes. Two Edges with the same nodes
   * are considered equal if they connect the same nodes, in either order.
   */
  class Edge : private totally_ordered<Edge> {
   public:
    /** Construct an invalid Edge. */
    Edge() {}

    /** If this Edge is valid, return a node of this Edge 
     */
    Node node1() const {
      assert(graph_ != NULL);
      assert(graph_->uid2index_[node1_uid_] < graph_->num_nodes());
      return graph_->node(graph_->uid2index_[node1_uid_]);
    }

    /** If this Edge is valid, return the other node of this Edge 
     *
     * Return Node() if the Edge is invalid
     */
    Node node2() const {
      assert(graph_ != NULL);
      assert(graph_->uid2index_[node2_uid_] < graph_->num_nodes());
      return graph_->node(graph_->uid2index_[node2_uid_]);
    }

    /** Test whether this edge and @a x are equal.
     *
     * Equal edges are from the same graph and have the same nodes.
     *
     * If either node is invalid, false is returned.
     */
    bool operator==(const Edge& x) const {
      if(graph_ == NULL || x.graph_ == NULL) return false;
      if(graph_ != x.graph_) return false;
      return (node_pair_uid(node1_uid_, node2_uid_) == 
              node_pair_uid(x.node1_uid_, x.node2_uid_));
    }

    /** Returns the edge with first the smallest graph and then
     * the smallest node_pair_uid
     *
     * @pre this is a valid edge and x is a valid edge. 
     */
    bool operator<(const Edge& x) const {
      if(graph_ != NULL && graph_ == x.graph_)
        return (node_pair_uid(node1_uid_, node2_uid_) < node_pair_uid(x.node1_uid_, x.node2_uid_));
      else
        return graph_< x.graph_;
    }

     /** Returns the length of the edge.
     *
     * @pre this is a valid edge.
     * @post the value is read only.
     */
    double length() const {
      assert(graph_ != NULL);
      return norm(node2().position() - node1().position());
    }


    /** Returns the value of the edge.
     *
     * @pre this is a valid edge.
     * @post the value is read only.
     */
    const edge_value_type& value() const {
    	assert(graph_ != NULL);
      return graph_->internal_edges_[index_].value_;  
    }

    /** Returns the value of the edge.
     *
     * @pre this is a valid edge.
     */
    edge_value_type& value() {
      assert(graph_ != NULL);
      return graph_->internal_edges_[index_].value_;  
    }

    /** Returns a unique index of the edge.
     *
     * @pre this is a valid edge.
     * @post edge a == edge b => a.index() = b.index()
     * @post edge a != edge b => a.index() != b.index()
     */
    size_type index() const {
      assert(graph_ != NULL);
      return index_;  
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;

    // Pointer back to Graph container
    graph_type* graph_;

    // UID for each individual node.
    // RI: uid2index_[node1_uid_], uid2index_[node2_uid_] < graph_->num_nodes()
    size_type node1_uid_;
    size_type node2_uid_;

    // Index where edge data is stored
    // RI: internal_edges_[index_] is valid
    size_type index_;

    /** Returns a node_pair_uid given the uid of two nodes.
     * @pre a,b refer to valid nodes in the same graph.
     * @pre a != b
     * @post node_pair_uid(a,b) is unique for unordered pair (a,b)
     *
     * Complexity: O(1).
     */
    static node_pair_uid_type node_pair_uid(size_type a, size_type b) {
      node_pair_uid_type left, right;
      if (a < b) {
        left = (node_pair_uid_type) a;
        right = (node_pair_uid_type) b;
      } else {
        left = (node_pair_uid_type) b;
        right = (node_pair_uid_type) a;
      }
      return (left << 32) + right; 
    }

    // Private Constructor
    Edge(const graph_type* graph, size_type node1_uid, size_type node2_uid, size_type index)
        : graph_(const_cast<graph_type*>(graph)), node1_uid_(node1_uid), node2_uid_(node2_uid), index_(index)  {
    }
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: O(1).
   */
  size_type num_edges() const {
    return num_edges_;
  }

  /** Add an edge to the graph, or return the current edge if it already exists.
   * @pre @a a and @a b are distinct valid nodes of this graph
   * @return an Edge object e with e.node1() == @a a and e.node2() == @a b
   * @post has_edge(@a a, @a b) == true
   * @post If old has_edge(@a a, @a b), new num_edges() == old num_edges().
   *       Else,                        new num_edges() == old num_edges() + 1.
   *
   * Can invalidate edge indexes -- in other words, old edge(@a i) might not
   * equal new edge(@a i).
   *
   * Complexity: O(d), where d is larger degree of the two nodes.
   */
  Edge add_edge(const Node& a, const Node& b, const edge_value_type& value = edge_value_type()) {
    // Sanity Checks
    assert(a.uid_ != b.uid_);
  	assert(a.uid_ < uid2index_.size() && b.uid_ < uid2index_.size());
  	assert(uid2index_[a.uid_] < num_nodes() && uid2index_[b.uid_] < num_nodes());
    
    // Check if the edge already exists
    size_type home_node, away_node;
    if (uid2index_[a.uid_] > uid2index_[b.uid_]) {
      away_node = a.uid_; home_node = b.uid_;
    }
    else {
      away_node = b.uid_; home_node = a.uid_;
    }

    size_type list_length = adj_list_[uid2index_[home_node]].size();
    for(size_type i = 0; i < list_length; ++i) {
      if(internal_edges_[adj_list_[uid2index_[home_node]][i]].uid2_ == away_node)
      	return Edge(this, a.uid_, b.uid_, adj_list_[uid2index_[home_node]][i]);
    }

		// Create new internal_edge
    internal_edge new_edge;

    // Grab index of edge 
    size_t index;
    if (!empty_edge_indicies_.empty()) {
  		index = empty_edge_indicies_.back();
  		empty_edge_indicies_.pop_back();
  	}
  	else {
  		index = internal_edges_.size();
  	  internal_edges_.push_back(new_edge);
  	}

  	// set values of internal edge.
    new_edge = {index, home_node, away_node, value};
    internal_edges_[index] = new_edge;

    // Update adjacency list
    adj_list_[uid2index_[a.uid_]].push_back(index);
    adj_list_[uid2index_[b.uid_]].push_back(index);

    // Update the first_edge_index_
    if(first_edge_index_ > home_node) 
      first_edge_index_ = uid2index_[home_node];

    ++num_edges_;

    return Edge(this, a.uid_, b.uid_, index);
  }

  /** Remove a edge from the graph.
   * @param[in] n1 node 1 of edge to be removed
   * @param[in] n2 node 2 of edge to be removed
   * @return 1 if old has_edge(n1, n2), 0 otherwise
   *
   * @post new size() == old size() - result.
   *
   * Can invalidate outstanding iterators. 
   *
   * If old has_edge(@a n1, @n2), then the edge and anything equivalent
   * become invalid.  All other edges remain valid.
   *
   * Complexity: O(d) : d is the maximum degree of the graph.
   */
  size_type remove_edge(const Node& n1, const Node& n2) {
    if(!has_edge(n1, n2)) return 0;

    num_edges_--;

    size_type index = uid2index_[n1.uid_];
    for(auto it = adj_list_[index].begin(); it != adj_list_[index].end(); ++it) {
      if(internal_edges_[(*it)].uid1_ == n2.uid_ || internal_edges_[(*it)].uid2_ == n2.uid_) {
          adj_list_[index].erase(it);
          break;
        }
    }

    index = uid2index_[n2.uid_];
    for(auto it = adj_list_[index].begin(); it != adj_list_[index].end(); ++it) {
      if(internal_edges_[(*it)].uid1_ == n1.uid_ || internal_edges_[(*it)].uid2_ == n1.uid_) {
          empty_edge_indicies_.push_back((*it));
          internal_edges_[(*it)].uid1_ = 0;
          internal_edges_[(*it)].uid2_ = 0; 
          adj_list_[index].erase(it);
          break;
        }
    }

    return 1;
  }

  /** Remove a edge from the graph.
   * @param[in] e edge to be removed
   * @return 1 if old has_edge(e.node1(), e.node2()), 0 otherwise
   *
   * @post new size() == old size() - result.
   *
   * Can invalidate outstanding iterators. 
   *
   * If old has_edge(@a e.node1(), @a e.node2()), then the edge and anything equivalent
   * become invalid.  All other edges remain valid.   *
   * Complexity: O(d) : d is the maximum degree of the graph.
   */
  size_type remove_edge(const Edge& e) {
    return remove_edge(e.node1(),  e.node2());
  }

  /** Remove a edge from the graph.
   * @param[in] it iterator pointing to edge to be removed
   * @return a iterator pointing to the edge of old ++it
   *          The returned iterator is guaranteed to be
   *          either new edge_end() or a valid iterator.
   *
   * @post new size() == old size() - 1. if it is valid
   *
   * Can invalidate outstanding iterators. 
   * If old has_edge(@a *e_it).node1(), (*e_it).node2()), 
   * then the edge becomes invalid, as do any other edges equal to it.
   * All other valid edges remain valid.
   *
   * Complexity: O(d) : d is the maximum degree of the graph.
   */
  edge_iterator remove_edge(edge_iterator e_it) {
    assert(e_it.graph_ == this);
    remove_edge((*e_it).node1(), (*e_it).node2());

    // If e_it is now invalid, find the next valid iterator
    if(e_it.away_index_ >= adj_list_[e_it.home_index_].size() ||
      uid2index_[internal_edges_[adj_list_[e_it.home_index_][e_it.away_index_]].uid2_] 
              <= e_it.home_index_) 
    {
      // Find the next edge in the adjacency list
      for(; e_it.home_index_ < num_nodes(); ++e_it.home_index_) { 
        size_type list_length = adj_list_[e_it.home_index_].size();
        for(; e_it.away_index_ < list_length; ++e_it.away_index_) {
          if (uid2index_[internal_edges_[adj_list_[e_it.home_index_][e_it.away_index_]].uid2_] 
              > e_it.home_index_)
            return e_it;
        }
        e_it.away_index_ = 0;
      }
    }
    return e_it;
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return true if, for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: O(d), where d is the largest degree of the graph
   */
  bool has_edge(const Node& a, const Node& b) const {
    // Sanity Check
    if(a.uid_ == b.uid_) return false;
		if(!(a.uid_ < uid2index_.size() && b.uid_ < uid2index_.size())) return false;
  	if(!(uid2index_[a.uid_] < num_nodes() && uid2index_[b.uid_] < num_nodes())) return false;
    
    size_type home_node, away_node;
    if (uid2index_[a.uid_] > uid2index_[b.uid_]) {
      away_node = a.uid_; home_node = b.uid_;
    }
    else {
      away_node = b.uid_; home_node = a.uid_;
    }

    // Check if we can find the edge
    size_type list_length = adj_list_[uid2index_[home_node]].size();
    for(size_type i = 0; i < list_length; ++i) {
      if(internal_edges_[adj_list_[uid2index_[home_node]][i]].uid2_ == away_node)
      {
        return true;
      }
    }
    return false;
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: O(num_edges())
   * Depreciated, use an iterator instead.
   */
  Edge edge(size_type i) const {
    edge_iterator it = edge_begin();
    for ( ; i != 0; --i) ++it;
    return *it;
  }

  ///////////////
  // Iterators //
  ///////////////

  /** @class Graph::NodeIterator
   * @brief Iterator class for nodes. A forward iterator. 
   */
  class NodeIterator : private totally_ordered<NodeIterator> {
   public:
    // These type definitions help us use STL's iterator_traits.
    /** Element type. */
    typedef Node value_type;
    /** Type of pointers to elements. */
    typedef Node* pointer;
    /** Type of references to elements. */
    typedef Node& reference;
    /** Iterator category. */
    typedef std::input_iterator_tag iterator_category;
    /** Difference between iterators */
    typedef std::ptrdiff_t difference_type;

    /** Construct an invalid NodeIterator. */
    NodeIterator() {}

    /** Returns the node at the current position
     * @pre The iterator must be valid
     */
    Node operator*() const {
      assert(graph_ != NULL);
      assert(index_ < graph_->num_nodes());
      return graph_->node(index_); 
    }

    /** Sets this iterator to the next position
     * and returns this iterator
     *
     * @pre The iterator must be valid
     * @post index_ is a valid internal node or num_nodes() if we're at the end.
     */
    NodeIterator& operator++() {
      assert(graph_ != NULL);
      assert(index_ < graph_->num_nodes());
      ++index_;
      return *this;
    }

    /** Checks to see if two iterators are equal
     *
     * Returns true if both iterators have a valid graph and
     * graph_ and index_ are equal.
     * Otherwise, returns false.
     */
    bool operator==(const NodeIterator& x) const {
      return (graph_ != NULL && graph_ == x.graph_ && index_ == x.index_);
    }

    NodeIterator& operator+=(int n){
      index_ += n;
      return *this;
    }
    
    NodeIterator& operator-=(int n){
      index_ -= n;
      return *this;
    }
    
    NodeIterator& operator[](int n){
      index_ = n;
      return *this;
    }
    
    int operator-(node_iterator a){
      return index_ - a.index_;
    }

    //node_iterator operator-(node_iterator a) {
    //  index_ -= a.index_;
    //  return *this;
    //}

   private:
    friend class Graph;

    // Pointer back to Graph container
    graph_type* graph_; 

    // The iterator's position in the graph's nodes
    // RI: index_ < graph_->num_nodes()
    size_type index_;     

    NodeIterator(const graph_type* graph, size_type index)
        : graph_(const_cast<graph_type*>(graph)), index_(index) {
    }

  };

  friend node_iterator operator+(node_iterator a, int n){
    return a+=n;
  }
  
  friend node_iterator operator+(int n, node_iterator a){
    return a+=n;
  }
  
  friend node_iterator operator-(node_iterator a, int n){
    return a-=n;
  }

  /** Returns a node_iterator pointing to the first node
   * of the graph.
   *
   * If the graph is empty, explicitly returns node_end();
   */
  node_iterator node_begin() const {
    if(num_nodes() > 0)
      return NodeIterator(this, 0);
    else
      return node_end();
  }

  /** Returns a node_iterator with a valid graph
   * and invalid position num_nodes() */
  node_iterator node_end() const {
    return NodeIterator(this, num_nodes());
  }


  /** @class Graph::EdgeIterator
   * @brief Iterator class for edges. A forward iterator. */
  class EdgeIterator : private totally_ordered<EdgeIterator>  {
   public:
    // These type definitions help us use STL's iterator_traits.
    /** Element type. */
    typedef Edge value_type;
    /** Type of pointers to elements. */
    typedef Edge* pointer;
    /** Type of references to elements. */
    typedef Edge& reference;
    /** Iterator category. */
    typedef std::input_iterator_tag iterator_category;
    /** Difference between iterators */
    typedef std::ptrdiff_t difference_type;

    /** Construct an invalid EdgeIterator. */
    EdgeIterator() {}

    /** Returns the edge at the current position
     * @pre The iterator must be valid
     */
    Edge operator*() const {
      assert(graph_ != NULL);
      assert(home_index_ < graph_->num_nodes());
      assert(away_index_ < graph_->adj_list_[home_index_].size());

      return Edge(graph_, 
      			graph_->internal_edges_[graph_->adj_list_[home_index_][away_index_]].uid1_, 
      			graph_->internal_edges_[graph_->adj_list_[home_index_][away_index_]].uid2_,
						graph_->adj_list_[home_index_][away_index_]);
    }

    /** Sets this iterator to the next position
     * and returns this iterator
     *
     * @pre The iterator must be valid
     * @post graph_->adj_list_[home_index_][away_index_] is a valid internal node 
     * or home_index_ = num_nodes() if we're at the end.
     */
    EdgeIterator& operator++() {
      assert(graph_ != NULL);
      assert(home_index_ < graph_->num_nodes());
      assert(away_index_ < graph_->adj_list_[home_index_].size());

      ++away_index_;

      // Find the next edge in the adjacency list
      for(; home_index_ < graph_->num_nodes(); ++home_index_) { 
        size_type list_length = graph_->adj_list_[home_index_].size();
        for(; away_index_ < list_length; ++away_index_) {
          if (graph_->uid2index_[graph_->internal_edges_[graph_->adj_list_[home_index_][away_index_]].uid2_] 
          		> home_index_)
            return *this;
        }
        away_index_ = 0;
      }

      // Explictly set implicit values
      home_index_ = graph_->num_nodes(); away_index_ = 0;
      return *this;
    }

    /** Checks to see if two iterators are equal
     *
     * Returns true if graph_, uid1_, and away_index_ are equal.
     * Otherwise, returns false.
     */
    bool operator==(const EdgeIterator& x) const {
      if(graph_ == NULL || x.graph_ == NULL) return false;
      return (graph_ == x.graph_ && home_index_ == x.home_index_ && away_index_ == x.away_index_);
    }

   private:
    friend class Graph;

    // Pointer back to Graph container
    graph_type* graph_; 

    // The index of the node with the smaller index.
    size_type home_index_;

    /**
     * The index in the list of neighbors of home_index_
     * that corresponds to the second node of the edge.
     * RI:  away_index_ < graph_->adj_list_[home_index_].size() 
     *      || (away_index_ == 0 && home_index_ == graph->num_nodes())
     */
    size_type away_index_;    

    // Private constructor
    EdgeIterator(const graph_type* graph, size_type home_index, size_type away_index)
        : graph_(const_cast<graph_type*>(graph)), home_index_(home_index), away_index_(away_index) {
    }
  };

  /** Returns an edge_iterator pointing to the first edge
   * of the adjacency list
   *
   * If the graph is empty, explicitly returns edge_end();
   */
  edge_iterator edge_begin() const {
    if(num_edges_ > 0)
      return EdgeIterator(this, first_edge_index_, 0);
    else
      return edge_end();
  }

  /** Returns an edge_iterator pointing to the end
   * of the adjacency list with uid1_ of num_edges_
   */
  edge_iterator edge_end() const {
    return EdgeIterator(this, num_nodes(), 0);
  }

  /** @class Graph::IncidentIterator
   * @brief Iterator class for edges incident to a node. A forward iterator. */
  class IncidentIterator : private totally_ordered<IncidentIterator> {
   public:
    // These type definitions help us use STL's iterator_traits.
    /** Element type. */
    typedef Edge value_type;
    /** Type of pointers to elements. */
    typedef Edge* pointer;
    /** Type of references to elements. */
    typedef Edge& reference;
    /** Iterator category. */
    typedef std::input_iterator_tag iterator_category;
    /** Difference between iterators */
    typedef std::ptrdiff_t difference_type;

    /** Construct an invalid IncidentIterator. */
    IncidentIterator() {}

    /** Returns the incident edge at the current position
     * @pre The iterator must be valid
     */
     Edge operator*() const {
      assert(graph_ != NULL);
      assert(graph_->uid2index_[home_uid_] < graph_->num_nodes());
      assert(away_index_ < graph_->adj_list_[graph_->uid2index_[home_uid_]].size());
      // Determine uid of neighbor node
      size_t away_uid_;
      if(graph_->internal_edges_[graph_->adj_list_[graph_->uid2index_[home_uid_]][away_index_]].uid1_ == home_uid_)
      	away_uid_ = graph_->internal_edges_[graph_->adj_list_[graph_->uid2index_[home_uid_]][away_index_]].uid2_;
      else
      	away_uid_ = graph_->internal_edges_[graph_->adj_list_[graph_->uid2index_[home_uid_]][away_index_]].uid1_;

      return Edge(graph_, 
      			home_uid_, 
      			away_uid_,
						graph_->adj_list_[graph_->uid2index_[home_uid_]][away_index_]);
    }

    /** Sets this iterator to the next position
     * and returns this iterator
     *
     * @pre The iterator must be valid
     * @post graph->node(graph_->adj_list_[uid1_][index_]) is a valid node 
     * or index_ = graph_->adj_list_[uid1_].size() if we're at the end.
     */
    IncidentIterator& operator++() {
      assert(graph_ != NULL);
      assert(graph_->uid2index_[home_uid_] < graph_->num_nodes());
      assert(away_index_ < graph_->adj_list_[graph_->uid2index_[home_uid_]].size());

      ++away_index_;
      return *this;
    }

    /** Checks to see if two iterators are equal
     *
     * Returns true if graph_, uid1_, and index_ are equal.
     * Otherwise, returns false.
     */
    bool operator==(const IncidentIterator& x) const {
      if(graph_ == NULL || x.graph_ == NULL) return false;
      return (graph_ == x.graph_ && home_uid_ == x.home_uid_ && away_index_ == x.away_index_);
    }

   private:
    friend class Graph;

    // Pointer back to Graph container
    graph_type* graph_; 

    // The index of the home node.
    size_type home_uid_;

    /**
     * The index in the list of neighbors of home_uid-
     * that corresponds to the second node of the edge.
     * RI:  away_index_ < graph_->adj_list_[home_uid_].size() 
     *      || (away_index_ == 0 && home_uid_ == graph->num_nodes())
     */
    size_type away_index_;  

    IncidentIterator(const graph_type* graph, size_type home_uid, size_type away_index)
        : graph_(const_cast<graph_type*>(graph)), home_uid_(home_uid), away_index_(away_index) {
    }
  };

 private:
  /**
   * Variables to keep track of the number and indicies of nodes and edges.
   * RI:  for all i < num_edges_, edge(i) returns a valid edge
   */      
  size_type num_edges_;

  // Struct for proxy representation of nodes.
  struct internal_node {
    size_type uid_;
    size_type index_;
    point_type position_;
    node_value_type value_;
  };

  // Struct for proxy representation of edges.
  struct internal_edge {
  	size_type index_;
  	size_type uid1_;
    size_type uid2_;
    edge_value_type value_;
  };

  /**
   * AF:  represents the first index_ that has an edge.
   * RI:  for all i < first_edge_index_, adj_list_[i].size() == 0
   *      adj_list_[first_edge_uid1].size() > 0
   */
  size_type first_edge_index_;

  /**
   *  AF(N) = N = {n_1, n_2, ...n_(num_nodes()-1)}
   *  RI:  internal_nodes_.size() == num_nodes();
   */
  std::vector<internal_node> internal_nodes_;

 /**
   *  Maps node uids to indicies
   *  RI:  n is a valid node -> uid2index[n.uid_] < num_nodes();
   */  std::vector<size_type> uid2index_;

  // Allows resuse of node UIDS.
  std::vector<size_type> empty_node_uids_;

  /** 
   * AF(E) = E = {{n_i, n_j} | index of j in adj_list_[i] && index of i in adj_list_[j]}
   * RI:  adj_list_.size() == num_node(); k nin adj_list_[k]
   */
  std::vector<std::vector<size_type>> adj_list_;

  // Holds data for each edge.
  std::vector<internal_edge> internal_edges_;

  // Allows resuse of internal edge indicies
  std::vector<size_type> empty_edge_indicies_;
};

#endif
