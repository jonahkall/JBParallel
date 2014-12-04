#ifndef CS207_GRAPH_HPP
#define CS207_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <vector>
#include <cassert>
#include <set>
#include <map>
#include <unordered_map>

#include "CS207/Util.hpp"
#include "Point.hpp"

/** @class Graph
 * @brief A template for 3D undirected graphs.
 *
 * Users can add and retrieve nodes and edges. Edges are unique (there is at
 * most one edge between any pair of distinct nodes).
 */
template <typename V, typename E>
class Graph {
 public:

  typedef V node_value_type;
  typedef E edge_value_type;

  /////////////////////////////
  // PUBLIC TYPE DEFINITIONS //
  /////////////////////////////

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
  Graph() {
    next_uid_ = 0;
    num_edges_ = 0;
  }
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
    return i2u_.size();
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    nodes_.clear();
    adj_list_.clear();
    next_uid_ = 0;
    num_edges_ = 0;
    i2u_.clear();
    edge_values_.clear();
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
    Node() {
      // We use these invalid values by convention.
      g_ = NULL;
      uid_ = -1;
    }

    int parent_size() const {
      return g_->i2u_.size();
    }

    /** Return this node's position. */
    const Point& position() const {
      return g_->nodes_[uid_].position;
    }

    Point& position() {
      return g_->nodes_[uid_].position;
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      return g_->nodes_[uid_].idx;
    }

    size_type u_index() const {
      assert(uid_ < g_->nodes_.size());
      return uid_;
    }
    
    /** Returns true if g2_ is the same as the current internal graph
     * pointer.
     */
    bool from_same_graph(Graph* g2_) const {
      return g2_ == g_;
    }

    /** Test whether this node and @a x are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& x) const {
       // Make sure the nodes have the same index and are from the same graph.
      if (uid_ != x.u_index() || !x.from_same_graph(g_)) 
        return false;
      return true;
    }

    /** Test whether this node is less than @a x in the global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any geometric meaning.
     *
     * The node ordering relation must obey trichotomy: For any two nodes x
     * and y, exactly one of x == y, x < y, and y < x is true.
     */
    bool operator<(const Node& x) const {
      if (g_ < x.g_) {
        return true;
      }
      return (x.from_same_graph(g_) && uid_ < x.u_index());
    } 
    /** 
     * Returns a reference to the current node's value, allowing modification.
     */
    node_value_type& value() {
      return g_->nodes_[uid_].value;
    };

    /** 
     * Returns a const reference to the current node's value, not allowing modification.
     */
    const node_value_type& value() const {
      return g_->nodes_[uid_].value;
    };

    /** 
     * Returns the number of adjacent nodes to the current node. O(1).
     */
    size_type degree() const {
      return g_->adj_list_[uid_].size();
    }

    /** 
     * Returns an incident_iterator to the first element of the set of adjacent nodes to the current node.
     */
    incident_iterator edge_begin() const {
      return incident_iterator(g_, uid_ , g_->adj_list_[uid_].begin());
    }

    /** 
     * Returns an incident_iterator to the end of the set of adjacent nodes to the current node. Do not 
     * dereference.
     */
    incident_iterator edge_end() const{
      return incident_iterator(g_, uid_ ,g_->adj_list_[uid_].end());
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    /** The private constructor for node -- can only be used by graph, sets
    * graph pointer to "g" and uid_ to "uid"
    */
    
    Node(const Graph* g, size_type uid) :
      uid_(uid), g_(const_cast<Graph*>(g)) {
    }

    size_type uid_;
    Graph* g_;
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
   * Complexity: O(1) amortized operations.
   */
  Node add_node(const Point& position, const node_value_type& value = node_value_type()) {
    node_info ni;
    ni.position = position;
    ni.value = value;
    ni.idx = i2u_.size();
    nodes_.push_back(ni);
    i2u_.push_back(next_uid_);
    adj_list_.push_back(std::set<size_type>());
    ++next_uid_;
    return Node(this, next_uid_ - 1);
  }

  /** Determine if this Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    return n.from_same_graph(this);
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    return Node(this, i2u_[i]);
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
  class Edge : private totally_ordered<Edge>{
   public:
    /** Construct an invalid Edge. */
    Edge() {
      node1_id_ = -1;
      node2_id_ = -1;
      edge_id_ = -1;
      g_ = NULL;
    }

    // ADD SPEC
    edge_value_type& value () {
      auto ord_pair = std::make_pair(std::min(node1_id_, node2_id_), std::max(node1_id_, node2_id_));
      return g_->edge_values_[ord_pair];
    }

    // ADD SPEC
    const edge_value_type& value () const {
      auto ord_pair = std::make_pair(std::min(node1_id_, node2_id_), std::max(node1_id_, node2_id_));
      return g_->edge_values_[ord_pair];
    }

    /** Return a node of this Edge */
    Node node1() const {
      return Node(g_, node1_id_);
    }

    /** Return the other node of this Edge */
    Node node2() const {
      return Node(g_, node2_id_);
    }

    /** Returns true if g and the edge's current internal 
     * graph pointer are the same.
     */
    bool from_same_graph(Graph* g) const {
      return g_ == g;
    }
    
    /** Test whether this edge and @a x are equal.
    *
    * Equal edges are from the same graph and have the same nodes.
    */
    bool operator==(const Edge& x) const {
      if (!x.from_same_graph(g_)){
        return false;
      }
      else {
        return (std::min(x.node1_id_, x.node2_id_) ==
            std::min(node1_id_, node2_id_) && std::max(x.node1_id_, x.node2_id_)
                == std::max(node1_id_, node2_id_));
      }
      return true;
    }

    /** Test whether this edge is less than @a x in the global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any geometric meaning.
     *
     * The edge ordering relation must obey trichotomy: For any two edges x
     * and y, exactly one of x == y, x < y, and y < x is true.
     */
    bool operator<(const Edge& x) const {
      if (g_ < x.g_) {
        return true;
      }

      size_type x_1 = std::min(x.node1_id_, x.node2_id_);
      size_type x_2 = std::max(x.node1_id_, x.node2_id_);

      // Compare the first node ids.
      if (std::min(node1_id_, node2_id_) != x_1) {
        return std::min(node1_id_, node2_id_) < x_1;
      }
      // If those are equal, compare the second node ids.
      else if (std::max(node1_id_, node2_id_) != x_2) {
        return std::max(node1_id_, node2_id_) < x_2;
      }
      // Both nodes have equal ids.
      else {
        return false;
      }
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    
    /** Private constructor so that Graph can construct edges. */
    Edge(const Graph* g, size_type node1_id, 
         size_type node2_id) {
      assert(node1_id != node2_id);
      g_ = const_cast<Graph*> (g);
      node1_id_ = node1_id;
      node2_id_ = node2_id;
    }

    Graph* g_;
    size_type node1_id_;
    size_type node2_id_;
    size_type edge_id_;
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: O(1)
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
   * equal new edge(@a i). Must not invalidate outstanding Edge objects.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge add_edge(const Node& a, const Node& b) {
    assert(a.u_index() != b.u_index());
    if (has_edge(a,b) || has_edge(b,a)) {
      return Edge(this, a.u_index(), b.u_index());
    }
    else {
      std::pair<size_type, size_type> to_insert;
      to_insert.first = std::min(a.u_index(), b.u_index());
      to_insert.second = std::max(a.u_index(), b.u_index());
      edge_values_.insert(std::make_pair(to_insert, edge_value_type()));
      ++num_edges_;
      int i = static_cast<int>(a.u_index());
      size_type j = static_cast<size_type>(b.u_index());
      adj_list_[i].insert(j);
      i = static_cast<int>(b.u_index());
      j = static_cast<size_type>(a.u_index());
      adj_list_[i].insert(j);
      return Edge(this, a.u_index(), b.u_index());
    }
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return true if, for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    int i = static_cast<int>(a.u_index());
    size_type j = static_cast<size_type>(b.u_index());
    if (static_cast<size_type>(i) < adj_list_.size()) {
      if (std::find(adj_list_[i].begin(), adj_list_[i].end(), j) != adj_list_[i].end())
        return true;
    }
    return false;
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: O(num_edges())
  */
  Edge edge(size_type i) const {
    edge_iterator it = edge_begin();
    for (; i != 0; --i) {
      ++it;
    }
    return *it;
  }


  ///////////////
  // Iterators //
  ///////////////

  /** @class Graph::NodeIterator
   * @brief Iterator class for nodes. A forward iterator. 
   * Really just a wrapper for an index: since the nodes are stored in a vector,
   * iterating is easy.
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
    NodeIterator() {
      g_ = NULL;
      node_index_ = -1;
    }

    // Returns the current index into the node_info vector.
    size_type index () const {
      return node_index_;
    }

    // Dereference.  Returns a Node with the current index on the fly.
    Node operator*() const {
      return Node(g_, g_->i2u_[node_index_]);
    }

    // Increment the NodeIterator.  Simply increments the current node index.
    // Returns the current node iterator after increment.
    NodeIterator& operator++() {
      ++node_index_;
      return *this;
    }

    NodeIterator& operator--() {
      assert(node_index_ > 0);
      --node_index_;
      return *this;
    }

    // Checks if two node iterators are equal by checking that they have the
    // same graph pointer and the same node index.
    bool operator==(const NodeIterator& ni) const {
      return (ni.index() == node_index_ && ni.g_ == g_);
    }

   private:
    friend class Graph;
    size_type node_index_;
    Graph* g_;
    NodeIterator(const Graph* g, size_type index = 0) :
        node_index_(index), g_(const_cast<Graph*>(g)){}
  };

  // Returns a node_iterator indexed to the first node in the graph.
  node_iterator node_begin() const {
    return NodeIterator(this, 0);
  }

  // Returns a node_iterator indexed to the last node in the graph.
  node_iterator node_end() const {
    return NodeIterator(this, i2u_.size());
  }

  /** @class Graph::EdgeIterator
   * @brief Iterator class for edges. A forward iterator. */
  class EdgeIterator : private totally_ordered<EdgeIterator> {
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
    EdgeIterator() {
      g_ = NULL;
      adj_index_ = -1;
    }

    // Helper function which checks if the iterator is at the end of the
    // current set in the adjacency list, and if so advances it to the start of
    // the next one.
    void reset_if_past_end() {
      if (adj_index_ == g_->adj_list_.size() -1)
        return;
      while (set_place_ == g_->adj_list_[adj_index_].end()) {
        ++adj_index_;
        set_place_ = g_->adj_list_[adj_index_].begin();
      }
    }

    // Dereference, returns an edge on the fly.
    Edge operator*() const {
      return Edge(g_, adj_index_, *set_place_);
    }

    // Increment, ensuring that we end up at an unexplored edge. This invariant
    // is maintained by only visiting an edge if the id of the node we're using
    // to index into the adj list is less than the id of the other node in the edge.
    // So if our adjacency list is: 0: {1,2,3}, 1: {0,3}, we would visit edges in
    // the order: [0,1], [0,2], [0,3], [1,3].
    EdgeIterator& operator++() {
      ++set_place_;
      reset_if_past_end();
      while (*set_place_ < adj_index_ && 
          set_place_ != g_->adj_list_[g_->adj_list_.size() - 1].end()) {
        ++set_place_;
        reset_if_past_end();
      }
      return *this;
    }

    // Check if the current EdgeIterator is equal to another EdgeIterator by checking
    // equality of all private members: graph pointer, index into adj list, iterator
    // into current set.
    bool operator==(const EdgeIterator& ei) const {
      return (ei.g_ == g_ && adj_index_ == ei.adj_index_ && set_place_ == ei.set_place_);
    }

   private:
    friend class Graph;
    Graph* g_;
    size_type adj_index_;
    std::set<size_type>::iterator set_place_;

    // EdgeIterator private constructor.
    EdgeIterator(const Graph* g, size_type adj_index,
        std::set<size_type>::iterator set_place) {
      g_ = const_cast<Graph*>(g);
      adj_index_ = adj_index;
      set_place_ = set_place;
      while (g_->adj_list_[adj_index_].size() == 0) {
        ++adj_index_;
      }
    }
  };

  // Returns the first edge in the graph which is valid to explore, that is
  // adj_index < *set_place_
  edge_iterator edge_begin() const {
    auto ei = EdgeIterator(this, 0, adj_list_[0].begin());
    return ei;
  };

  // Returns an edge_iterator to the end of the graph's edges in
  // adjacency list format.
  edge_iterator edge_end() const {
    return EdgeIterator(this, adj_list_.size() - 1, adj_list_[adj_list_.size() - 1].end());
  };


  /** @class Graph::IncidentIterator
   * @brief Iterator class for edges incident to a node. A forward iterator.
   * Because adjacent nodes are stored in a std::set in the adjacency list,
   * this iterator is essentially a wrapper for the std::set iterator.
   */
  class IncidentIterator : private totally_ordered<IncidentIterator>{
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
    IncidentIterator() {
      g_ = NULL;
      curr_node_ = -1;
    }

    // Returns an edge on the fly with the current graph pointer, and
    // appropriate node indices.
    Edge operator*() const {
      return Edge(g_, curr_node_, *set_iter_);
    }

    // Increments by incrementing the set iterator.
    IncidentIterator& operator++() {
      ++set_iter_;
      return *this;
    }

    // Checks for equality between IncidentIterators by checking equality of all
    // private members.
    bool operator==(const IncidentIterator& i_i) const {
      return (i_i.set_iter_ == set_iter_ && i_i.g_ == g_ &&
          i_i.curr_node_ == curr_node_);
    }

   private:
    friend class Graph;
    Graph* g_;
    size_type curr_node_;
    std::set<size_type>::iterator set_iter_;

    // IncidentIterator private constructor.
    IncidentIterator(const Graph* g, size_type curr_node,
        std::set<size_type>::iterator set_iter) :
        g_(const_cast<Graph*>(g)), curr_node_(curr_node), set_iter_(set_iter) {}
  };

  // Extra little helper function for testing code by printing out the current state of the
  // member adjacency list.
  void print_adj_list() {
    for (size_type i = 0; i < adj_list_.size(); ++i) {
      std::cout << "Node  " << i << "adjacent to: " << "\n";
      for (std::set<size_type>::iterator it = adj_list_[i].begin();
          it != adj_list_[i].end(); ++it) {
        std::cout << *it << " ";
      }
      std::cout << "\n";
    }
  }

  int remove_node(const Node& n) {
    if (n.index() >= i2u_.size()) {
      return -1;
    }
    auto it = i2u_.erase(i2u_.begin() + nodes_[n.uid_].idx);
    nodes_[n.uid_].idx = -1;
    for (; it != i2u_.end(); ++it) {
      --nodes_[*it].idx;
    }
    num_edges_ -= adj_list_[n.uid_].size();
    for (auto ii_t = n.edge_begin(); ii_t != n.edge_end(); ++ii_t) {
      size_type adj_uid = (*ii_t).node2().uid_;
      adj_list_[adj_uid].erase(adj_list_[adj_uid].find(n.uid_));
    }
    adj_list_[n.uid_].clear();
    return 1;
  }

  node_iterator remove_node(node_iterator n_it) {
    return NodeIterator(this, remove_node(*n_it));
  }

  size_type remove_edge(const Node& a, const Node& b) {
    if (!has_edge(a, b)) {
      return 0;
    }
    --num_edges_;
    auto u1 = a.u_index();
    auto u2 = b.u_index();
    adj_list_[u1].erase(adj_list_[u1].find(u2));
    adj_list_[u2].erase(adj_list_[u2].find(u1));
    auto value_pair = std::make_pair(std::min(u1,u2), std::max(u1,u2));
    edge_values_.erase(edge_values_.find(value_pair));
    return 1;
  }

  size_type remove_edge(const Edge& e) {
    return remove_edge(e.node1(), e.node2());
  }

  edge_iterator remove_edge(edge_iterator e_it) {
    size_type adj_ind = (*e_it).node1().u_index();
    remove_edge(*e_it);
    return EdgeIterator(this, adj_ind, 0);
  }

 private:
  // Struct for holding info about a node, including position, index into
  // nodes_ (uid), and a value.
  struct node_info {
    Point position;
    size_type idx;
    node_value_type value;
  };
  std::vector<node_info> nodes_;
  std::vector<size_type> i2u_;
  std::vector<std::set<size_type>> adj_list_;
  size_type num_edges_;
  size_type next_uid_;
  std::map<std::pair<size_type, size_type>, edge_value_type> edge_values_;
};

#endif
