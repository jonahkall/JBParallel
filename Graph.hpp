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
#include <functional>

#include <boost/iterator/transform_iterator.hpp>
#include <boost/iterator/filter_iterator.hpp>

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
  /** Synonym for NodeIterator */
  class NodeIterator;
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

  // Struct to hold data about edges in the adjacency list, contains
  // field to hold edge_value data as well as the uid of the node adjacent
  // to n where we are in the adjacency list of n.uid_.
  struct edge_info {
    mutable edge_value_type value;
    size_type uid;
    bool operator< (const edge_info& ei) const {
      return uid < ei.uid;
    }
    bool operator == (const edge_info& ei) const {
      return (uid == ei.uid && ei.value == value);
    }
    edge_info(const size_type u) {
      uid = u;
    };
  };

  // Define the type of our iterator over each entry in the adjacency list.
  typedef typename std::set<edge_info>::iterator std_edge_iterator;

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

    /** Return this node's position.*/ 
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

    incident_iterator edge_begin() {
      return incident_iterator(g_, uid_ , g_->adj_list_[uid_].begin());
    }

    /** 
     * Returns an incident_iterator to the end of the set of adjacent nodes to the current node. Do not 
     * dereference.
     */
    incident_iterator edge_end() const{
      return incident_iterator(g_, uid_ ,g_->adj_list_[uid_].end());
    }

    incident_iterator edge_end() {
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
    adj_list_.push_back(std::set<edge_info>());
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
    assert (i < i2u_.size());
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

    // Return a reference to the value of the current edge.
    // Note: only adj_list_[min(a.u_index(), b.u_index())] contains the correct
    // value.
    // Complexity: O(d)
    edge_value_type& value () {
      for (auto it = g_->adj_list_[std::min(node1_id_, node2_id_)].begin();
          it != g_->adj_list_[std::min(node1_id_, node2_id_)].end(); ++it) {
        if ((*it).uid == std::max(node1_id_, node2_id_)) {
          return (*it).value;
        }
      }
      // We should never get here: if so our adjacency list is missing entries
      // and we should fail.
      assert(false);
    }

    // Return a (non-modifiable) const reference to the current edge's value.
    // Note: if a and b are the two nodes of an edge, only min(a.u_index(), b.u_index())) will have the
    // up to date value: the value in the other entry into the adjacency list is WRONG.
    // Complexity: O(d)
    const edge_value_type& value () const {
      for (auto it = g_->adj_list_[std::min(node1_id_, node2_id_)].begin();
          it != g_->adj_list_[std::min(node1_id_, node2_id_)].end(); ++it) {
        if ((*it).uid == std::max(node1_id_, node2_id_)) {
          return (*it).value;
        }
      }
      // We should never get here: if so our adjacency list is missing entries
      // and we should fail.
      assert(false);
    }

    /** Return a node of this Edge */
    Node node1() const {
      return Node(g_, node1_id_);
    }

    /** Return the other node of this Edge */
    Node node2() const {
      return Node(g_, node2_id_);
    }

    double length() const {
      return norm(g_->nodes_[node1_id_].position - g_->nodes_[node2_id_].position);
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
      ++num_edges_;
      int i = static_cast<int>(a.u_index());
      size_type j = static_cast<size_type>(b.u_index());
      edge_info ei(j);
      adj_list_[i].insert(ei);
      i = static_cast<int>(b.u_index());
      j = static_cast<size_type>(a.u_index());
      edge_info ei_r(j);
      adj_list_[i].insert(ei_r);
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
      for (size_t k = 0; k < adj_list_[i].size(); ++k) {
        for (auto it = adj_list_[i].begin(); it != adj_list_[i].end(); ++it) {
          if ((*it).uid == j)
            return true;
        }
      }
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
		
		NodeIterator& operator+=(int n){
			node_index_ += n;
			return *this;
		}
		
		NodeIterator& operator-=(int n){
			node_index_ -= n;
			return *this;
		}
		
		bool operator<(const NodeIterator& ni) const{
			if (g_ < ni.g_){
				return true;
			}
			if (g_ == ni.g_ && node_index_ < ni.node_index_){
				return true;
			}
			return false;
		}
		
		NodeIterator& operator[](int n){
			node_index_ = n;
			return *this;
		}
		
<<<<<<< HEAD
		node_iterator operator-(node_iterator a){
			return NodeIterator(this->g_, node_index_ - a->node_index_);
		}
=======
		int operator-(node_iterator a){
			return node_index_ - a.node_index_;
		}

    //node_iterator operator-(node_iterator a){
    //  return NodeIterator(this->g_, node_index_ - a->node_index_);
   // }
>>>>>>> dfa1aa42332cdafc40f019416c9cfdeaae37e5fb
	

   private:
    friend class Graph;
    size_type node_index_;
    Graph* g_;
    NodeIterator(const Graph* g, size_type index = 0) :
        node_index_(index), g_(const_cast<Graph*>(g)){}
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

  // Returns a node_iterator indexed to the first node in the graph.
  node_iterator node_begin() const {
    return NodeIterator(this, 0);
  }

  // Returns a node_iterator indexed to the last node in the graph.
  node_iterator node_end() const {
    return NodeIterator(this, i2u_.size());
  }
	
	friend node_iterator operator+(node_iterator a, int n){
		return a+=n;
	}
	
	friend node_iterator operator+(int n, node_iterator a){
		return a+=n;
	}
	
	friend node_iterator operator-(node_iterator a, int n){
		return a-=n;
	}
	

  struct uid2Node {
    Graph* g_;
    Node operator()(const int& uid) {
      return Node(g_, uid);
    }
    uid2Node(const Graph* g) : g_(const_cast<Graph*>(g)){};
    friend class Graph;
  };

  //std::function<Node()>

  // std::function<Node(size_type)> uid2Node;

  // using NodeIterator = boost::transform_iterator<uid2Node,
  //     std::vector<size_type>::const_iterator>;

  // // Returns a node_iterator indexed to the first node in the graph.
  // NodeIterator node_begin() const {
  //   uid2Node = [] (size_type uid) {return Node(this, uid);};
  //   auto u = boost::make_transform_iterator(i2u_.begin(), uid2Node);
  //   return u;
  // }

  // // Returns a node_iterator indexed to the last node in the graph.
  // NodeIterator node_end() const {
  //   uid2Node = [] (size_type uid) {return Node(this, uid);};
  //   auto u = boost::make_transform_iterator(i2u_.end(), uid2Node);
  //   return u;  
  // }

  //boost::transform_iterator<uid2Node,
   //   std::vector<size_type>::const_iterator>

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
      return Edge(g_, adj_index_, (*set_place_).uid);
    }

    // Increment, ensuring that we end up at an unexplored edge. This invariant
    // is maintained by only visiting an edge if the id of the node we're using
    // to index into the adj list is less than the id of the other node in the edge.
    // So if our adjacency list is: 0: {1,2,3}, 1: {0,3}, we would visit edges in
    // the order: [0,1], [0,2], [0,3], [1,3].
    EdgeIterator& operator++() {
      ++set_place_;
      reset_if_past_end();
      while ((*set_place_).uid < adj_index_ && 
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
    std_edge_iterator set_place_;

    // EdgeIterator private constructor.
    EdgeIterator(const Graph* g, size_type adj_index,
        std_edge_iterator set_place) {
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
      return Edge(g_, curr_node_, (*set_iter_).uid);
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
    std_edge_iterator set_iter_;

    // IncidentIterator private constructor.
    IncidentIterator(const Graph* g, size_type curr_node,
        std_edge_iterator set_iter) :
        g_(const_cast<Graph*>(g)), curr_node_(curr_node), set_iter_(set_iter) {}
  };

  // Extra little helper function for testing code by printing out the current state of the
  // member adjacency list.
  void print_adj_list() {
    for (size_type i = 0; i < adj_list_.size(); ++i) {
      std::cout << "Node  " << i << "adjacent to: " << "\n";
      for (std_edge_iterator it = adj_list_[i].begin();
          it != adj_list_[i].end(); ++it) {
        std::cout << (*it).uid << " ";
      }
      std::cout << "\n";
    }
  }


  /** Remove a node from the graph.
   * @pre @a n is a Node object.
   * @return -1 if we are trying to delete an out-of-range node, 1 if we succeed.
   * @post There is no index i such that @a i2u_[i] = n.uid_
   * @post Invalidates all node and edge iterators except node iterators whose index
   * is less than the index being deleted.
   * @post Invalidates incident iterators to std::set<edge_info> entries in adjacency list
   * at index
   * @post Invalidates all edges with node1() or node2() equal to n
   * @post Node @a n is no longer a valid node in the graph.
   * @post maintains graph invariants.
   *
   *
   * Complexity: O(num_nodes())
   */
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
      for (auto siter = adj_list_[adj_uid].begin(); siter != adj_list_[adj_uid].end(); ++siter) {
        if ((*siter).uid == n.uid_) {
          adj_list_[adj_uid].erase(siter);
          break;
        }
      }
    }
    adj_list_[n.uid_].clear();
    return 1;
  }

  // Wrapper around remove node above which takes a node_iterator, and returns
  // a node_iterator.
  /*node_iterator remove_node(node_iterator n_it) {
    return NodeIterator(this, remove_node(*n_it));
  }*/

  /** Remove an edge from the graph.
   * @pre @a a, @a b are valid nodes in the graph.
   * @return 0 if there is no edge between a and b, 1 upon successful removal.
   * @post adj_list_[a.u_index()] does not contain any edge_info objects with a uid of b.u_index()
   * @post adj_list_[b.u_index()] does not contain any 
   * @post (summary of two conditions above) (a,b) is not a valid edge in the graph.
   * @post Invalidates all edge_iterators pointing within the adjacency list entries at
   * a.u_index() or b.u_index()
   * @post Invalidates incident iterators within adj_list_[b.u_index()] or
   * adj_list_[a.u_index()]
   *
   * Complexity: O(d)
   */
  size_type remove_edge(const Node& a, const Node& b) {
    if (!has_edge(a, b)) {
      return 0;
    }
    --num_edges_;
    auto u1 = a.u_index();
    auto u2 = b.u_index();
    for (std_edge_iterator it = adj_list_[u1].begin(); it != adj_list_[u1].end(); ++it) {
      if ((*it).uid == u2) {
        adj_list_[u1].erase(it);
      }
    }
    for (std_edge_iterator it2 = adj_list_[u2].begin(); it2 != adj_list_[u2].end(); ++it2) {
      if ((*it2).uid == u1) {
        adj_list_[u2].erase(it2);
      }
    }
    return 1;
  }

  // Wrapper around the above remove_edge function which simply
  // returns remove_edge of the two nodes of an edge.
  size_type remove_edge(const Edge& e) {
    return remove_edge(e.node1(), e.node2());
  }

  // Wrapper around remove_Edge which takes an edge iterator and
  // returns an edge_iterator to the beginning of the adjacency list.
  edge_iterator remove_edge(edge_iterator e_it) {
    remove_edge(*e_it);
    return EdgeIterator(this, adj_list_[0].begin(), 0);
  }

 private:
  // Struct for holding info about a node, including position, index into
  // i2u_ (idx), and a value.
  struct node_info {
    Point position;
    size_type idx;
    node_value_type value;
  };
  std::vector<node_info> nodes_;
  std::vector<size_type> i2u_;
  std::vector<std::set<edge_info>> adj_list_;
  size_type num_edges_;
  size_type next_uid_;
};


#endif
