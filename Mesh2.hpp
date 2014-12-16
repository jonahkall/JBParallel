/**
* CS207 HW4
* Editors: Wenshuai Ye && Yuhao Zhu 
**/

#include <algorithm>
#include <vector>
#include <cassert>
#include <cmath>

#include "CS207/Util.hpp"
#include "Point.hpp"
using namespace std;

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
class Mesh{
public:

	/////////////////////////////
	// PUBLIC TYPE DEFINITIONS //
	/////////////////////////////

	typedef N node_value_type;

	typedef E edge_value_type;

	typedef T triangle_value_type;

	/** Type of this graph. */
	typedef Mesh mesh_type;

	/** Predeclaration of Node type. */
	class Node;
	/** Synonym for Node (following STL conventions). */
	typedef Node node_type;

	/** Predeclaration of Edge type. */
	class Edge;
	/** Synonym for Edge (following STL conventions). */
	typedef Edge edge_type;

	/** Predeclaration of Triangle type*/
	class Triangle;
	/** Synonym for Edge (following STL conventions). */
	typedef Triangle triangle_type;

	/** Type of indexes and sizes.
	Return type of Graph::Node::index(), Graph::num_nodes(),
	Graph::num_edges(), and argument type of Graph::node(size_type) */
	typedef unsigned size_type;

	/** Type of node iterators, which iterate over all graph nodes. */
	class node_iterator;

	/** Type of edge iterators, which iterate over all graph edges. */
	class edge_iterator;

	/** Type of incident iterators, which iterate incident edges to a node. */
	class IncidentIterator;

	/** Synonym for IncidentIterator */
	typedef IncidentIterator incident_iterator;

	/** Type of triangle iterators, which iterate over all triangles. */
	class triangle_iterator;

	/** Type of triangle incident iterators, 
	which iterate incident triangles to a node. */
	class nincident_tri_iterator;


	/** Type of triangle incident iterators, 
	which iterate incident triangles to a triangle. */
	class eincident_tri_iterator;


	/** Return the number of nodes in the mesh.
	*
	* Complexity: O(1).
	*/
	size_type size() const {
		return node_i2u.size();
	}

	/** Remove all nodes and edges from this mesh.
	* @post num_nodes() == 0 && num_edges() == 0
	*
	* Invalidates all outstanding Node and Edge objects.
	*/
	void clear() {
		nodes_.clear();
		edges_.clear();
		node_i2u.clear();
		edge_values.clear();
		edge_size_=0;
	}


	Mesh(): edge_size_(0),nodes_(), edges_(), triangles_() {}

	~Mesh() = default;

	class Node: private totally_ordered<Node>{
	public:
		Node() {
		}
      
		/** Return this Node's position. 
		* Complexity: O(1)
		*/
		Point& position() const {
			return m_->nodes_[node_uid_].points;
		}

		/** Return this node's index, a number in the range [0, mesh_size). 
		* Complexity: O(1)
		*/
		size_type index() const {
			return m_->nodes_[node_uid_].index;
		}

		/** Test whether this node and @a x are equal.
		*Return true if @a x and this node have the same mesh and the same uid.
		* Complexity: O(1)
		*/
		bool operator==(const Node& x) const {
			if (node_uid_ == x.node_uid_ && m_ == x.m_) return true;
			(void) x;
			return false;
		}

		/** Test whether this node is less than @a x in the global order.
		*Return true if uid < @a x.uid
		* Complexity: O(1)
		*/
		bool operator<(const Node& x) const {
			if (node_uid_ < x.node_uid_) return true;
			(void) x;
			return false;
		}

		/** Return this Node's value
		* Complexity: O(1)
		*/
		node_value_type& value(){
			return m_-> nodes_[node_uid_].value;
		}

		/** Return this Node's value
		* Complexity: O(1)
		*/
		const node_value_type& value() const{
			return m_->nodes_[node_uid_].value;
		}

		/** Return the number of incident edges it has to the node
		*
		* Complexity O(1)
		*/
		size_type edge_degree() const {
			return m_->edges_[node_uid_].size();
		}

		/* Return the incident_iterator to the current node's first incident edge
		* Complexity O(1)
		**/
		incident_iterator edge_begin() const {
			return incident_iterator(m_, node_uid_, 0);
		}

		/* Return the incident_iterator to the current node's 
		*    last incident edge (invalid one)
		* Complexity O(1)
		**/
		incident_iterator edge_end() const {
			return incident_iterator(m_, node_uid_, m_->edges_[node_uid_].size());
		}

		/** Return the number of incident triangles it has to the node
		*
		* Complexity O(1)
		*/
		size_type triangle_degree() const{
			return m_->nodes_[node_uid_].adj_triangles_.size();
		}
      
		/* Return the nincident_tri_iterator to the current node's 
		* first incident triangle
		* Complexity O(1)
		**/
		nincident_tri_iterator triangle_begin() const{
			return nincident_tri_iterator(m_,node_uid_,0);
		}

		/* Return the nincident_tri_iterator to the current node's last
		* incident triangle (invalid triangle)
		* Complexity O(1)
		**/
		nincident_tri_iterator triangle_end() const{
			return nincident_tri_iterator(m_,node_uid_,m_->nodes_[node_uid_].adj_triangles_.size());
		}

	private:
		// Allow Graph to access Node's private member data and functions.
		Node(const Mesh* m, size_type uid) 
			: node_uid_(uid), m_(const_cast<Mesh*>(m)){}
		friend class Mesh;
		size_type node_uid_;
		Mesh *m_;
	};

	/** Return the number of nodes in the mesh.
	* Complexity: O(1).
	*/
	size_type num_nodes() const {
		return size();
	}

	/** Determine if this Node belongs to this mesh
	* @return True if @a n is currently a Node of this mesh
	* Complexity: O(1).
	*/
	bool has_node(const Node& n) const {
		return (n.m_ == this && node_i2u[nodes_[n.node_uid_].index] == n.node_uid_);
	}

	/** Add a node to the mesh, returning the added node.
	* @param[in] position The new node's position
	* @param[in] value The value corresponding to the node
	* @post new size() == old size() + 1
	* @post result_node.index() == old size()
	* @return the node with corresponding input position
	* Complexity: O(1) amortized operations.
	*/
	Node add_node(const Point& position, const node_value_type& value = node_value_type()) {
		node_elem new_elem;
		new_elem.value = value;
		new_elem.points = position;
		new_elem.index = node_i2u.size();
		node_i2u.push_back(nodes_.size());
		nodes_.push_back(new_elem);
		edges_.push_back(vector<edge_elem>());
		return Node(this, nodes_.size() - 1);
		(void) position;
	}

	/** Return the node with index @a i.
	* @pre 0 <= @a i < num_nodes()
	* @post result_node.index() == i
	* Complexity: O(1).
	*/
	Node node(size_type i) const {
		return Node(this, node_i2u[i]);
		(void) i;
	}

	class Edge: private totally_ordered<Edge>{
	public:
		/** Construct an invalid Edge. */
		Edge(){}
       
		/** Return a node of this Edge */
		Node node1() const {
			return Node(m_, uid1_);
		}

		/** Return the other node of this Edge */
		Node node2() const {
			return Node(m_, uid2_);
		}

		/** Return this Edge's value
		* Complexity: O(uid1.edge_degree()) worst case
		*/
		edge_value_type& value(){
			size_type edge_index_;
			for (size_type i =0; i != m_->edges_[uid1_].size(); ++i){
				if (m_->edges_[uid1_][i].uid == uid2_){
					edge_index_ = m_->edges_[uid1_][i].edge_index;
					return m_->edge_values[edge_index_];
				}
			}
		}

		/** Return this Edge's value
		* Complexity: O(uid1.edge_degree()) worst case
		*/
		const edge_value_type& value() const{
			size_type edge_index_;
			for (size_type i =0; i != m_->edges_[uid1_].size(); ++i){
				if (m_->edges_[uid1_][i].uid == uid2_){
					edge_index_ = m_->edges_[uid1_][i].edge_index;
					return m_->edge_values[edge_index_];
				}
			}
		}

		/** Return whether this edge and @a x are equal (they are from the same 
		* mesh and have the same nodes).
		* Complexity: O(1)
		*/
		bool operator==(const Edge& x) const {
			return (m_ == x.m_ && uid1_ == x.uid1_ && uid2_ == x.uid2_);
			(void) x;
			return false;
		}

		/** Return whether this edge is less than @a x in the global order 
		*  (if uid is  less than de uid of @a x).
		* Complexity: O(1)
		*/
		bool operator<(const Edge& x) const {
			if (m_ != x.m_) return true;
			else return (uid1_<x.uid1_ && uid2_<x.uid2_);
			(void) x;
		}

		/** Return this Edge's length.
		* Complexity: O(1)
		*/
		double length() const{
			return norm(node1().position() - node2().position());
		}

		/* Return the eincident_tri_iterator to the current 
		*    triangle's first incident triangle
		* Complexity O(adj_edges) worst case
		**/
		eincident_tri_iterator triangle_begin() const{
			size_type uid2_index_;
			for (size_type i =0; i != m_->edges_[uid1_].size(); ++i){
				if (m_->edges_[uid1_][i].uid == uid2_)
					uid2_index_ = i;
			}
			return eincident_tri_iterator(m_,uid1_,uid2_index_,0);
		}

		/* Return the eincident_tri_iterator to the current 
		*    triangle's first incident triangle
		* Complexity O(adj_edges/2) on average
		**/
		eincident_tri_iterator triangle_end() const{
			size_type uid2_index_;
			for (size_type i =0; i != m_->edges_[uid1_].size(); ++i){
				if (m_->edges_[uid1_][i].uid == uid2_)
					uid2_index_ = i;
			}
			return eincident_tri_iterator(m_,uid1_,uid2_index_,m_->edges_[uid1_][uid2_index_].adj_triangles_.size());
		}

		/** Return the number of incident triangles it has to 
		*      the current edge
		* Complexity O(1)
		*/
		size_type degree() const{
			size_type uid2_index_;
			for (size_type i =0; i != m_->edges_[uid1_].size(); ++i){
				if (m_->edges_[uid1_][i].uid == uid2_)
					uid2_index_ = i;
			}
			return m_->edges_[uid1_][uid2_index_].adj_triangles_.size();
		}

		Triangle eAdj_Triangle(size_type i) const{
			assert(i>=1 && i<=2);
			size_type index_;
			for (size_type j =0; j != m_->edges_[uid1_].size(); ++j){
				if (m_->edges_[uid1_][j].uid == uid2_)
					index_ = j;
			}	
			return Triangle(m_,m_->edges_[uid1_][index_].adj_triangles_[i-1]);
		}

	private:
		// Allow Graph to access Edge's private member data and functions.
		friend class Mesh;
		size_type uid1_;
		size_type uid2_;
		Edge(const Mesh* m,size_type uid1, size_type uid2) : uid1_(uid1), uid2_(uid2), m_(const_cast<Mesh*>(m)){}
		Mesh* m_;
	};

	/** Return the number of edges this mesh has. 
	* Complexity: O(1)
	*/
	size_type num_edges() const {
		return edge_size_;
	}

	/** Test whether two nodes are connected by an edge.
	* @param[in] a, b The nodes that might be connected by an edge
	* @pre @a a and @a b are valid nodes of this mesh
	* @return true if, for some @a i, edge(@a i) connects @a a and @a b.
	*
	* Complexity: No more than O(a.degree()) 
	*/
	bool has_edge(const Node& a, const Node& b) const {
		for (size_type i=0; i<edges_[a.node_uid_].size(); ++i){
			if (b.node_uid_ == edges_[a.node_uid_][i].uid) return true;
		}
		(void) a; (void) b;
		return false;
	}

	/** Add an edge to the mesh, or return the current edge if it already exists.
	* @pre @a a and @a b are distinct valid nodes of this mesh
	* @param[in] a, b The nodes that might or might not be connected by an edge
	* @param[in] value The corresponding value to the edge
	* @post has_edge(@a a, @a b) == true
	* @post If old has_edge(@a a, @a b), new num_edges() == old num_edges().
	*       Else,                        new num_edges() == old num_edges() + 1.
	*
	* @return an Edge object e with e.node1() == @a a and e.node2() == @a b
	* Complexity: O(1) amortized time
	*/
	Edge add_edge(const Node& a, const Node& b, const edge_value_type& value = edge_value_type()) {
		if (!has_edge(a,b)){
			edge_elem elem1;
			edge_elem elem2;
			elem1.uid = a.node_uid_;
			elem1.edge_index = edge_size_;
			elem2.uid = b.node_uid_;
			elem2.edge_index = edge_size_;
			edges_[a.node_uid_].push_back(elem2);
			edges_[b.node_uid_].push_back(elem1);
			edge_values.push_back(value);
			++edge_size_;
		}
		return Edge(this, a.node_uid_, b.node_uid_);
		(void) a, (void) b;
	}

	/** Return the node with index @a i.*/
	Edge edge(size_type i) const{
		edge_iterator it=edge_begin();
		for (;i!=0;--i) ++it;
		return *it;
	}

	class Triangle: private totally_ordered<Triangle>{
	public:
		/** Construct an invalid Edge. */
		Triangle(){}

		/** Return a node of this triangle based on the predefined 
		*    index within a triangle
		* @param[in] i The index of the Node.
		* @pre 0<i<=3
		* @post Node.node_uid_==m_->triangles_[tri_uid_].node_uidi_
		*
		* Complexity: O(1)
		* triangleNode(i) is the node opposite to the edge triangleEdge(i)
		*/
		Node triangle_node(size_type i) const{
			assert(0<i && i<=3);
			if (i==1) return Node(m_,m_->triangles_[tri_uid_].node_uid1);
			else if (i==2) return Node(m_,m_->triangles_[tri_uid_].node_uid2);
			else return Node(m_,m_->triangles_[tri_uid_].node_uid3);
		}

		/** Return an edge of this triangle based on the predefined 
		*    index within a triangle
		* @param[in] i The index of the edge.
		* @pre 0<i<=3
		* @post Node.edge_uid_==m_->triangles_[tri_uid_].edge_uidi_
		*
		* Complexity: O(1)
		* triangleEdge(i) is the edge opposite to the edge triangleNode(i)
		*/
		Edge triangle_edge(size_type i) const{
			assert(0<i && i <= 3);
			if (i==1) return Edge(m_,m_->triangles_[tri_uid_].node_uid2,m_->triangles_[tri_uid_].node_uid3);
			else if (i==2) return Edge(m_,m_->triangles_[tri_uid_].node_uid1,m_->triangles_[tri_uid_].node_uid3);
			else return Edge(m_,m_->triangles_[tri_uid_].node_uid1,m_->triangles_[tri_uid_].node_uid2);
		}

		/** Return this Triangle's area
		* Complexity: O(1)
		*/
		double area(){
			double length1 = triangle_edge(1).length();
			double length2 = triangle_edge(2).length();
			double length3 = triangle_edge(3).length();
			double s = (length1 + length2 + length3)/2;
			return sqrt(s*(s-length1)*(s-length2)*(s-length3));
		}

		/** Return the outward normal vector of this triangle 
		*      corresponding to an edge
		* @param[i] i The index of the edge
		* @pre 0<i<=3
		* Complexity: O(1)
		*/
		Point normal_vector(size_type i){
			Point p = triangle_edge(i).node1().position() - triangle_edge(i).node2().position();
			double nx = p.y;
			double ny = -p.x;
			Point checkP = triangle_edge(i).node1().position() - triangle_node(i).position();
			double ln = sqrt(nx*nx + ny*ny);
			double le = triangle_edge(i).length();
			if (nx*checkP.x + ny*checkP.y<0){
				nx = -nx/ln*le;
				ny = -ny/ln*le;
			}
			else{
				nx = nx/ln*le;
				ny = ny/ln*le;
			}
			return Point(nx,ny,0);
		}

		/** Return this Triangle's value
		* Complexity: O(1)
		*/
		triangle_value_type& value(){
			return m_->triangles_[tri_uid_].value;
		}

		/** Return this Triangle's value
		* Complexity: O(1)
		*/
		const triangle_value_type& value() const{
			return m_->triangles_[tri_uid_].value;
		}

		/** Return whether this triangle and @a t are equal 
		*  (they are from the same mesh and have the same nodes).
		* Complexity: O(1)
		*/
		bool operator==(const Triangle& t) const {
			return (m_==t.m_ && tri_uid_ == t.tri_uid_);
		}

		/** Return the adjacent triangle associated with the index.
		* @param[in] i the index of the adjacent triangle
		* @pre 0 < i <= 3
		*
		* Complexity: O(1) amortized time*/

		Triangle triAdj_triangle(size_type i) const{
			assert(i>0 && i<=3);
			if (triangle_edge(i).degree() == 1) return triangle_edge(i).eAdj_Triangle(1);
			else{
				if (triangle_edge(i).eAdj_Triangle(1) == *this)
					return triangle_edge(i).eAdj_Triangle(2);
				else
					return triangle_edge(i).eAdj_Triangle(1);
			}
		}
      
	private:
		friend class Mesh;
		Mesh* m_; 
		size_type tri_uid_; 
        
		// Private constructor
		Triangle(const Mesh* mesh, size_type uid) : m_(const_cast<Mesh*>(mesh)), tri_uid_(uid) {}
	};

	/** Return the total number of triangles in the mesh.
	* Complexity: O(1)
	*/
	size_type num_triangles() const {
		return triangles_.size();
	}

	/** Test whether there exist a triangle that contains the 3 nodes.
	* @param[in] a,b,c The input nodes that might have formed a triangle
	* @pre @a a, @a b and @a c are valid nodes of the mesh
	* @return True if, for some @a i, there exists a triangle(@a i) 
	*    whose nodes are @a a, @a b and @a c.
	*
	* Complexity: O(num_triangles()) worst case
	*/
	bool has_triangle(const Node& a, const Node& b, const Node& c) const {
		assert(has_node(a) && has_node(b) && has_node(c));
		for (size_type t;t<triangles_.size();++t){
			if (a.node_uid_ == triangles_[t].node_uid1 && 
				b.node_uid_ == triangles_[t].node_uid2 && 
					c.node_uid_ == triangles_[t].node_uid3) return true;
			if (a.node_uid_ == triangles_[t].node_uid2 && 
				b.node_uid_ == triangles_[t].node_uid1 && 
					c.node_uid_ == triangles_[t].node_uid3) return true;
			if (a.node_uid_ == triangles_[t].node_uid3 && 
				b.node_uid_ == triangles_[t].node_uid1 && 
					c.node_uid_ == triangles_[t].node_uid2) return true;
			if (a.node_uid_ == triangles_[t].node_uid1 && 
				b.node_uid_ == triangles_[t].node_uid3 && 
					c.node_uid_ == triangles_[t].node_uid2) return true;
			if (a.node_uid_ == triangles_[t].node_uid2 && 
				b.node_uid_ == triangles_[t].node_uid3 && 
					c.node_uid_ == triangles_[t].node_uid1) return true;
			if (a.node_uid_ == triangles_[t].node_uid3 && 
				b.node_uid_ == triangles_[t].node_uid2 && 
					c.node_uid_ == triangles_[t].node_uid1) return true;
		}
		return false;
	}

	/** Add a triangle to the mesh, or return the current triangle 
	*    if it already exists.
	* @param[in] a,b,c The input nodes that might have formed a triangle
	* @param[in] value The value corresponding to the triangle
	* @pre @a a, @a b and @a c are distinct valid nodes of this mesh 
	*    (They are not on one line)
	* @post has_triangle(@a a, @a b, @a c) == true
	* @post If old has_triangle(@a a, @a b, @a c), 
	*          new num_triangles() == old num_triangles();
	*       Else, new num_triangles() == old num_triangles() + 1 
	*          and result_node.uid() == old num_triangles()
	*
	* @return a Triangle t with t.node1(), t.node2(), t.node3()==@a a, @a b, @a c
	* Complexity: O(a.edge_degree() + b.edge_deree() + c.edge_degree()) worst case
	*/
	Triangle add_triangle(const Node& a, const Node& b, const Node& c, const triangle_value_type& value = triangle_value_type()) {
		if (!has_triangle(a,b,c)){
			nodes_[a.node_uid_].adj_triangles_.push_back(num_triangles());
			nodes_[b.node_uid_].adj_triangles_.push_back(num_triangles());
			nodes_[c.node_uid_].adj_triangles_.push_back(num_triangles());

			add_edge(a,b);
			add_edge(b,c);
			add_edge(a,c);
  	  
			size_type indexac;
			size_type indexab;
			size_type indexbc;
			size_type indexca;
			size_type indexba;
			size_type indexcb;

			for (size_type i =0; i != edges_[a.node_uid_].size(); ++i){
				if (edges_[a.node_uid_][i].uid == b.node_uid_)
					indexab = i;
				if (edges_[a.node_uid_][i].uid == c.node_uid_)
					indexac = i;
			}

			for (size_type i =0; i != edges_[b.node_uid_].size(); ++i){
				if (edges_[b.node_uid_][i].uid == c.node_uid_)
					indexbc = i;
				if (edges_[b.node_uid_][i].uid == a.node_uid_)
					indexba = i;
			}

			for (size_type i =0; i != edges_[c.node_uid_].size(); ++i){
				if (edges_[c.node_uid_][i].uid == a.node_uid_)
					indexca = i;
				if (edges_[c.node_uid_][i].uid == b.node_uid_)
					indexcb = i;
			}

			edges_[a.node_uid_][indexab].adj_triangles_.push_back(num_triangles());
			edges_[a.node_uid_][indexac].adj_triangles_.push_back(num_triangles());
			edges_[b.node_uid_][indexba].adj_triangles_.push_back(num_triangles());
			edges_[b.node_uid_][indexbc].adj_triangles_.push_back(num_triangles());
			edges_[c.node_uid_][indexca].adj_triangles_.push_back(num_triangles());
			edges_[c.node_uid_][indexcb].adj_triangles_.push_back(num_triangles());

			triangle_elem new_triangle;
			new_triangle.value = value;
			new_triangle.node_uid1 = a.node_uid_;
			new_triangle.node_uid2 = b.node_uid_;
			new_triangle.node_uid3 = c.node_uid_;
			new_triangle.triangle_uid = triangles_.size();
			triangles_.push_back(new_triangle);	  
		}
		return Triangle(this, triangles_.size()-1);
	}

	/** Return the triangle with index @a i.
	* @pre 0 <= @a i < num_triangles()
	* Complexity: O(1)
	*/
	Triangle triangle(size_type i) const {
		assert(i<num_triangles());
		return Triangle(this, i);
	}


	///////////////
	// Iterators //
	///////////////

	/** @class Mesh::NodeIterator
	* @brief Iterator class for nodes. A forward iterator. */
	class node_iterator : private totally_ordered<node_iterator>{
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
		node_iterator() {
		}
		/** Dereferencing operator
		* @return the node it points to
		* Complexity: O(1)
		*/ 
		Node operator*() const{
			return Node(m_, m_->node_i2u[index]);
		}

		/** Increment operator
		* @return the new iterator that points to the next node
		* Complexity: O(1)
		*/
		node_iterator& operator++(){
			++index;
			return *this;
		}

		/** Equal operator
		* @param[in] x a node iterator
		* @return true if they are equal
		* Complexity: O(1)
		*/
		bool operator==(const node_iterator& x) const{
			return (m_ == x.m_ && index == x.index);
		}

	private:
		friend class Mesh;
		node_iterator(const Mesh* mesh, size_type index_) 
			: m_(const_cast<Mesh*>(mesh)),index(index_){}
		Mesh* m_;
		size_type index;
	};

	/** Get the beginning node
	* @return the iterator that points to the first node of the graph
	* Complexity: O(1)
	*/
	node_iterator node_begin() const{
		return node_iterator(this, 0);
	}

	/** Get the last node
	* @return the iterator that points to the last node of the graph 
	*   (an invalid node)
	* Complexity: O(1)
	*/
	node_iterator node_end() const{
		return node_iterator(this,node_i2u.size());
	}

	/** @class Mesh::EdgeIterator
	* @brief Iterator class for edges. A forward iterator. */
	class edge_iterator : private totally_ordered<edge_iterator>{
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
		edge_iterator() {
		}

		/** Dereferencing operator
		* @return the edge it points to
		* Complexity: O(1)
		*/ 
		Edge operator*() const{
			return Edge(m_, uid1_, m_->edges_[uid1_][index].uid);
		}

		/** Increment operator
		* @return the new iterator that points to the next edge
		* Complexity: Best case O(1)
		*/
		edge_iterator& operator++(){
			do
			{
				if (++index >= m_->edges_[uid1_].size()){
					index = 0;
					if (++uid1_ == m_->edges_.size())
						break;
				}
			}
			while (uid1_ > m_->edges_[uid1_][index].uid);
			return *this;
		}

		/** Equal operator
		* @param[in] x an edge iterator
		* @return true if they are equal
		* Complexity: O(1)
		*/
		bool operator==(const edge_iterator& x) const{
			return (m_==x.m_ && uid1_ == x.uid1_ && index == x.index);
		}

	private:
		friend class Mesh;
		Mesh* m_;
		size_type uid1_;
		size_type index;
		edge_iterator(const Mesh* mesh, size_type uid1, size_type index_) 
			: m_(const_cast<Mesh*>(mesh)),uid1_(uid1), index(index_){}

	};

	/** Get the beginning edge
	* @return the iterator that points to the first edge of the graph
	* Complexity: O(1)
	*/
	edge_iterator edge_begin() const{
		return edge_iterator(this,0,0);
	}
   
	/** Get the last edge
	* @pre caller is a valid graph
	* @return the iterator that points to the last edge of 
	*   the graph (an invalid edge).
	* Complexity: O(1)
	*/
	edge_iterator edge_end() const{
		return edge_iterator(this,edges_.size(),0);
	}


	/**@class Mesh::IncidentIterator
	* @brief Iterator class for edges incident to a node. A forward iterator. */
	class IncidentIterator: private totally_ordered<IncidentIterator>{
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
		}

		/** Dereferencing operator
		* @return the edge it points to
		* Complexity: O(1)
		*/ 
		Edge operator*() const{
			return Edge(m_, node_uid1_, m_->edges_[node_uid1_][index].uid);
		}

		/** Increment operator
		* @return the new iterator that points to the next edge 
		*   in the adjacency list
		* Complexity: O(1)
		*/
		IncidentIterator& operator++(){
			++index;
			return *this;
		}

		/** Equal operator
		* @param[in] x an edge iterator
		* @return true if they are in the same graph and point to the same edge 
		* Complexity: O(1)
		*/
		bool operator==(const IncidentIterator& x) const{
			return (m_ == x.m_ && node_uid1_ == x.node_uid1_ && index == x.index); 
		}

	private:
		friend class Mesh;
		size_type node_uid1_;
		size_type index;
		Mesh* m_;
		IncidentIterator(const Mesh* mesh, size_type uid1, size_type index_) 
			: m_(const_cast<Mesh*>(mesh)),node_uid1_(uid1),index(index_){}
	};

	class triangle_iterator: private totally_ordered<triangle_iterator>{
	public:
		// These type definitions help us use STL's iterator_traits.
		/** Element type. */
		typedef Triangle value_type;
		/** Type of pointers to elements. */
		typedef Triangle* pointer;
		/** Type of references to elements. */
		typedef Triangle& reference;
		/** Iterator category. */
		typedef input_iterator_tag iterator_category;
		/** Difference between iterators */
		typedef ptrdiff_t difference_type;

		/** Construct an invalid triangleIterator. */
		triangle_iterator() {
		}

		/** Dereferencing operator
		* @return the triangle it points to
		* Complexity: O(1)
		*/ 
		Triangle operator*() const{
			return Triangle(m_,tri_uid_);
		}

		/** Increment operator
		* @return the new iterator that points to the next triangle
		* Complexity: O(1)
		*/
		triangle_iterator& operator++(){
			++tri_uid_;
			return *this;
		}
		/** Equal operator
		* @param[in] x a triangle iterator
		* @return true if they are equal
		* Complexity: O(1)
		*/
		bool operator==(const triangle_iterator& x) const{
			return (m_==x.m_ && tri_uid_==x.tri_uid_);
		}
		
		// all triangle functions from this point on implemented by us
    triangle_iterator& operator[](int n){
      assert(n < m_->triangles_.size() && n >= 0);
      tri_uid_ = n;
      return *this;
    }
    
    int operator-(triangle_iterator a){
      return tri_uid_ - a.tri_uid_;
    }

    triangle_iterator& operator+=(int n){
      tri_uid_ += n;
      return *this;
    }
    
    triangle_iterator& operator-=(int n){
      tri_uid_ -= n;
      return *this;
    }

	private:
		friend class Mesh;
		Mesh* m_;
		size_type tri_uid_;
		triangle_iterator(const Mesh* mesh, size_type uid): m_(const_cast<Mesh*>(mesh)),tri_uid_(uid){}
	};
	
  friend triangle_iterator operator+(triangle_iterator a, int n){
    return a+=n;
  }
  
  friend triangle_iterator operator+(int n, triangle_iterator a){
    return a+=n;
  }
  
  friend triangle_iterator operator-(triangle_iterator a, int n){
    return a-=n;
  }
	// end of our additions
	

	/* Return the first triangle interator
	* Complexity O(1)
	**/
	triangle_iterator triangle_begin(){
		return triangle_iterator(this,0);
	}


	/* Return the last triangle interator
	* Complexity O(1)
	**/
	triangle_iterator triangle_end(){
		return triangle_iterator(this,triangles_.size());
	}

	class nincident_tri_iterator: private totally_ordered<nincident_tri_iterator>{
	public:
		// These type definitions help us use STL's iterator_traits.
		/** Element type. */
		typedef Triangle value_type;
		/** Type of pointers to elements. */
		typedef Triangle* pointer;
		/** Type of references to elements. */
		typedef Triangle& reference;
		/** Iterator category. */
		typedef input_iterator_tag iterator_category;
		/** Difference between iterators */
		typedef ptrdiff_t difference_type;

		/** Construct an invalid nIncidentTriIterator. */
		nincident_tri_iterator() {}

		/** Dereferencing operator
		* @return the triangle it points to
		* Complexity: O(1)
		*/ 
		Triangle operator*() const{
			return Triangle(m_, m_->nodes_[uid_].adj_triangles_[adjtri_idx_]);
		}

		/** Increment operator
		* @return the triangle incident to node iterator that
		*    points to the next triangle incident to the given node
		* Complexity: O(1)
		*/
		nincident_tri_iterator& operator++(){
			++adjtri_idx_;
			return *this;
		}

		/** Equal operator
		* @param[in] x an triangle iterator
		* @return true if they are equal
		* Complexity: O(1)
		*/
		bool operator==(const nincident_tri_iterator& x) const{
			return (m_ == x.m_ && uid_==x.uid_ && adjtri_idx_ == x.adjtri_idx_);
		}

	private:
		friend class Mesh;
		Mesh* m_;
		size_type uid_;
		size_type adjtri_idx_;
		nincident_tri_iterator(const Mesh* mesh, size_type uid, size_type adjtri_idx)
			: m_(const_cast<Mesh*>(mesh)),uid_(uid), adjtri_idx_(adjtri_idx){}
	};

	class eincident_tri_iterator: private totally_ordered<eincident_tri_iterator>{
	public:
		// These type definitions help us use STL's iterator_traits.
		/** Element type. */
		typedef Triangle value_type;
		/** Type of pointers to elements. */
		typedef Triangle* pointer;
		/** Type of references to elements. */
		typedef Triangle& reference;
		/** Iterator category. */
		typedef input_iterator_tag iterator_category;
		/** Difference between iterators */
		typedef ptrdiff_t difference_type;

		/** Construct an invalid eincidentTriIterator. */
		eincident_tri_iterator() {}

		/** Dereferencing operator
		* @return the triangle it points to
		* Complexity: O(1)
		*/ 
		Triangle operator*() const{
			return Triangle(m_, m_->edges_[uid1_][uid2_index_].adj_triangles_[adjtri_idx_]);
		}

		/** Increment operator
		* @return the triangle incident to triangle iterator that points
		*     to the next triangle incident to the given triangle
		* Complexity: O(1)
		*/
		eincident_tri_iterator& operator++(){
			++adjtri_idx_;
			return *this;
		}

		/** Equal operator
		* @param[in] x an triangle iterator
		* @return true if they are equal
		* Complexity: O(1)
		*/
		bool operator==(const eincident_tri_iterator& x) const{
			return (m_ == x.m_ && uid1_==x.uid1_ && uid2_index_==x.uid2_index_ && adjtri_idx_ == x.adjtri_idx_);
		}

	private:
		friend class Mesh;
		Mesh* m_;
		size_type uid1_;
		size_type uid2_index_;
		size_type adjtri_idx_;
		eincident_tri_iterator(const Mesh* mesh, size_type uid1, size_type uid2_index, size_type adjtri_idx)
			: m_(const_cast<Mesh*>(mesh)),uid1_(uid1), uid2_index_(uid2_index),adjtri_idx_(adjtri_idx){}
	};

private:
	struct node_elem{
		Point points;
		size_type index;
		node_value_type value;
		vector<size_type> adj_triangles_;
	};

	struct edge_elem{
		size_type uid; 
		size_type edge_index;
		vector<size_type> adj_triangles_;
	};

	struct triangle_elem{
		// uid of 3 Nodes of this triangle
		size_type node_uid1; 
		size_type node_uid2;
		size_type node_uid3;
		// triangle id
		size_type triangle_uid;
		triangle_value_type value;
	};


	size_type edge_size_;
	vector<node_elem> nodes_;
	vector<vector<edge_elem>> edges_;
	vector<edge_value_type> edge_values;
	vector<triangle_elem> triangles_;
	vector<size_type> node_i2u;

};
