/*
 * Graph.h
 *
 *  Created on: Sep 30, 2010
 *      Author: andrew.gelfand
 */

#ifndef GRAPH_H_
#define GRAPH_H_

#include <algorithm>
#include <vector>
#include <iostream>
#include <list>
#include <map>
#include <queue>
#include <utility>
#include <cassert>
#include <limits>
#include <set>
#include <sstream>
#include <math.h>
#include "binaryheap.h"
#include "myRandom.h"
#include "Util.h"

using namespace std;

struct WeightedEdge {
	/**
	 * Two vertices comprising this edge object
	 */
	int m_u;
	int m_v;

	/**
	 * The weight of edge u-v
	 */
	int m_weight;

	WeightedEdge(int u, int v, int weight) {
		m_u = u;
		m_v = v;
		m_weight = weight;
	}
};

class CompareWeightedEdges {
private:
	bool m_isIncreasing;
public:
	CompareWeightedEdges(const bool& isIncreasing = true) {
		m_isIncreasing = isIncreasing;
	};

	int operator() (WeightedEdge& a, WeightedEdge& b) {
		if (m_isIncreasing) {
			return a.m_weight > b.m_weight;
		}
		else {
			return a.m_weight < b.m_weight;
		}
	}
};

struct WeightedVertex {
protected:
	/**
	 * Vertex id
	 */
	int m_vertex;

	/**
	 * Vertex weight
	 */
	int m_weight;

	uint32_t m_heap_pos;
public:

	WeightedVertex(int v, int weight) {
		m_vertex = v;
		m_weight = weight;
		m_heap_pos = 0;
	}

	int getVertex() {
		return m_vertex;
	}

	void setWeight(int weight) {
		m_weight = weight;
	}

	int getWeight() {
		return m_weight;
	}

	bool operator<(WeightedVertex& other) {
	  return (m_weight > other.getWeight());
	}

	uint32_t& heap_pos() { return m_heap_pos; }
};

class CompareWeightedVertices {
public:
	int operator() (WeightedVertex* a, WeightedVertex* b) {
		return (a->getWeight() > b->getWeight());
	}
};

struct UpdateWeightedVertexPos {
  void operator()(WeightedVertex *a, uint32_t pos) {
    a->heap_pos() = pos;
  }
};

struct Graph {
protected:

	/**
	 * Indicates if the graph is weighted or not - i.e. has weighted edges
	 */
	bool m_isWeighted;

    /**
     * An adjacency list backing this graph object.  The adjacency lists
     * representation is <i>n</i> lists such that the <i>i</i>-th list contains
     * a sequence of neighbors of vertex <i>i</i>.
     */
    vector<vector<int> > m_adjList;

    /**
     * The number of vertices
     */
    int m_numVertices;

    /**
     * The number of edges
     */
    int m_numEdges;

    /**
     * The weights of the edges in the graph. Each edge i-j in the graph is associated
     * with the following unique key:
     * 		key = n*(i) + j (if i < j) and key = n*(j) + i (if i > j)
     * where n is the number of vertices in the graph.
     */
    map<int,int> m_edgeIDToWeightMap;

	/**
	 * Computes the cost of removing the specified vertex, where cost is
	 * specified in terms of the number of fill-in edges added
	 */
	inline int computeMinFillRemovalCost(int v) {
		// Get the adjacency list of this vertex
		vector<int>::iterator cur, last, inner;
		cur = m_adjList[v].begin();
		last = m_adjList[v].end();

		int edgesAdded = 0;

		while (cur != last) {
			inner = cur;
			inner++;
			while (inner != last) {
				if (!containsEdge(*cur, *inner)) {
					edgesAdded++;
				}
				inner++;
			}
			cur++;
		}

		return edgesAdded;
	}

    /**
	 * Removes vertex <i>i</i> from the graph and also removes all edges
	 * connecting to/from this vertex.
	 * @param i
	 */
	inline void removeVertex(int i) {
		// Remove arcs involving the variable to be deleted
		vector<int>::iterator current, last, it;
		current = m_adjList[i].begin();
		last = m_adjList[i].end();
		while (current != last) {
			vector<int> newList;
			for (it = m_adjList[*current].begin(); it != m_adjList[*current].end(); it++) {
				if ((*it) != i) {
					newList.push_back(*it);
				}
			}
			m_adjList[*current].clear();
			m_adjList[*current] = newList;
			current++;
		}
		m_adjList[i].clear();
	}

	/**
	 * Recursively removes any dangling edges from this graph.  A dangling edge is an
	 * edge containing a vertex with degree 1.
	 */
	void removeDanglingEdges() {
		// The set of vertices with degree 1
		vector<int> degreeOneVertices;
		for (int i = 0; i < m_numVertices; i++) {
			if (getDegreeOf(i) == 1) {
				degreeOneVertices.push_back(i);
			}
		}

		while (!degreeOneVertices.empty()) {
			int v = degreeOneVertices.back();
			degreeOneVertices.pop_back();

			// Remove the edge this vertex is involved with
			if (getDegreeOf(v) == 1) {
				// Get the neighboring vertex
				int u = getAdjacentVertices(v)[0];

				// Remove the edge between v and u
				removeEdge(v,u);

				// Add u to list if it is singly connected
				if (getDegreeOf(u) == 1) {
					degreeOneVertices.push_back(u);
				}
			}
		}
	}

	/**
	 * Compares two items that are pairs of integers.  Comparison is based upon the
	 * second item in the pair (i.e. a.second) and returns true if (a.second > b.second).
	 * If used by a sort operation or priority queue, this operator will place items in
	 * ascending order by size.
	 */
	class CompareDistances {
	public:
		int inline operator() (const pair<int,int>& a, const pair<int,int>& b) const {
			return (a.second > b.second);
		}
	};

	/**
	 * Recursively searches for a path of length >= path length that forms a cycle.  The
	 * search begins from vertex v.
	 */
	inline bool searchPath(int v, int pathLength, vector<bool>& visited, vector<int>& path) {
		visited[v] = true;
		path.push_back(v);

		if (((int)path.size()) == pathLength) {
			if (containsEdge(path.front(),path.back())) {
				// Indicate that a cycle has been found
				return true;
			}

			path.pop_back();
			visited[v] = false;
			return false;
		}

		// Search all neighbors that have not been visited
		for (int i = 0; i < (int)m_adjList[v].size(); i++) {
			int u = m_adjList[v][i];
			if (!visited[u]) {
				bool foundCycle = searchPath(u, pathLength, visited, path);

				if (foundCycle) {
					return foundCycle;
				}
			}
		}

		visited[v] = false;
		path.pop_back();

		return false;
	}

	inline void recursiveDFS(int root, int v, int cutoffDepth, vector<bool>& visited, vector<int>& path, set<int>& vertices) {
		visited[v] = true;
		path.push_back(v);

		if (v >= root) {
			vertices.insert(v);
		}
		if (((int)path.size()) == cutoffDepth) {
			path.pop_back();
			visited[v] = false;
			return;
		}

		// Search all neighbors that have not been visited
		for (int i = 0; i < (int)m_adjList[v].size(); i++) {
			int u = m_adjList[v][i];
			if (!visited[u]) {
				recursiveDFS(root, u, cutoffDepth, visited, path, vertices);
			}
		}
		visited[v] = false;
		path.pop_back();
		return;
	}

	/**
	 * Recursively searches for a path of length >= path length. If <tt>isCyclesOnly</tt> is true, this path
	 * must form a cycle.  The search begins from vertex v. Only paths over unique sets of vertices are
	 * retained in <tt>cycleVertices</tt>.
	 * @param cycleVertices - Contains the set of vertices in each cycle.  The order of the path is not preserved
	 */
	inline void searchPath(int v, int pathLength, vector<bool>& visited, vector<int>& path, vector<vector<int> >& cycleVertices, bool isCyclesOnly) {
		visited[v] = true;
		path.push_back(v);
		if (((int)path.size()) == pathLength) {
			if ((!isCyclesOnly) || containsEdge(path.front(),path.back())) {

				// Add this cycle to the the current set of cycles
				vector<int> temp = path;
				sort(temp.begin(), temp.end());

				bool isRedundant = false;
				for (int i = 0; i < (int)cycleVertices.size(); i++) {
					if (temp.size() == cycleVertices[i].size()) {
						isRedundant = true;
						for (int j = 0; j < (int)cycleVertices[i].size(); j++) {
							if (temp[j] != cycleVertices[i][j]) {
								isRedundant = false;
								break;
							}
						}
					}
					if (isRedundant) {
						break;
					}
				}

				// Add the cycle if it is not redundant
				if (!isRedundant) {
					cycleVertices.push_back(temp);
				}
			}
			path.pop_back();
			visited[v] = false;
			return;
		}

		// Search all neighbors that have not been visited
		for (int i = 0; i < (int)m_adjList[v].size(); i++) {
			int u = m_adjList[v][i];
			if (!visited[u]) {
				searchPath(u, pathLength, visited, path, cycleVertices, isCyclesOnly);
			}
		}
		visited[v] = false;
		path.pop_back();

		return;
	}

	vector<vector<int> > m_maximalCliques;

	/**
	 * The bron-kerbosch algorithm
	 */
	void BronKerbosch(vector<int> R, vector<int> P, vector<int> X) {
		if (P.empty() && X.empty()) {
			m_maximalCliques.push_back(R);
		}
		else {
			while (!P.empty()) {
				int v = P.front();
				vector<int> Rnew;
				addInt(R,v,Rnew);
				vector<int> Pnew;
				getIntersectionInt(P, m_adjList[v], Pnew);
				vector<int> Xnew;
				getIntersectionInt(X, m_adjList[v], Xnew);
				BronKerbosch(Rnew, Pnew, Xnew);

				removeInt(P,v,P);
				addInt(X,v,X);
			}
		}
	}

	/**
	 * Variables used for checking two edge connectivity
	 */
	map<int,int> m_pre;
	map<int,int> m_low;
	int m_cnt;
	vector<bool> m_included;
	vector<pair<int,int> > m_bridgeEdges;
	void dfs_bridge(int u, int v) {
		m_pre[v] = m_cnt++;
        m_low[v] = m_pre[v];
        for (int i = 0; i < (int)m_adjList[v].size(); i++) {
        	int w = m_adjList[v][i];
            if (m_included[w] && (m_pre[w] == -1)) {
                dfs_bridge(v, w);
                m_low[v] = min(m_low[v], m_low[w]);
                if (m_low[w] == m_pre[w]) {
                    pair<int,int> edge(v,w);
                    m_bridgeEdges.push_back(edge);
                }
            }
            else if (m_included[w] && (w != u)) {
                m_low[v] = min(m_low[v], m_pre[w]);
        	}
    	}
	}

	struct Node {
	protected:
		Node* m_parent;
		int m_vertex;

	public:
		Node(Node* parent, int vertex) {
			m_parent = parent;
			m_vertex = vertex;
		}

		Node* getParent() {
			return m_parent;
		}

		int getVertex() {
			return m_vertex;
		}
	};

public:

    /**
     * Constructs a new Graph object for the specified number of vertices. The graph
     * is initialized to contain no edges.
     * @param numVertices -
     */
    Graph(int numVertices) {
        m_numVertices = numVertices;
        m_adjList = vector<vector<int> > (m_numVertices);
        m_numEdges = 0;
        m_isWeighted = false;
        m_cnt = 0;
    }

    /**
	 * Constructs a new weighted Graph object for the specified number of vertices. The graph
	 * is initialized to contain no edges.
	 * @param numVertices -
	 * @param isWeighed  -
	 */
	Graph(int numVertices, bool isWeighted) {
		m_numVertices = numVertices;
		m_adjList = vector<vector<int> > (m_numVertices);
		m_numEdges = 0;
		m_isWeighted = isWeighted;
		m_cnt = 0;
	}

    /**
     * Copy constructor.
     * @param Graph graph -
     */
    Graph(Graph& graph) {
        m_numVertices = graph.m_numVertices;

        // Initialize the adjacency lists and adjacency matrix
        m_adjList = vector<vector<int> > (m_numVertices);

        // And copy from graph object to the new object
        for (int i = 0; i < m_numVertices; i++) {
        	vector<int>::iterator current, last;
        	current = graph.m_adjList[i].begin();
        	last = graph.m_adjList[i].end();
        	while (current != last) {
        		m_adjList[i].push_back(*current);
        		current++;
        	}
        }
        m_numEdges = graph.m_numEdges;

        m_isWeighted = graph.m_isWeighted;
        m_edgeIDToWeightMap = graph.m_edgeIDToWeightMap;
        m_cnt = 0;
    }

    /**
     * Merges this graph and the specified graph (i.e. takes the union of the
     * edges in the two graphs).
     */
    void merge(Graph* graph) {
    	assert(m_numVertices == graph->m_numVertices);
    	assert(!m_isWeighted);

    	// Perform the union of the edges
		for (int i = 0; i < m_numVertices; i++) {
			vector<int>::iterator current, last;
			current = graph->m_adjList[i].begin();
			last = graph->m_adjList[i].end();
			while (current != last) {
				addEdge(i, *current);
				current++;
			}
		}
    }

    /**
     * Generates a string output of the current graph
     */
    string toString() {
    	std::stringstream ss;
    	for (int i = 0; i < m_numVertices; i++) {
    		if (!m_adjList[i].empty()) {
				vector<int>::iterator current, last;
				current = m_adjList[i].begin();
				last = m_adjList[i].end();
				ss << i << ":";
				while (current != last) {
					ss << " " << *current;
					current++;
				}
				ss << endl;
    		}
    	}
    	return ss.str();
    }

    /**
     * Adds an edge to the graph.
     *
     * NOTE:  Checking is done to prevent addition of duplicate edges in the
     * adjacency lists
     * @param i
     * @param j
     */
    void addEdge(int i, int j) {
    	assert(i >= 0); assert(j >= 0);
    	assert(i < m_numVertices); assert(j < m_numVertices);
    	if (i == j) {
    		return;
    	}
    	assert(!m_isWeighted);

    	// If this edge does not exist...add it
    	if (!containsEdge(i,j)) {
    		m_adjList[i].push_back(j);
    		std::sort(m_adjList[i].begin(),m_adjList[i].end());
			m_adjList[j].push_back(i);
			sort(m_adjList[j].begin(),m_adjList[j].end());

			m_numEdges++;
    	}
    }

    /**
     * Adds an edge to the graph.
     *
     * NOTE:  Checking is done to prevent addition of duplicate edges in the
     * adjacency lists
     * @param i
     * @param j
     * @param weight - the weight of this edge
     */
    void addEdge(int i, int j, int weight) {
    	assert(i >= 0); assert(j >= 0);
		assert(i < m_numVertices); assert(j < m_numVertices);
    	if (i == j) {
    		return;
    	}
    	assert(m_isWeighted);

    	// If this edge does not exist...add it
    	if (!containsEdge(i,j)) {
    		m_adjList[i].push_back(j);
    		std::sort(m_adjList[i].begin(),m_adjList[i].end());
			m_adjList[j].push_back(i);
			sort(m_adjList[j].begin(),m_adjList[j].end());

			if (i < j) {
				int key = m_numVertices*i + j;
				m_edgeIDToWeightMap[key] = weight;
			}
			else {
				int key = m_numVertices*j + i;
				m_edgeIDToWeightMap[key] = weight;
			}

			m_numEdges++;
    	}
    }

    /**
     * Returns <tt>true</tt> if the edge from i to j exists in this graph and
     * <tt>false</tt> otherwise. If <tt>true</tt> updates the weight of the
     * specified edge too.
     */
    inline bool containsEdge(int i, int j, int& weight) {
    	if (binary_search(m_adjList[i].begin(), m_adjList[i].end(), j)) {
    		if (i < j) {
    			weight = m_edgeIDToWeightMap[(m_numVertices*i + j)];
    		}
    		else {
    			weight = m_edgeIDToWeightMap[(m_numVertices*j + i)];
    		}
    		return true;
    	}
    	return false;
    }

    /**
     * Returns <tt>true</tt> if the edge from i to j exists in this graph and
     * <tt>false</tt> otherwise. Since the graph is undirected, the arc from i to j is
     * equivalent to the arc from j to i.
     * @param int i -
     * @param int j -
     */
    inline bool containsEdge(int i, int j) {
    	return binary_search(m_adjList[i].begin(), m_adjList[i].end(), j);
    }

    /**
     * Removes the specified edge if it exists
     */
    void removeEdge(int i, int j) {
    	if (i == j) {
    		return;
    	}
    	assert(!m_isWeighted);

    	if (containsEdge(i,j)) {
    		// Remove vertex i from vertex j's adjacency list
    		int idx = -1;
    		for (int k = 0; k < (int)m_adjList[i].size(); k++) {
    			if (m_adjList[i][k] == j) {
    				idx = k;
    				break;
    			}
    		}
    		assert(idx != -1);
    		m_adjList[i].erase(m_adjList[i].begin() + idx);

    		// Remove vertex j from vertex i's adjacency list
			idx = -1;
			for (int k = 0; k < (int)m_adjList[j].size(); k++) {
				if (m_adjList[j][k] == i) {
					idx = k;
					break;
				}
			}
			assert(idx != -1);
			m_adjList[j].erase(m_adjList[j].begin() + idx);

    		m_numEdges--;
    	}
    }

    /**
     * Get the set of vertices with degree > 0
     */
    vector<int> getActiveVertices() {
    	vector<int> returnSet;
    	for (int i = 0; i < m_numVertices; i++) {
    		if (!m_adjList[i].empty()) {
    			returnSet.push_back(i);
    		}
    	}
    	return returnSet;
    }

    /**
     * Gets the degree of vertex i in the graph.
     * @param i
     */
    inline int getDegreeOf(int i) {
        int degree = m_adjList[i].size();
        return degree;
    }

    /**
     * Gets the connectivity of vertex i, where
     * connectivity is defined as the number of
     * vertices neighboring vertex i that are
     * also adjacent.
     */
    inline int getConnectivityOf(int i) {
    	int connectivity = 0;
    	for (int j = 0; j < (int)m_adjList[i].size(); j++) {
    		for (int k = j + 1; k < (int)m_adjList[i].size(); k++) {
    			if (containsEdge(m_adjList[i][j],m_adjList[i][k])) {
    				connectivity++;
    			}
    		}
    	}
    	return connectivity;
    }

    /**
     * Gets the number of edges in the graph
     */
    int getNumberOfEdges() {
    	return m_numEdges;
    }

    /**
     * Gets the number of vertices in the graph
     */
    int getNumberOfVertices() {
    	return m_numVertices;
    }

    /**
     * Gets the set of vertices adjacent to i
     */
    inline const vector<int>& getAdjacentVertices(int i) {
    	return m_adjList[i];
    }

    /**
     * Determines if the specified set of vertices form an edge-connected
     * subgraph. Returns <tt>true</tt> if the graph is edge-connected;
     * <tt>false</tt> otherwise.
     */
    bool isEdgeConnected(vector<int>& vertices) {
    	if (((int)vertices.size()) < 3) {
    		return false;
    	}

    	m_low.clear();
    	m_pre.clear();
    	m_included = vector<bool>(m_numVertices);
    	fill(m_included.begin(), m_included.end(), false);
    	for (int i = 0; i < (int)vertices.size(); i++) {
    		assert((vertices[i] >= 0) && (vertices[i] < m_numVertices));
    		m_included[vertices[i]] = true;
    		m_low[vertices[i]] = -1;
    		m_pre[vertices[i]] = -1;
    	}

    	m_cnt = 0;
    	m_bridgeEdges.clear();
    	for (int i = 0; i < (int)vertices.size(); i++) {
    		if (m_pre[vertices[i]] == -1) {
    			dfs_bridge(vertices[i],vertices[i]);
    		}
    	}

    	return (m_bridgeEdges.empty());
    }

    /**
	 * Gets the shortest path from vertex s to any terminal
	 * vertex.  <tt>terminal</tt> is assumed to be a vector
	 * of length <tt>getNumberOfVertices()</tt>
	 * where terminal[i] is <tt>true</tt> if vertex i is a
	 * terminal vertex and <tt>false</tt> otherwise. Does not
	 * permit cycles.
	 * @return vector<int> -- Returns an empty vector if no such path can be found.
	 */
	vector<int> getShortestPath(int s, vector<bool> terminal) {
		assert(((int)terminal.size()) == m_numVertices);
		assert(s >= 0);
		assert(s < m_numVertices);

		// Create the queue
		queue<Node*> fifo;

		// Create a data structure for clean-up too
		vector<Node*> nodes;

		bool foundTerminalNode = false;

		// Do a breadth-first search from vertex s, stopping at the first
		// terminal vertex encountered
		Node* root = new Node(NULL, s);
		fifo.push(root);
		nodes.push_back(root);

		Node* next;
		while (!fifo.empty()) {

			next = fifo.front();
			fifo.pop();

			int v = next->getVertex();

			// Goal check
			if (terminal[v]) {
				foundTerminalNode = true;
				break;
			}

			// Generate children of vertex v
			for (int i = 0; i < (int)m_adjList[v].size(); i++) {
				int u = m_adjList[v][i];

				// Don't consider obvious loops
				if ((u == s) || ((next->getParent() != NULL) && (next->getParent()->getVertex() == u))) {
					continue;
				}

				Node* node = new Node(next, u);
				fifo.push(node);
				nodes.push_back(node);
			}
		}

		// Construct the path
		vector<int> path;
    	if (foundTerminalNode) {
    		path.push_back(next->getVertex());
    		next = next->getParent();
    		while (next != NULL) {
    			path.push_back(next->getVertex());
    			next = next->getParent();
    		}
    	}

    	// Error check: Path must not contain a cycle
    	for (int i = 0; i < (int)path.size(); i++) {
    		for (int j = i + 1; j < (int)path.size(); j++) {
    			assert(path[i] != path[j]);
    		}
    	}

    	// Clean up
    	for (int i = 0; i < (int)nodes.size(); i++) {
    		delete nodes[i];
    	}
    	nodes.clear();

    	return path;
	}

    /**
     * Gets a cycle of smallest length from vertex v.
     */
    vector<int> getShortestCycleFromV(int v) {
    	assert(v >= 0);
    	assert(v < m_numVertices);

    	bool foundCycle = false;
    	vector<int> path;

    	if (m_numVertices == 3) {
    		if (containsEdge(0,1) && containsEdge(0,2) && containsEdge(1,2)) {
    			path.push_back(0); path.push_back(1); path.push_back(2);
    			return path;
    		}
    	}

    	// Start search with paths of length 3
    	for (int i = 3; i < m_numVertices; i++) {
    		vector<bool> visited = vector<bool>(m_numVertices);
    		fill(visited.begin(),visited.end(),false);
    		path.clear();

    		foundCycle = searchPath(v,i,visited,path);

    		if (foundCycle) {
    			break;
    		}
    	}

    	return path;
    }

    /**
     * Enumerates all paths of the specified length originating from vertex <tt>v</tt>.
     */
    void getPathsFromV(int v, int pathLength, vector<vector<int> >& vertices) {
    	assert((v >= 0) && (v < m_numVertices));
    	vector<bool> visited = vector<bool>(m_numVertices);
		fill(visited.begin(), visited.end(), false);
		vector<int> path;
		searchPath(v, pathLength, visited, path, vertices, false);
    }

    /**
     * Gets all vertices u that are within k number of edges of v (excluding v).
     */
    vector<int> getkNeighborhood(int v, int k) {
    	assert((v >= 0) && (v < m_numVertices));

    	set<int> vertices;
    	vector<int> path;
		vector<bool> visited = vector<bool>(m_numVertices);
		fill(visited.begin(), visited.end(), false);

    	recursiveDFS(v, v, k+1, visited, path, vertices);

    	vector<int> neighborhood;
    	set<int>::iterator it;
    	for (it = vertices.begin(); it != vertices.end(); it++) {
    		if (*it != v) {
    			neighborhood.push_back((*it));
    		}
    	}
    	return neighborhood;
    }

    /**
     * Enumerates all simple cycles of the specified length originating from <tt>v</tt>. Only
     * the unique sets of vertices are retained in <tt>cycleVertices</tt>. So, if two different
     * cycles traverse the same set of vertices, cycleVertices will only contain one of the vertex sets.
     */
    void getCyclesFromV(int v, int cycleLength, vector<vector<int> >& cycleVertices) {
    	if ((cycleLength < 3) || (cycleLength > m_numVertices)) {
    		return;
    	}
    	assert((v >= 0) && (v < m_numVertices));

    	vector<bool> visited = vector<bool>(m_numVertices);
    	fill(visited.begin(), visited.end(), false);
    	vector<int> path;
    	searchPath(v, cycleLength, visited, path, cycleVertices, true);
    }

    /**
     * Gets the number of connected components in this graph.
     */
    int getNumberOfConnectedComponents() {
    	vector<bool> visited = vector<bool>(m_numVertices);
    	fill(visited.begin(),visited.end(),false);
    	int numVisited = 0;

    	int numComponents = 0;

    	// Must visit all vertices
    	while (numVisited != m_numVertices) {

    		// Add the first non-visited vertex to the stack
			vector<int> mystack;
			int idx = -1;
			for (int i = 0; i < m_numVertices; i++) {
				if (!visited[i]) {
					idx = i;
					break;
				}
			}
			assert(idx >= 0);
			mystack.push_back(idx);

			numComponents++;

			// Perform a depth first search
			while (!mystack.empty()) {

				idx = mystack.back();
				mystack.pop_back();
				if (visited[idx]) {
					continue;
				}
				numVisited++;
				visited[idx] = true;

				// Push unvisited neighbors of idx onto the stack
				for (int j = 0; j < (int)getAdjacentVertices(idx).size(); j++) {
					if (!visited[getAdjacentVertices(idx)[j]]) {
						mystack.push_back(getAdjacentVertices(idx)[j]);
					}
				}
			}
    	}

    	return numComponents;
    }

    /**
     * Gets the vertices in each connected component of this graph
     */
    vector<vector<int> > getConnectedComponents() {
    	vector<bool> visited = vector<bool>(m_numVertices);
    	fill(visited.begin(),visited.end(),false);
    	int numVisited = 0;

    	vector<vector<int> > components;

    	// Must visit all vertices
    	while (numVisited != m_numVertices) {

    		vector<int> curComponent;

    		// Add the first non-visited vertex to the stack
			vector<int> mystack;
			int idx = -1;
			for (int i = 0; i < m_numVertices; i++) {
				if (!visited[i]) {
					idx = i;
					break;
				}
			}
			assert(idx >= 0);
			mystack.push_back(idx);

			// Perform a depth first search
			while (!mystack.empty()) {

				idx = mystack.back();
				mystack.pop_back();
				if (visited[idx]) {
					continue;
				}
				numVisited++;
				visited[idx] = true;
				curComponent.push_back(idx);

				// Push unvisited neighbors of idx onto the stack
				for (int j = 0; j < (int)getAdjacentVertices(idx).size(); j++) {
					if (!visited[getAdjacentVertices(idx)[j]]) {
						mystack.push_back(getAdjacentVertices(idx)[j]);
					}
				}
			}

			sort(curComponent.begin(), curComponent.end());
			components.push_back(curComponent);
    	}

    	return components;
    }


    /**
     * Returns the cyclomatic number, mu, which is equal to mu=|E|-|V|+c, where |E|
     * is the number of edges in the graph, |V| is the number of vertices in the graph
     * and c is the number of connected components.
     */
    int getCycleSpaceDimension() {
    	// Determine the number of cycles in this graph
		int c = getNumberOfConnectedComponents();
		int mu = m_numEdges - m_numVertices + c;
		return mu;
    }

    /**
     * Gets a fundamental cycle basis of this graph object.  The number of cycles
     * is equal to mu=|E|-|V|+c, where |E| is the number of edges in the graph, |V| is
     * the number of vertices in the graph and c is the number of connected components.
     * Returns a vector cycles of length mu, where cycles[i] is the ith cycle in the basis.
     * If all the connected components of the graph are not 2-edge connected, then a set of
     * 'bridge' edges will be returned in <tt>edges</tt>.  For example, imagine we have the
     * following barbell graph:
     *   0     3
     *   |>1-2<|
     *   4     5
     * The cyclomatic number of this graph is 7 - 6 + 1 = 2, so cycles would contain
     * cycles[0]=<0,1,4> and cycles[1]=<2,3,5>. In addition, since edge 1-2 is a bridge edge,
     * edges would contain edges[0]=<1,2>.
     *
     */
    void getFundamentalCycleBasis(vector<vector<int> >& cycles, vector<vector<int> >& edges, myRandom* random) {
    	// Determine the number of cycles in this graph
    	int c = getNumberOfConnectedComponents();
    	int mu = m_numEdges - m_numVertices + c;
    	//cout << "mu=" << mu << "=" << m_numEdges << "-" << m_numVertices << "+" << c << endl;

    	// Clear the return object
		cycles.clear();
		edges.clear();

    	if (mu == 0) {
			return;
    	}

    	// Populate the cycle basis vector
    	cycles = vector<vector<int> >(mu);

    	// Keep track of the set of unvisited edges
    	Graph* unvisitedEdgeGraph = new Graph(*this);

    	// Remove any dangling edges
    	unvisitedEdgeGraph->removeDanglingEdges();

    	// Keep track of the visited vertices
    	vector<bool> visited = vector<bool>(m_numVertices);
    	for (int i = 0; i < m_numVertices; i++) {
    		if (unvisitedEdgeGraph->getDegreeOf(i) == 0) {
    			visited[i] = true;
    		}
    		else {
    			visited[i] = false;
    		}
    	}

    	// Form an initial cycle from some unvisited vertex
		vector<int> verts;
		for (int i = 0; i < m_numVertices; i++) {
			if (!visited[i]) {
				verts.push_back(i);
			}
		}
		int v = verts[random->getInt((int)verts.size())];
		assert(v >= 0);

		// Get a neighbor of vertex v
		int u = unvisitedEdgeGraph->getAdjacentVertices(v)[0];
		unvisitedEdgeGraph->removeEdge(v, u);
		visited[v] = true;
		cycles[0] = unvisitedEdgeGraph->getShortestPath(u,visited);
		assert(!cycles[0].empty());

		// Keep track of the vertices visited in this component
		set<int> visitedInCurComponent;
		for (int i = 0; i < (int)cycles[0].size(); i++) {
			visited[cycles[0][i]] = true;
			visitedInCurComponent.insert(cycles[0][i]);
			if (i > 0) {
				unvisitedEdgeGraph->removeEdge(cycles[0][i-1], cycles[0][i]);
			}
		}

		// Mark this is a cycle rather than ear
		cycles[0].push_back(cycles[0].front());

    	// Now begin the search
    	int numCycles = 1;
    	int componentNum = 1;
    	bool usePrevVert = false;
    	while (numCycles < mu) {

    		// If the previous ear was a 'bridge' edge. Follow the bridge until we find another cycle.
    		// Note that all dangling edges have been removed.
    		if (usePrevVert) {
    			v = u;
    			assert (unvisitedEdgeGraph->getDegreeOf(v) > 0);
				u = unvisitedEdgeGraph->getAdjacentVertices(v)[0];
				usePrevVert = false;
    		}
    		else {
				// While visited vertices in this component have unvisited edges
				u = -1;
				set<int>::iterator it;
				for (it = visitedInCurComponent.begin(); it != visitedInCurComponent.end(); it++) {
					v = (*it);

					if (unvisitedEdgeGraph->getDegreeOf(v) > 0) {
						u = unvisitedEdgeGraph->getAdjacentVertices(v)[0];
						break;
					}
				}
    		}

    		// If there is an edge to start an ear...
    		if (u >= 0) {

    			// Remove the edge v and u
				unvisitedEdgeGraph->removeEdge(v, u);

				// If both v and u have been visited...the ear is a chord
				if (visited[u]) {
					cycles[numCycles].push_back(u);
					cycles[numCycles].push_back(v);
				}
				else {
					// Otherwise, try to get a path from vertex u to a previously visited vertex
					vector<int> path = unvisitedEdgeGraph->getShortestPath(u,visited);

					// If the returned path is empty, then the graph is not two-edge connected
					// so this edge is a bridge edge
					if (path.empty()) {
						vector<int> bridgeEdge;
						bridgeEdge.push_back(v); bridgeEdge.push_back(u);
						edges.push_back(bridgeEdge);
						visited[u] = true;
						usePrevVert = true;
						numCycles--;
					}
					else {
						cycles[numCycles] = path;
						cycles[numCycles].push_back(v);
						for (int j = 0; j < (int)path.size(); j++) {

							// Mark each vertex along this path as visited and remove all the
							// edges on this path
							if (j > 0) {
								unvisitedEdgeGraph->removeEdge(path[j-1],path[j]);
							}
							visited[path[j]] = true;
							visitedInCurComponent.insert(path[j]);
						}
					}
				}
    		}
    		// Find a new cycle
    		else {
    			componentNum++;

    			verts.clear();
				for (int i = 0; i < m_numVertices; i++) {
					if (!visited[i]) {
						verts.push_back(i);
					}
				}
				v = verts[random->getInt((int)verts.size())];
				assert(v >= 0);

				u = unvisitedEdgeGraph->getAdjacentVertices(v)[0];
				unvisitedEdgeGraph->removeEdge(v, u);
				visited[v] = true;
				cycles[numCycles] = unvisitedEdgeGraph->getShortestPath(u,visited);
				assert(!cycles[numCycles].empty());

				// Keep track of the vertices visited in this component
				visitedInCurComponent.clear();
				for (int i = 0; i < (int)cycles[numCycles].size(); i++) {
					visited[cycles[numCycles][i]] = true;
					visitedInCurComponent.insert(cycles[numCycles][i]);
					if (i > 0) {
						unvisitedEdgeGraph->removeEdge(cycles[numCycles][i-1], cycles[numCycles][i]);
					}
				}

				// Mark this as a cycle
				cycles[numCycles].push_back(cycles[numCycles].front());
    		}
    		numCycles++;
    	}

    	// Sanity Check
    	assert(unvisitedEdgeGraph->getNumberOfEdges() == 0);

    	// Now complete cycles for the ears we have found
    	for (int i = 0; i < mu; i++) {
    		// If this is a cycle...
    		if (cycles[i].front() == cycles[i].back()) {

    			// Mark all edges on this cycle as visited
    			for (int j = 1; j < (int)cycles[i].size(); j++) {
    				unvisitedEdgeGraph->addEdge(cycles[i][j-1],cycles[i][j]);
    			}
    			cycles[i].pop_back();
    		}
    		// Otherwise...this is an ear
    		else {
    			// Get the start and end point of this ear path
				int pathLength = cycles[i].size();
				assert(pathLength >= 2); // Path must be at least an edge
				int u = cycles[i][0];
				int v = cycles[i][pathLength-1];

				// Get the shortest path from u to v on the set of edges already
				// in the graph
				vector<bool> terminal = vector<bool>(m_numVertices);
				fill(terminal.begin(),terminal.end(),false);
				terminal[v] = true;
				vector<int> path = unvisitedEdgeGraph->getShortestPath(u,terminal);

				// Add the ear edges to the visitedGraph
				for (int j = 1; j < (int)cycles[i].size(); j++) {
					unvisitedEdgeGraph->addEdge(cycles[i][j-1],cycles[i][j]);
				}

				// Create the cycle. Note that the first and last element of the path are
				// v and u and are therefore redundant - i.e. already appear in ears[i]
				for (int j = 1; j < (int)path.size() - 1; j++) {
					cycles[i].push_back(path[j]);
				}
    		}
    	}

    	delete unvisitedEdgeGraph;

    	// Quick validation
    	for (int i = 0; i < mu; i++) {
    		int cycleLength = cycles[i].size();
    		for (int j = 1; j < cycleLength; j++) {
    			assert(containsEdge(cycles[i][j-1],cycles[i][j]));
    		}
    		assert(containsEdge(cycles[i][0],cycles[i][cycleLength-1]));
    	}
    }

    /**
     * Gets a min-fill ordering of the vertices in this graph with degree > 0.  Assumes that the
     * graph of such vertices is connected. Populates <tt>order</tt> with the discovered elimination order.
     */
    static vector<int> getMinFillOrdering(Graph* inGraph, vector<int> evidence, myRandom* random, vector<vector<int> >& clusters, int pool) {
    	assert(pool >= 0);

    	// Get a copy of the graph object
    	Graph* graph = new Graph(*inGraph);

    	vector<int> order;
    	clusters.clear();

    	vector<bool> processed = vector<bool>(graph->m_numVertices);
    	fill(processed.begin(),processed.end(),false);
    	for (int i = 0; i < (int)evidence.size(); i++) {
    		assert((evidence[i] >= 0) && (evidence[i] < graph->m_numVertices));
    		processed[evidence[i]] = true;
    	}

    	vector<WeightedVertex*> vertices = vector<WeightedVertex*>(graph->m_numVertices);
    	vector<WeightedVertex*> verticesHeap;

    	// 1) Mark all unconnected vertices as processed;
    	// 2) Pre-compute the removal costs of connected vertices;
    	// 3) Identify any simplical vertices
        int numInOrdering = 0;
        for (int i = 0; i < graph->m_numVertices; i++) {
        	if (graph->getDegreeOf(i) == 0) {
        		processed[i] = true;
        	}
        	else if ((graph->getDegreeOf(i) <= 1) && (!processed[i]))  {

        		vertices[i] = new WeightedVertex(i,0);
        		verticesHeap.push_back(vertices[i]);
        		numInOrdering++;
        	}
        	else if (!processed[i]) {

        		vertices[i] = new WeightedVertex(i,graph->computeMinFillRemovalCost(i));
        		verticesHeap.push_back(vertices[i]);
        		numInOrdering++;
        	}
        }

        make_heap(verticesHeap.begin(), verticesHeap.end(), CompareWeightedVertices(), UpdateWeightedVertexPos());

        // All vertices in the graph must be included in the ordering
        for (int i = 0; i < numInOrdering; i++) {

        	// ID of vertex to be deleted
            int vertexToBeDeleted = -1;
			bool isSimplical = false;

			// If the vertex at the top of the heap is simplical, remove it immediately
			if (verticesHeap.front()->getWeight() == 0) {
				isSimplical = true;
				vertexToBeDeleted = verticesHeap.front()->getVertex();

				pop_heap(verticesHeap.begin(), verticesHeap.end(), CompareWeightedVertices(), UpdateWeightedVertexPos());
				verticesHeap.pop_back();
			}
			// Otherwise, get the set of min cost vertices and select a vertex randomly from this set
			else {
				int minCost = verticesHeap.front()->getWeight();

				// Get the set of min cost vertices
				vector<WeightedVertex*> minCostVertices;
				while (!verticesHeap.empty() && (verticesHeap.front()->getWeight() == minCost)) {
					minCostVertices.push_back(verticesHeap.front());
					pop_heap (verticesHeap.begin(), verticesHeap.end(), CompareWeightedVertices(), UpdateWeightedVertexPos());
					verticesHeap.pop_back();
				}

				int cnt = 0;
				while (!verticesHeap.empty() && (cnt < pool)) {
					minCostVertices.push_back(verticesHeap.front());
					pop_heap (verticesHeap.begin(), verticesHeap.end(), CompareWeightedVertices(), UpdateWeightedVertexPos());
					verticesHeap.pop_back();
					cnt++;
				}

				int minIdx = 0;
				if (random != NULL) {
					minIdx = random->getInt(minCostVertices.size());
				}

				vertexToBeDeleted = minCostVertices[minIdx]->getVertex();

				for (int i = 0; i < (int)minCostVertices.size(); i++) {
					if (i != minIdx) {
						verticesHeap.push_back(minCostVertices[i]);
						push_heap(verticesHeap.begin(), verticesHeap.end(), CompareWeightedVertices(), UpdateWeightedVertexPos());
					}
				}
			}

            // Clique of the vertex being deleted
			vector<int> cluster;
			set<int> updateSet;

            // If this wasn't a simplical vertex, must add appropriate fill-in edges
            // and update removal costs
            if (!isSimplical) {

				// Prior to its removal, add an edge between non-adjacent neighbors
				// of the vertex to be deleted; record the vertices in the cluster
				// while doing so and also keep track of neighbors and neighbors of
            	// each neighbor as must update removal costs of these nodes
				vector<int>::iterator current, last, inner;
				current = graph->m_adjList[vertexToBeDeleted].begin();
				last = graph->m_adjList[vertexToBeDeleted].end();
				while (current != last) {
					// Add neighbor to current cluster
					cluster.push_back(*current);

					// Add neighbor to update list
					updateSet.insert(*current);

					// Add necessary edges
					inner = current;
					inner++;
					while (inner != last) {
						graph->addEdge(*current, *inner);
						inner++;
					}

					// Add neighbors of each neighbor to update list
					vector<int>::iterator it1, it2;
					it1 = graph->m_adjList[*current].begin();
					it2 = graph->m_adjList[*current].end();
					while (it1 != it2) {
						if (*it1 != vertexToBeDeleted) {
							updateSet.insert(*it1);
						}
						it1++;
					}

					current++;
				}
            }
            else {
            	// Update the removal costs of the neighbors only
            	vector<int>::iterator current, last;
				current = graph->m_adjList[vertexToBeDeleted].begin();
				last = graph->m_adjList[vertexToBeDeleted].end();
				while (current != last) {
					cluster.push_back(*current);
					updateSet.insert(*current);
					current++;
				}
            }

            // Add the vertex to be deleted to the current cluster
			cluster.push_back(vertexToBeDeleted);
			clusters.push_back(cluster);

            // Remove the vertex to be deleted from the graph, mark it as
			// processed and place it in the ordering
			processed[vertexToBeDeleted] = true;
			order.push_back(vertexToBeDeleted);
			graph->removeVertex(vertexToBeDeleted);

			// Update removal costs
			set<int>::iterator cur;
			cur = updateSet.begin();
			while (cur != updateSet.end()) {
				if (!processed[*cur]) {

					int newCost = graph->computeMinFillRemovalCost(*cur);
					if (newCost != vertices[*cur]->getWeight()) {
						vertices[*cur]->setWeight(newCost);
						update_heap_pos (verticesHeap.begin(), verticesHeap.end(), verticesHeap.begin() + vertices[*cur]->heap_pos(), CompareWeightedVertices(), UpdateWeightedVertexPos());
					}
				}
				cur++;
			}
        }

        delete graph;

        for (int i = 0; i < (int)evidence.size(); i++) {
        	order.push_back(evidence[i]);
        }
        return order;
    }

    /**
     * Gets a min-fill ordering of the vertices in this graph with degree > 0.  Assumes that the
     * graph of such vertices is connected. Populates <tt>order</tt> with the discovered elimination order.
     */
    static vector<int> getMinFillOrdering(Graph* inGraph, myRandom* random, vector<vector<int> >& clusters, int pool) {
    	vector<int> evidence;
    	return getMinFillOrdering(inGraph, evidence, random, clusters, pool);
    }

    /**
     * Gets the set of maximal cliques in this graph using the Bron-Kerbosch algorithm.
     */
    vector<vector<int> > getMaximalCliques() {
    	m_maximalCliques.clear();
    	vector<int> R;
    	vector<int> P = vector<int>(m_numVertices);
    	for (int i = 0; i < m_numVertices; i++) {
    		P[i] = i;
    	}
    	vector<int> X;
    	BronKerbosch(R,P,X);

    	return m_maximalCliques;
    }

    /**
     * Gets an elimination order using the min-fill heuristic.  The elimination order is returned in <tt>order</tt>
     * and the set of clusters created along this ordering is returned in <tt>clusters</tt>.
     * @param random - Random number generator
     * @param pool - The number of vertices added to the pool of min-fill vertices the next vertex to be removed is
     *               drawn from.  So, if the min-fill vertex set includes vertices {0,3,9} all with a cost of 1, and
     *               pool is 2, then the next two lowest cost vertices will be added to the pool, increasing the pool
     *               of vertices the vertex to be removed is drawn from 3 to 5.
     */
    int getMinFillOrdering(vector<int>& order, vector<vector<int> >& clusters, myRandom* random, int pool) {
		int width = 0;
		assert(pool >= 0);
		order.clear();
		clusters.clear();

		// Keep track of the variables that have been processed/added to the ordering
		vector<bool> processed = vector<bool>(m_numVertices);
		std::fill(processed.begin(), processed.end(), false);

		// Pre-compute the removal costs
		vector<WeightedVertex*> vertices;
		vector<WeightedVertex*> verticesHeap;
		for (int i = 0; i < m_numVertices; i++) {
			vertices.push_back(new WeightedVertex(i, computeMinFillRemovalCost(i)));
			verticesHeap.push_back(vertices[i]);
		}

		make_heap(verticesHeap.begin(), verticesHeap.end(), CompareWeightedVertices(), UpdateWeightedVertexPos());

		// All vertices in the graph must be included in the ordering
		for (int i = 0; i < m_numVertices; i++) {

			// ID of vertex to be deleted
			int vertexToBeDeleted = -1;
			bool isSimplical = false;

			// Get the vertex at the top of the heap
			// If the vertex is simplical, remove it immediately
			if (verticesHeap.front()->getWeight() == 0) {
				isSimplical = true;
				vertexToBeDeleted = verticesHeap.front()->getVertex();

				pop_heap(verticesHeap.begin(), verticesHeap.end(), CompareWeightedVertices(), UpdateWeightedVertexPos());
				verticesHeap.pop_back();
			}
			// Otherwise, get the set of min cost vertices and select a vertex randomly from this set
			else {
				int minCost = verticesHeap.front()->getWeight();

				// Get the set of min cost vertices
				vector<WeightedVertex*> minCostVertices;
				while (!verticesHeap.empty() && (verticesHeap.front()->getWeight() == minCost)) {
					minCostVertices.push_back(verticesHeap.front());
					pop_heap (verticesHeap.begin(), verticesHeap.end(), CompareWeightedVertices(), UpdateWeightedVertexPos());
					verticesHeap.pop_back();
				}

				int cnt = 0;
				while (!verticesHeap.empty() && (cnt < pool)) {
					minCostVertices.push_back(verticesHeap.front());
					pop_heap (verticesHeap.begin(), verticesHeap.end(), CompareWeightedVertices(), UpdateWeightedVertexPos());
					verticesHeap.pop_back();
					cnt++;
				}

				int minIdx = 0;
				if (random != NULL) {
					minIdx = random->getInt(minCostVertices.size());
				}

				vertexToBeDeleted = minCostVertices[minIdx]->getVertex();

				for (int i = 0; i < (int)minCostVertices.size(); i++) {
					if (i != minIdx) {
						verticesHeap.push_back(minCostVertices[i]);
						push_heap(verticesHeap.begin(), verticesHeap.end(), CompareWeightedVertices(), UpdateWeightedVertexPos());
					}
				}
			}

			// Update the width of the current ordering given deletion of the
			// vertexToBeDeleted. NOTE that adjacency lists do not include the
			// vertex being deleted
			if (int(m_adjList[vertexToBeDeleted].size()) > width) {
				width = int(m_adjList[vertexToBeDeleted].size());
			}

			// Clique of the vertex being deleted
			vector<int> cluster;
			set<int> updateSet;

			// If this wasn't a simplical vertex, must add appropriate fill-in edges
			// and update removal costs
			if (!isSimplical) {

				// Prior to its removal, add an edge between non-adjacent neighbors
				// of the vertex to be deleted; record the vertices in the cluster
				// while doing so and also keep track of neighbors and neighbors of
				// each neighbor as must update removal costs of these nodes
				vector<int>::iterator current, last, inner;
				current = m_adjList[vertexToBeDeleted].begin();
				last = m_adjList[vertexToBeDeleted].end();
				while (current != last) {
					// Add neighbor to current cluster
					cluster.push_back(*current);

					// Add neighbor to update list
					updateSet.insert(*current);

					// Add necessary edges
					inner = current;
					inner++;
					while (inner != last) {
						addEdge(*current, *inner);
						inner++;
					}

					// Add neighbors of each neighbor to update list
					// too
					vector<int>::iterator it1, it2;
					it1 = m_adjList[*current].begin();
					it2 = m_adjList[*current].end();
					while (it1 != it2) {
						if (*it1 != vertexToBeDeleted) {
							updateSet.insert(*it1);
						}
						it1++;
					}

					current++;
				}
			}
			else {
				// Update the removal costs of the neighbors only
				vector<int>::iterator current, last;
				current = m_adjList[vertexToBeDeleted].begin();
				last = m_adjList[vertexToBeDeleted].end();
				while (current != last) {
					cluster.push_back(*current);
					updateSet.insert(*current);
					current++;
				}
			}

			// Add the vertex to be deleted to the current cluster
			cluster.push_back(vertexToBeDeleted);
			clusters.push_back(cluster);

			// Remove the vertex to be deleted from the graph, mark it as
			// processed and place it in the ordering
			processed[vertexToBeDeleted] = true;
			order.push_back(vertexToBeDeleted);
			removeVertex(vertexToBeDeleted);

			// Update removal costs
			set<int>::iterator cur;
			cur = updateSet.begin();
			while (cur != updateSet.end()) {
				if (!processed[*cur]) {

					// Compute the new cost and update if needed
					int newCost = computeMinFillRemovalCost(*cur);
					if (vertices[*cur]->getWeight() != newCost) {
						vertices[*cur]->setWeight(newCost);
						update_heap_pos (verticesHeap.begin(), verticesHeap.end(), verticesHeap.begin() + vertices[*cur]->heap_pos(), CompareWeightedVertices(), UpdateWeightedVertexPos());
					}
				}
				cur++;
			}
		}

		return width;
    }

    /**
     * Computes the width of the specified elimination ordering and also determines
     * the sequence of clusters produced by the specified ordering.
     * @param order vector<int> --
     * @param clusters vector<vector<int> > -
     * @return width -
     */
    int computeWidth(vector<int>& order, vector<vector<int> >& clusters) {
    	int width = 0;
    	clusters.clear();

    	for (int i = 0; i < (int)order.size(); i++) {
    		int vertexToBeDeleted = order[i];

			// Update the width of the current ordering given deletion of the
			// vertexToBeDeleted. NOTE that adjacency lists do not include the
			// vertex being deleted
			if (int(m_adjList[vertexToBeDeleted].size()) > width) {
				width = int(m_adjList[vertexToBeDeleted].size());
			}

			vector<int> cluster;

			// Prior to its removal, add an edge between non-adjacent neighbors
			// of the vertex to be deleted; record the vertices in the cluster
			// while doing so
			vector<int>::iterator current, last, inner;
			current = m_adjList[vertexToBeDeleted].begin();
			last = m_adjList[vertexToBeDeleted].end();
			while (current != last) {
				// Add neighbor to current cluster
				cluster.push_back(*current);

				// Add necessary edges
				inner = current;
				inner++;
				while (inner != last) {
					addEdge(*current, *inner);
					inner++;
				}

				current++;
			}

            // Add the vertex to be deleted to the current cluster
            cluster.push_back(vertexToBeDeleted);
            clusters.push_back(cluster);

            // Remove the vertex to be deleted from the graph
			removeVertex(vertexToBeDeleted);
    	}

    	return width;
    }

    /**
     * Creates a sorted edge set for the specified cycle
     */
    vector<int> getEdgeSet(vector<int> cycle) {

    	vector<int> edgeSet;

    	int u = cycle[0];
		int v = cycle[((int)cycle.size())-1];

		if (u < v) {
			edgeSet.push_back(u*m_numVertices + v);
		}
		else {
			edgeSet.push_back(v*m_numVertices + u);
		}

    	for (int i = 1; i < (int)cycle.size(); i++) {
    		u = cycle[i-1];
    		v = cycle[i];

    		if (u < v) {
				edgeSet.push_back(u*m_numVertices + v);
			}
			else {
				edgeSet.push_back(v*m_numVertices + u);
			}
    	}

    	sort(edgeSet.begin(), edgeSet.end());

    	return edgeSet;
    }

    /**
     * Gets a maximum weight spanning tree on the specified graph using Prim's algorithm.
     * Assumes that the graph is connected.
     * @param root - vertex to grow spanning tree from
     * @param componentSize - number of vertices in the connected component containing root
     * @return vector<pair<int,int> > - A set of edges comprising the max weight spanning tree
     */
    vector<pair<int,int> > getMaxWeightSpanningTree(int root, int componentSize) {
    	assert((root >= 0) && (root < m_numVertices));
    	assert(m_isWeighted);

    	// Get the maximum weight edge in this graph
    	map<int,int>::iterator it;
    	it = m_edgeIDToWeightMap.begin();
    	int maxWeight = it->second;
    	it++;
    	while (it != m_edgeIDToWeightMap.end()) {
    		if (it->second > maxWeight) {
    			maxWeight = it->second;
    		}
    		it++;
    	}

    	// Create the spanning tree object
    	vector<pair<int,int> > spanningTree;

    	// Create the priority queue data structure
    	priority_queue<WeightedEdge, vector<WeightedEdge >, CompareWeightedEdges > priorityQueue;

    	// Mark all vertices as unvisited
    	vector<bool> visited = vector<bool>(m_numVertices);
    	fill(visited.begin(),visited.end(),false);
    	int numAdded = 0;

    	// Build spanning tree starting from the root
    	visited[root] = true;
    	numAdded++;
    	// Get each edge from 0
    	for (int i = 0; i < (int)m_adjList[root].size(); i++) {
    		// Get the weight of edge root-v
    		int v = m_adjList[root][i];

    		int key = m_numVertices*root + v;

    		if (v < root) {
    			key = m_numVertices*v + root;
    		}

    		int weight = maxWeight - m_edgeIDToWeightMap[key];

    		// Add this edge to the priority queue
    		WeightedEdge edge(root,v,weight);

//cout << "Pushing edge " << root << "-" << v << ":" << weight << " to queue" << endl;
    		priorityQueue.push(edge);
    	}

    	// While all vertices have not been visited
    	while (numAdded < componentSize) {
    		// Ensure the graph is connected
    		assert(!priorityQueue.empty());

    		// Get the next min weight edge from some visited node, u, to some unvisited node v
    		while(visited[priorityQueue.top().m_v]) {
//cout << "Ignoring edge " << priorityQueue.top().m_u << "-" << priorityQueue.top().m_v << endl;
    			priorityQueue.pop();

    			// If the queue is empty...we are done
    			assert(!priorityQueue.empty());
    		}

    		// Get the vertices of this edge and remove the edge from the queue
    		int u = priorityQueue.top().m_u;
    		int v = priorityQueue.top().m_v;
//cout << "Adding edge " << u << "-" << v << ":" << priorityQueue.top().m_weight << " from queue" << endl;
    		priorityQueue.pop();

    		// Mark node v as visited
    		visited[v] = true;
    		numAdded++;

    		// Add the edge to the return object
    		pair<int,int> stedge(u,v);
    		spanningTree.push_back(stedge);

    		// Add each edge v-s to the priority queue so long as s has not been visited
    		for (int i = 0; i < (int)m_adjList[v].size(); i++) {
    			int s = m_adjList[v][i];
    			if (!visited[s]) {
    				int key = v*m_numVertices + s;
    				if (s < v) {
    					key = s*m_numVertices + v;
    				}

    	    		int weight = maxWeight - m_edgeIDToWeightMap[key];

//cout << "Pushing edge " << v << "-" << s << ":" << weight << " to queue" << endl;
    	    		// Add this edge to the priority queue
    	    		WeightedEdge edge(v,s,weight);
    	    		priorityQueue.push(edge);
    			}
    		}
    	}

    	return spanningTree;
    }
};

#endif /* GRAPH_H_ */
