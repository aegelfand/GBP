/*
 * JoinGraph.h
 *
 *  Created on: Sep 26, 2011
 *      Author: agelfand
 */

#ifndef JOINGRAPH_H_
#define JOINGRAPH_H_

#include <algorithm>
#include <sstream>
#include <utility>
#include <vector>

#include "LogFunction.h"
#include "Graph.h"
#include "myRandom.h"
#include "Util.h"
#include "Variable.h"

using namespace std;

struct JGEdge;

struct JGNode {
protected:
	/**
	 * The set of variables in this JG node
	 */
	vector<Variable*> m_variables;

	/**
	 * The set of original functions assigned to this JG node
	 */
	vector<LogFunction*> m_originalFunctions;

	/**
	 * The set of edges to/from this JG node
	 */
	vector<JGEdge*> m_edges;

    /**
     * ID of this cluster
     */
    int m_ID;

    /**
     * Count of clusters
     */
    static int s_count;

public:

	/**
	 * Default constructor
	 * NOTE:  Assumes that the variables are sorted in ascending order by ID.
	 */
	JGNode(vector<Variable*> variables) {
		m_variables = variables;
		m_ID = s_count++;
	}

	JGNode(Variable* var) {
		m_variables.push_back(var);
		m_ID = s_count++;
	}

	/**
	 * Adds the specified set of variables to this JG node
	 */
	void addVariables(vector<Variable*> variables) {
		getUnionVar(m_variables,variables,m_variables);
	}

    /**
     * Gets the variables in this JG node
     */
    const vector<Variable*>& getVariables() {
        return m_variables;
    }

    /**
     * Gets the set of edges of this JG node
     */
    const vector<JGEdge*>& getEdges() {
    	return m_edges;
    }
	/**
	 * Override of the equals operator. Comparison is based on JG node ID.
	 */
	bool operator==(const JGNode& other) { return (m_ID == other.m_ID); }

    /**
     * Assigns a function to this region.
     */
    void addFunction(LogFunction* function) {
    	// Ensure function scope is a subseteq m_variables
    	vector<Variable*> isect;
    	getIntersectionVar(m_variables, function->getScope(), isect);
    	assert(isect.size() == function->getScope().size());

    	m_originalFunctions.push_back(function);
    }

    /**
     * Gets the set of functions assigned to this JG node
     */
    const vector<LogFunction*>& getOriginalFunctions() {
    	return m_originalFunctions;
    }

    int getID() {
    	return m_ID;
    }

    /**
     * Adds an edge to this JG node
     */
    void addEdge(JGEdge* edge);

    /**
     * Removes the JG edge from this JG node
     */
    void removeEdge(JGEdge* edge);

};

struct JGEdge {
protected:
    /**
     * JG node i
     */
    JGNode* m_node_i;

    /**
     * JG node j
     */
    JGNode* m_node_j;

	/**
	 * The set of variables this edge is defined over (separator).
	 */
	vector<Variable*> m_separator;

    /**
     * ID of this edge
     */
    int m_ID;

    /**
     * Count of edges
     */
    static int s_count;

public:

    /**
	 * Default constructor
	 */
	JGEdge(JGNode* node_i, JGNode* node_j, vector<Variable*> separator) {
		if (node_i->getVariables().size() >= node_j->getVariables().size()) {
			m_node_i = node_i;
			m_node_j = node_j;
		}
		else {
			m_node_j = node_i;
			m_node_i = node_j;
		}

		m_node_i->addEdge(this);
		m_node_j->addEdge(this);
		m_ID = s_count++;

		// Validate the separator
		m_separator = separator;
		assert(!m_separator.empty());
		vector<Variable*> isect;
		getIntersectionVar(m_node_i->getVariables(), m_separator, isect);
		assert(isect.size() == m_separator.size());
		getIntersectionVar(m_node_j->getVariables(), m_separator, isect);
		assert(isect.size() == m_separator.size());
	}

    /**
     * Switches nodes i and j.
     */
    void swap() {
    	JGNode* temp = m_node_i;
    	m_node_i = m_node_j;
    	m_node_j = temp;
    }

	/**
	 * Override of the equals operator. Comparison is based on edge ID.
	 */
	bool operator==(const JGEdge& other) { return (m_ID == other.m_ID); }
	bool operator !=(const JGEdge& other) { return (m_ID != other.m_ID); }

	/**
	 * Gets JG node i
	 */
	JGNode* getNode_i() {
		return m_node_i;
	}

	/**
	 * Gets JG node j
	 */
	JGNode* getNode_j() {
		return m_node_j;
	}

	/**
	 * Gets the variables on this edge (i.e. the separator)
	 */
	const vector<Variable*>& getVariables() {
		return m_separator;
	}

	/**
	 * Replaces node_i with the specified node.
	 */
	void replaceNode_i(JGNode* node_i) {
		m_node_i->removeEdge(this);
		m_node_i = node_i;
		m_node_i->addEdge(this);
	}

	/**
	 * Replaces node_j with the specified node.
	 */
	void replaceNode_j(JGNode* node_j) {
		m_node_j->removeEdge(this);
		m_node_j = node_j;
		m_node_j->addEdge(this);
	}
};

/**
 * Implementation of a message passed along the Join Graph.
 * There are two types of JG messages.  The first is a just a
 * function assigned to a bucket; the other is an actual message
 * sent from a previous bucket
 */
struct JGMessage {
protected:
	/**
	 * The node sending the message
	 */
	JGNode* m_node;

	/**
	 * The index of the bucket from which the message was sent
	 */
	int m_sendingBucketID;

	const static int ASSIGNED_FUNCTION_ID = -1;

	/**
	 * The scope of this message
	 */
	vector<Variable*> m_scope;

	/**
	 * The function associated with this message
	 */
	LogFunction* m_function;

public:

	/**
	 * Creates a JGMessage object for an assigned function
	 */
	JGMessage(LogFunction* function) {
		m_sendingBucketID = ASSIGNED_FUNCTION_ID;
		m_function = function;
		m_scope = m_function->getScope();
		m_node = NULL;
	}

	/**
	 * Creates a JGMessage object
	 */
	JGMessage(JGNode* sendingNode, int indexOfSendingBucket, vector<Variable*> msgScope) {
		m_function = NULL;
		m_node = sendingNode;
		m_sendingBucketID = indexOfSendingBucket;
		m_scope = msgScope;

	}

	bool isAssignedFunction() {
		return (m_sendingBucketID == ASSIGNED_FUNCTION_ID);
	}

	const vector<Variable*>& getScope() {
		return m_scope;
	}

	int scopeSize() {
		return getScope().size();
	}

	JGNode* getSendingNode() {
		return m_node;
	}

	LogFunction* getFunction() {
		return m_function;
	}

	string toString() {
		std::stringstream ss;
		if (isAssignedFunction()) {
			ss << "Function " << varsToString(m_scope);
		}
		else {
			ss << "Msg over " << varsToString(m_scope) << " from cluster " << varsToString(m_node->getVariables()) << " in bucket[" << m_sendingBucketID << "]";
		}
		return ss.str();
	}

	bool operator==(const JGMessage& other) {
		return (m_scope.size() == other.m_scope.size());
	}
	bool operator<(const JGMessage& other) {
		return (m_scope.size() > other.m_scope.size());
	}
};

static bool compareByScopeSize(JGMessage* A, JGMessage* B) {
	return (A->getScope().size() > B->getScope().size());
}

struct JoinGraph {
protected:
	/**
	 * The set of edges in this join graph
	 */
	vector<JGEdge*> m_edges;

	/**
	 * The set of nodes in this join graph
	 */
	vector<JGNode*> m_nodes;

public:

	/**
	 * Default join graph constructor;
	 */
	JoinGraph(vector<JGNode*> nodes, vector<JGEdge*> edges) {
		m_nodes = nodes;
		m_edges = edges;
	}

	const vector<JGNode*>& getNodes() {
		return m_nodes;
	}

	const vector<JGEdge*>& getEdges() {
		return m_edges;
	}

	string toString() {
		std::stringstream ss;
		for (int i = 0; i < (int)m_nodes.size(); i++) {
			ss << "Cluster[" << i << "]: " << varsToString(m_nodes[i]->getVariables()) << endl;
		}
		for (int i = 0; i < (int)m_edges.size(); i++) {
			ss << "Edge[" << i << "]: " << varsToString(m_edges[i]->getNode_i()->getVariables());
			ss << "-" << varsToString(m_edges[i]->getNode_j()->getVariables()) << " ";
			ss << varsToString(m_edges[i]->getVariables()) << endl;
		}
		return ss.str();
	}

	/**
	 * Simplifies the join graph by removing redundant clusters.
	 */
	void simplify() {

		// In reverse order along the bucket tree, look for separators and
		// clusters defined on the same set of variables. If a cluster and
		// separator on the same set of variables is found, the cluster is
		// redundant and can be removed.
		for (int i = (int)m_nodes.size() - 1; i >= 0; i--) {
			for (int j = (int)m_edges.size() - 1; j >= 0; j--) {

				// If cluster i is redundant, merge it w/ cluster at other end of edge j
				if (((m_edges[j]->getNode_i() == m_nodes[i]) ||
					 (m_edges[j]->getNode_j() == m_nodes[i])) &&
					(isEqualVar(m_edges[j]->getVariables(),m_nodes[i]->getVariables()))) {

					JGNode* otherCluster;
					if (m_edges[j]->getNode_i() == m_nodes[i]) {
						otherCluster = m_edges[j]->getNode_j();
					}
					else {
						assert(m_edges[j]->getNode_j() == m_nodes[i]);
						otherCluster = m_edges[j]->getNode_i();
					}

					for (int m = (int)m_nodes[i]->getEdges().size() - 1; m >= 0; m--) {

						JGEdge* edge = m_nodes[i]->getEdges()[m];

						if (edge != m_edges[j]) {

							if (edge->getNode_i() == m_nodes[i]) {
								edge->replaceNode_i(otherCluster);
							}
							else {
								assert(edge->getNode_j() == m_nodes[i]);
								edge->replaceNode_j(otherCluster);
							}

						}
						else {
							otherCluster->removeEdge(edge);
							m_nodes[i]->removeEdge(edge);
						}
					}

					for (int m = 0; m < (int)m_nodes[i]->getOriginalFunctions().size(); m++) {
						otherCluster->addFunction(m_nodes[i]->getOriginalFunctions()[m]);
					}

					// remove edge j
					m_edges.erase(m_edges.begin() + j);

					// remove cluster i
					m_nodes.erase(m_nodes.begin() + i);

					break;

				} // if cluster i is redundant
			} // for each edge
		} // for each cluster
	}

	static JoinGraph* runMiniBucketSchematic(vector<int> order, vector<Variable*>& variables, vector<LogFunction*>& functions, int iBound) {
		vector<JGEdge*> edges;
		vector<JGNode*> nodes;

		int numVariables = variables.size();

		//
		// Assign functions to buckets.
		//
		// Each bucket is associated with a variable and buckets are ordered
		// according to the specified elimination order. Place each function i
		// in the first bucket j where variable j occurs in function i's scope.
		vector<vector<JGMessage*> > messagesInBucket = vector<vector<JGMessage*> >(numVariables);
		for (int i = 0; i < (int)functions.size(); i++) {
			bool assigned = false;
			for (int j = 0; j < numVariables; j++) {
				if (containsVar(functions[i]->getScope(), variables[order[j]])) {
					JGMessage* newJGMsg = new JGMessage(functions[i]);

					messagesInBucket[j].push_back(newJGMsg);
					assigned = true;
					break;
				}
			}
			assert(assigned);
		}

		//
		// The mini-bucket tree structure
		//
		// The bucket tree has one slot for each variable in the ordering. Each bucket can contain
		// a set of clusters, whrere each cluster is a set of variables.
		vector<vector<JGNode*> > clustersInBucket = vector<vector<JGNode*> >(numVariables);

		// At each bucket i, partition the functions assigned to that bucket
		// and any messages from previous buckets into mini-buckets (clusters)
		// of at most iBound variables.  Perform this partitioning in a greedy
		// fashion based upon message size. Then simulate passing of messages
		// from each mini-bucket to some future bucket
		for (int i = 0; i < numVariables; i++) {

			int numMsgsInBucket_i = messagesInBucket[i].size();

			// Sort the messages in bucket i largest to smallest by scope size
			sort(messagesInBucket[i].begin(),messagesInBucket[i].end(),compareByScopeSize);

			//
			// Begin the partitioning process. Greedily add functions to the current mini-bucket
			// until we surpass the specified iBound
			//
			int curCluster = -1;
			vector<bool> assignedToCluster = vector<bool>(numMsgsInBucket_i);
			fill(assignedToCluster.begin(),assignedToCluster.end(),false);
			int numMsgsProcessed = 0;
			while (numMsgsProcessed < numMsgsInBucket_i) {

				// Determine the index of the first message in bucket i that has not been processed
				int idx = 0;
				for (idx = 0; idx < numMsgsInBucket_i; idx++) {
					if (!assignedToCluster[idx]) {
						break;
					}
				}

				JGNode* newNode = new JGNode(messagesInBucket[i][idx]->getScope());
				if (messagesInBucket[i][idx]->isAssignedFunction()) {
					newNode->addFunction(messagesInBucket[i][idx]->getFunction());
				}
				else {
					JGEdge* newEdge = new JGEdge(messagesInBucket[i][idx]->getSendingNode(), newNode, messagesInBucket[i][idx]->getScope());
					edges.push_back(newEdge);
				}

				// Add the new cluster to bucket i
				clustersInBucket[i].push_back(newNode);
				nodes.push_back(newNode);
				curCluster++;
				assignedToCluster[idx] = true;
				numMsgsProcessed++;

				//
				// Greedily add the remaining, available function to the current cluster
				//
				for (int j = 1; j < numMsgsInBucket_i; j++) {

					if (assignedToCluster[j]) {
						continue;
					}

					vector<Variable*> newScope;
					getUnionVar(clustersInBucket[i][curCluster]->getVariables(),messagesInBucket[i][j]->getScope(),newScope);

					if (iBound >= ((int)newScope.size())) {

						clustersInBucket[i][curCluster]->addVariables(messagesInBucket[i][j]->getScope());

						if (messagesInBucket[i][j]->isAssignedFunction()) {
							clustersInBucket[i][curCluster]->addFunction(messagesInBucket[i][j]->getFunction());
						}
						else {
							JGEdge* newEdge = new JGEdge(messagesInBucket[i][j]->getSendingNode(), clustersInBucket[i][curCluster], messagesInBucket[i][j]->getScope());
							edges.push_back(newEdge);
						}

						numMsgsProcessed++;
						assignedToCluster[j] = true;
					}
				}
			} // while messages need to be processed

			// Now perform the forwarding operation. A message from each cluster is generated
			// by removing the variable in bucket i.  This message is forwarded to the first
			// bucket j that appears in this message scope
			for (int k = 0; k < (int)clustersInBucket[i].size(); k++) {
				for (int j = i + 1; j < numVariables; j++) {
					if (containsVar(clustersInBucket[i][k]->getVariables(), variables[order[j]])) {

						// Create the message
						vector<Variable*> msgScope;
						removeVar(clustersInBucket[i][k]->getVariables(), variables[order[i]], msgScope);

						// Don't propagate empty messages
						if (!msgScope.empty()) {
							JGMessage* newJGMsg = new JGMessage(clustersInBucket[i][k], i, msgScope);

							// And add to bucket j
							messagesInBucket[j].push_back(newJGMsg);
						}

						break;
					}
				}
			}

		} // for each bucket (i = 0; i < numVariables)

		// Last connect the clusters in each bucket
		for (int i = 0; i < numVariables; i++) {
			for (int j = 1; j < (int)clustersInBucket[i].size(); j++) {
				vector<Variable*> sep;
				sep.push_back(variables[order[i]]);
				JGEdge* newEdge = new JGEdge(clustersInBucket[i][j-1],clustersInBucket[i][j],sep);
				edges.push_back(newEdge);
			}
		}

		return new JoinGraph(nodes, edges);
	}

};


#endif /* JOINGRAPH_H_ */
