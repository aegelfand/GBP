/*
 * Partition.h
 *
 *  Created on: Feb 6, 2013
 *      Author: root
 */

#ifndef PARTITION_H_
#define PARTITION_H_

#include <algorithm>
#include <sstream>
#include <string>
#include <vector>

#include "LogFunction.h"
#include "Util.h"
#include "Variable.h"

struct ClusterEdge;

struct Cluster {
protected:
	/**
	 * The set of variables in this region (in ascending order by variable ID)
	 */
	vector<Variable*> m_variables;

	/**
	 * The set of functions assigned to this region.
	 */
	vector<LogFunction*> m_functions;

	/**
	 * The set of edges connecting this cluster to other clusters
	 */
	vector<ClusterEdge*> m_edges;

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
	 * Compares clusters a and b based upon their IDs.
	 */
	static inline bool compareClusterByID(Cluster* a, Cluster* b) {
		return (a->m_ID < b->m_ID);
	}

	/**
	 * Gets the intersection of sets A and B and returns the result in C - i.e. C = A cap B
	 * Assumes that both <tt>A</tt> and <tt>B</tt> are sorted.
	 * @return bool -- Returns <tt>true</tt> if the intersection is not empty; <tt>false</tt> otherwise;
	 */
	static inline bool getIntersectionCluster(const vector<Cluster*>& A, const vector<Cluster*>& B, vector<Cluster*>& C) {
		vector<Cluster*> temp;
		temp.resize(A.size());
		vector<Cluster*>::iterator curr,last;
		last = temp.end();
		curr = set_intersection(A.begin(), A.end(), B.begin(), B.end(), temp.begin(), compareClusterByID);
		temp.erase(curr,last);
		C = temp;
		return !C.empty();
	}

	/**
	 * Gets the union of sets A and B and returns the result in C - i.e. C = A cup B
	 * Assumes that both <tt>A</tt> and <tt>B</tt> are sorted.
	 */
	static inline void getUnionCluster(const vector<Cluster*>& A, const vector<Cluster*>& B, vector<Cluster*>& C) {
		vector<Cluster*> temp;
		temp.resize(A.size() + B.size());
		vector<Cluster*>::iterator curr,last;
		last = temp.end();
		curr = set_union(A.begin(), A.end(), B.begin(), B.end(), temp.begin(), compareClusterByID);
		temp.erase(curr,last);
		C = temp;
	}

	/**
	 * Adds the specified cluster B to the set A
	 */
	static inline void addCluster(vector<Cluster*>& A, Cluster* B) {
		vector<Cluster*> temp;
		temp.push_back(B);
		getUnionCluster(A,temp,A);
	}

	/**
	 * Creates a new cluster
	 */
    Cluster(LogFunction* function) {
    	m_variables = vector<Variable*>(function->getScope());
		m_functions.push_back(function);
		m_ID = s_count++;
    }

    /**
     * Creates an empty cluster over the specified set of variables
     */
    Cluster(vector<Variable*> variables) {
    	m_variables = vector<Variable*>(variables);
		m_ID = s_count++;
    }

    /**
     * Creates a new cluster by merging the two specified clusters
     */
    Cluster(Cluster& a, Cluster& b) {
    	getUnionVar(a.m_variables, b.m_variables, m_variables);
    	for (int i = 0; i < (int)a.m_functions.size(); i++) {
    		m_functions.push_back(a.m_functions[i]);
    	}
    	for (int i = 0; i < (int)b.m_functions.size(); i++) {
			m_functions.push_back(b.m_functions[i]);
		}
    	m_ID = s_count++;
    }

    /**
     * Adds a function to this cluster
     */
    void addFunction(LogFunction* func) {
    	assert(isSubsetVar(func->getScope(), m_variables));
    	m_functions.push_back(func);
    }

    /**
     * Adds a set of variables to this cluster
     */
    void addVariables(vector<Variable*> variables) {
    	getUnionVar(m_variables, variables, m_variables);
    }


    /**
     * Gets the set of variables defining this cluster
     */
    const vector<Variable*>& getVariables() {
    	return m_variables;
    }

    /**
     * Gets the set of functions assigned to this cluster
     */
    const vector<LogFunction*>& getFunctions() {
    	return m_functions;
    }

    /**
     * Adds an edge to this cluster
     */
    void addEdge(ClusterEdge* edge);

    /**
	 * Removes a neighbor to this cluster
	 */
	void removeEdge(ClusterEdge* edge);

    /**
     * Gets the neighbors of this cluster
     */
    const vector<ClusterEdge*>& getEdges() {
    	return m_edges;
    }

    /**
     * Gets the neighbors of this cluster whose cluster
     * size is greater than the specified iBound
     */
    vector<ClusterEdge*> getEdges(int iBound);

    /**
     * Gets the ID of this cluster
     */
    int getID() {
    	return m_ID;
    }

    /**
     * Gets a string representation of this cluster and its neighbors
     */
    string toString();

    /**
     * Returns true if cluster u is a neighbor of this cluster
     */
    bool isNeighbor(Cluster* u);

    /**
     * Merges the specified cluster into this cluster.
     * If isMergeEdges is <tt>true</tt> the edges of
     * this cluster and other cluster are merged too.
     */
    void merge(Cluster* other, bool isMergeEdges);

};

/**
 * Implementation of a cluster edge structure
 */
struct ClusterEdge {
protected:
	/**
	 * The current score of this edge
	 */
	double m_score;

	/**
	 * Cluster i on this edge
	 */
	Cluster* m_cluster_i;

	/**
	 * Cluster j on this edge
	 */
	Cluster* m_cluster_j;

	/**
	 * Flag indicating that this edge needs to be updated
	 */
	bool m_updateNeeded;

	/**
	 * ID of this cluster edge
	 */
	int m_ID;

	/**
	 * Count of clusters
	 */
	static int s_count;

public:

	/**
	 * Constructor
	 */
	ClusterEdge(Cluster* cluster_i, Cluster* cluster_j) {
		m_cluster_i = cluster_i;
		m_cluster_j = cluster_j;
		m_updateNeeded = true;
		m_score = -1;
		m_ID = s_count++;
		m_cluster_i->addEdge(this);
		m_cluster_j->addEdge(this);
	}

	/**
	 * Gets the ID of this cluster edge
	 */
	int getID() {
		return m_ID;
	}

	/**
	 * Gets the score on this edge
	 */
	double getScore() {
		return m_score;
	}

	/**
	 * Sets the score on this edge
	 */
	void setScore(double score) {
		m_score = score;
	}

	/**
	 * Returns <tt>true</tt> if this edge needs to be updated;
	 * <tt>false</tt> otherwise.
	 */
	bool isUpdateNeeded() {
		return m_updateNeeded;
	}

	/**
	 * Sets the update status flag
	 */
	void setUpdateStatus(bool updateNeeded) {
		m_updateNeeded = updateNeeded;
	}

	/**
	 * Gets cluster i in this edge
	 */
	Cluster* getCluster_i() {
		return m_cluster_i;
	}

	/**
	 * Gets cluster j in this edge
	 */
	Cluster* getCluster_j() {
		return m_cluster_j;
	}

	/**
	 * Gets the other cluster on this edge, i.e. given cluster
	 * u on edge (u,v) returns cluster v.
	 * Returns <tt>NULL</tt> if u is not a cluster in this edge.
	 */
	Cluster* getOtherCluster(Cluster* u) {
		if (m_cluster_i->getID() == u->getID()) {
			return m_cluster_j;
		}
		else if(m_cluster_j->getID() == u->getID()) {
			return m_cluster_i;
		}
		return NULL;
	}

	/**
	 * Swaps cluster u on this edge for cluster t, i.e. given edge
	 * (u,v) and cluster's u and t creates the edge (t,v).
	 * Returns <tt>true</tt> if edge contains cluster u; <tt>false</tt>
	 * otherwise.
	 */
	bool swap(Cluster* u, Cluster* t) {
		if (m_cluster_i->getID() == u->getID()) {
			m_cluster_i = t;
			u->removeEdge(this);
			t->addEdge(this);
			return true;
		}
		else if (m_cluster_j->getID() == u->getID()) {
			m_cluster_j = t;
			u->removeEdge(this);
			t->addEdge(this);
			return true;
		}
		return false;
	}

	/**
	 * Returns <tt>true</tt> if this edge contains cluster u;
	 * <tt>false</tt> otherwise.
	 */
	bool containsCluster(Cluster* u) {
		if (m_cluster_i->getID() == u->getID()) {
			return true;
		}
		else if (m_cluster_j->getID() == u->getID()) {
			return true;
		}
		return false;
	}

	/**
	 * Gets a string representation of this edge
	 */
	string toString() {
		stringstream ss;
		if (m_cluster_i->getID() < m_cluster_j->getID()) {
			ss << m_cluster_i->getID() << ":(" << varsToString(m_cluster_i->getVariables()) << ") - ";
			ss << m_cluster_j->getID() << ":(" << varsToString(m_cluster_j->getVariables()) << ")";
		}
		else {
			ss << m_cluster_j->getID() << ":(" << varsToString(m_cluster_j->getVariables()) << ") - ";
			ss << m_cluster_i->getID() << ":(" << varsToString(m_cluster_i->getVariables()) << ")";
		}
		return ss.str();
	}

};

static bool compareClustersBySize(Cluster* a, Cluster* b) {
	return (a->getVariables().size() > b->getVariables().size());
}

// Task Type
typedef enum {
	MI, TC, GAP_L1, GAP_LINF
} SCORING_CRITERION;

struct ClusterSet{
protected:
	/**
	 * Typedef of the vector backing this cycle set
	 */
	typedef vector<Cluster*> ClusterList;

	/**
	 * The set of clusters in this set
	 */
	ClusterList m_clusters;

	/**
	 * The set of cluster edges
	 */
	vector<ClusterEdge*> m_clusterEdges;

	/**
	 * Minimium score when computing the GAP or MI measures
	 */
	const static double MIN_SCORE = -10000.0;

	/**
	 * Computes the normalized, mutual information between cluster a and cluster b
	 */
	double computeMI(Cluster* a, Cluster* b) {
		// Boundary condition...the clusters have empty intersection
		vector<Variable*> isect;
		getIntersectionVar(a->getVariables(), b->getVariables(), isect);
		if (isect.empty()) {
			return 0;
		}

		// Get the union of the variables in these two clusters
		vector<Variable*> variables;
		getUnionVar(a->getVariables(), b->getVariables(), variables);

		// Get the union of the functions in these two clusters
		vector<LogFunction*> functions;
		for (int i = 0; i < (int)a->getFunctions().size(); i++) {
			functions.push_back(a->getFunctions()[i]);
		}
		for (int i = 0; i < (int)b->getFunctions().size(); i++) {
			functions.push_back(b->getFunctions()[i]);
		}

		// Determine the joint
		vector<Variable*> elim;
		LogFunction* belief_a_cup_b = new LogFunction();
		belief_a_cup_b->multiplyAndMarginalizeOut(functions, elim);
		belief_a_cup_b->normalize();

		// Determine cluster a
		getDifferenceVar(variables, a->getVariables(), elim);
		LogFunction* belief_a = belief_a_cup_b->marginalize(elim);
		belief_a->normalize();

		// Determine cluster b
		getDifferenceVar(variables, b->getVariables(), elim);
		LogFunction* belief_b = belief_a_cup_b->marginalize(elim);
		belief_b->normalize();

		// Compute the Mutual Information (MI)
		double MI = belief_a->getEntropy() + belief_b->getEntropy() - belief_a_cup_b->getEntropy();
		assert(MI >= 0);

		// Clean-up
		delete belief_a;
		delete belief_b;
		delete belief_a_cup_b;

		return MI;
	}

	/**
	 * Computes the total correlation between the variables in cluster a and cluster b
	 */
	double computeTC(Cluster* a, Cluster* b) {
		// Boundary condition...the clusters have empty intersection
		vector<Variable*> isect;
		getIntersectionVar(a->getVariables(), b->getVariables(), isect);
		if (isect.empty()) {
			return 0;
		}

		// Get the union of the variables in these two clusters
		vector<Variable*> variables;
		getUnionVar(a->getVariables(), b->getVariables(), variables);

		// Get the union of the functions in these two clusters
		vector<LogFunction*> functions;
		for (int i = 0; i < (int)a->getFunctions().size(); i++) {
			functions.push_back(a->getFunctions()[i]);
		}
		for (int i = 0; i < (int)b->getFunctions().size(); i++) {
			functions.push_back(b->getFunctions()[i]);
		}

		// Determine the joint
		vector<Variable*> elim;
		LogFunction* belief_a_cup_b = new LogFunction();
		belief_a_cup_b->multiplyAndMarginalizeOut(functions, elim);
		belief_a_cup_b->normalize();

		// Get the joint entropy
		double TC = -belief_a_cup_b->getEntropy();

		// Get the marginal entropy of every variable
		for (int i = 0; i < (int)variables.size(); i++) {
			elim.clear();
			removeVar(variables,variables[i],elim);
			LogFunction* marg_i = belief_a_cup_b->marginalize(elim);
			marg_i->normalize();
			TC += marg_i->getEntropy();
			delete marg_i;
		}

		// Clean-up
		delete belief_a_cup_b;

		return TC;
	}

	/**
	 * Computes the gap between the projection of the union of clusters a and b
	 * onto cluster c
	 */
	double computeGap(Cluster* a, Cluster* b, Cluster* c, SCORING_CRITERION criteria) {
		assert(criteria != MI);

		// Get the union of the variables in these two clusters
		vector<Variable*> variables;
		getUnionVar(a->getVariables(), b->getVariables(), variables);

		// Boundary condition
		vector<Variable*> isect;
		getIntersectionVar(variables, c->getVariables(), isect);
		if (isect.empty()) {
			return 0;
		}

		// Get the union of the functions in these two clusters
		vector<LogFunction*> functions;
		for (int i = 0; i < (int)a->getFunctions().size(); i++) {
			functions.push_back(a->getFunctions()[i]);
		}
		for (int i = 0; i < (int)b->getFunctions().size(); i++) {
			functions.push_back(b->getFunctions()[i]);
		}

		// Compute the message sent from ab to c
		vector<Variable*> elim;
		getDifferenceVar(variables, c->getVariables(), elim);
		LogFunction* m_ab_to_c = new LogFunction();
		m_ab_to_c->multiplyAndMarginalizeOut(functions, elim);

		// Compute the message sent from a to c
		elim.clear();
		getDifferenceVar(a->getVariables(), c->getVariables(), elim);
		LogFunction* m_a_to_c = new LogFunction();
		m_a_to_c->multiplyAndMarginalizeOut(a->getFunctions(), elim);

		// Compute the message sent from b to c
		elim.clear();
		getDifferenceVar(b->getVariables(), c->getVariables(), elim);
		LogFunction* m_b_to_c = new LogFunction();
		m_b_to_c->multiplyAndMarginalizeOut(b->getFunctions(), elim);

		// Compute the belief at c given message sent from ab to c
		elim.clear();
		LogFunction* c_given_ab = new LogFunction();
		c_given_ab->multiplyAndMarginalizeOut(c->getFunctions(), elim);
		c_given_ab->multiplyBy(m_ab_to_c);
		c_given_ab->normalize();

		// Compute the belief at c given messages sent from a to c and from b to c
		elim.clear();
		LogFunction* c_given_a_b = new LogFunction();
		c_given_a_b->multiplyAndMarginalizeOut(c->getFunctions(), elim);
		c_given_a_b->multiplyBy(m_a_to_c);
		c_given_a_b->multiplyBy(m_b_to_c);
		c_given_a_b->normalize();

		// Compare the difference between sending information independently
		// from a->c and b->c versus sending joint information from {a \cup b}->c
		double diff = 0;

		if (criteria == GAP_L1) {
			diff = LogFunction::getDistanceL1(c_given_ab, c_given_a_b);
		}
		else {
			diff = LogFunction::getDistanceLinf(c_given_ab, c_given_a_b);
		}

		// Clean-up
		delete m_ab_to_c;
		delete m_a_to_c;
		delete m_b_to_c;
		delete c_given_ab;
		delete c_given_a_b;

		return diff;
	}

	/**
	 * Removes the specified cluster from this cluster set, as well as all
	 * edges to/from this cluster
	 */
	inline bool removeCluster(Cluster* a, bool isRemoveEdges = true) {
		int idx = -1;
		for (int i = 0; i < (int)m_clusters.size(); i++) {
			if (m_clusters[i]->getID() == a->getID()) {
				idx = i;
				break;
			}
		}
		if (idx >= 0) {
			m_clusters.erase(m_clusters.begin() + idx);

			if (isRemoveEdges) {

				// Remove each edge involving cluster a
				for (int i = 0; i < (int)a->getEdges().size(); i++) {
					idx = -1;
					for (int j = 0; j < (int)m_clusterEdges.size(); j++) {
						if (m_clusterEdges[j]->getID() == a->getEdges()[i]->getID()) {
							idx = j;
							break;
						}
					}
					assert(idx >= 0);

					ClusterEdge* clusterEdge = m_clusterEdges[idx];
					m_clusterEdges.erase(m_clusterEdges.begin() + idx);
					delete clusterEdge;
				}
			}

			return true;
		}
		return false;
	}

	void clear() {
		for (int i = 0; i < (int)m_clusters.size(); i++) {
			delete m_clusters[i];
		}
		m_clusters.clear();
		for (int i = 0; i < (int)m_clusterEdges.size(); i++) {
			delete m_clusterEdges[i];
		}
		m_clusterEdges.clear();
	}

public:

	/**
	 * Creates an empty cluster set
	 */
	ClusterSet() { }

	/**
	 * Creates a cluster from the specified set of functions
	 */
	ClusterSet(const vector<LogFunction*>& functions) {

		// Sort the set of functions by size
		vector<LogFunction*> functions_sorted = vector<LogFunction*>(functions);
		sort(functions_sorted.begin(), functions_sorted.end(), LogFunction::compareBySize);

		// Create a cluster for each non-trivial function
		for	(int i = 0; i < (int)functions_sorted.size(); i++) {

			// Don't create clusters for empty factors
			if (functions_sorted[i]->getNumValues() > 1) {

				// Don't create a cluster for functions subsumed by other functions
				bool addedToExistingFactor = false;
				for (int j = 0; j < (int)m_clusters.size(); j++) {
					if (isSubsetVar(functions_sorted[i]->getScope(),
							        m_clusters[j]->getVariables())) {
						addedToExistingFactor = true;
						m_clusters[j]->addFunction(functions_sorted[i]);
						break;
					}
				}
				if (addedToExistingFactor) {
					continue;
				}

				// Create a new cluster and add it to the cluster set
				Cluster* cluster = new Cluster(functions_sorted[i]);
				m_clusters.push_back(cluster);
			}
		}

		// Determine the neighbors of each cluster, where a neighboring
		// cluster is a cluster with non-empty intersection
		for (int i = 0; i < (int)m_clusters.size(); i++) {
			for (int j = i + 1; j < (int)m_clusters.size(); j++) {
				vector<Variable*> isect;
				getIntersectionVar(m_clusters[i]->getVariables(), m_clusters[j]->getVariables(), isect);
				if (!isect.empty()) {
					ClusterEdge* clusterEdge = new ClusterEdge(m_clusters[i], m_clusters[j]);
					m_clusterEdges.push_back(clusterEdge);
				}
			}
		}
	}

	~ClusterSet() {
		for (int i = 0; i < (int)m_clusterEdges.size(); i++) {
			delete m_clusterEdges[i];
		}
		m_clusterEdges.clear();

		for (int i = 0; i < (int)m_clusters.size(); i++) {
			delete m_clusters[i];
		}
		m_clusters.clear();
	}

	/**
	 * Computes the size of the cluster produced by merging clusters u and v
	 */
	static inline int sizeOfMerge(Cluster* u, Cluster* v) {
		// Get the union of the variables in these two clusters
		vector<Variable*> variables;
		getUnionVar(u->getVariables(), v->getVariables(), variables);

		return (int)variables.size();
	}

	/**
	 * Computes the size of the cluster produced by merging clusters on this edge
	 */
	static inline int sizeOfMerge(ClusterEdge* edge) {
		vector<Variable*> variables;
		getUnionVar(edge->getCluster_i()->getVariables(), edge->getCluster_j()->getVariables(), variables);

		return (int)variables.size();
	}

	/**
	 * Merges the two clusters on this edge.
	 */
	void merge(ClusterEdge* clusterEdge) {
		// Locate this edge
		bool foundEdge = false;
		for (int i = 0; i < (int)m_clusterEdges.size(); i++) {
			if (m_clusterEdges[i]->getID() == clusterEdge->getID()) {
				foundEdge = true;
				break;
			}
		}
		assert(foundEdge);

		Cluster* cluster_i = clusterEdge->getCluster_i();
		Cluster* cluster_j = clusterEdge->getCluster_j();
		cluster_i->merge(cluster_j, true);
		removeCluster(cluster_j);
		delete cluster_j;

		// If this merge yielded clusters that are subsumed by cluster_i,
		// continue to apply the merge operation. For example, if
		// cluster_i = {1,2,3} and cluster_j = {2,4}, the merge will
		// produce cluster_i' = {1,2,3,4}. Now if there is some cluster
		// cluster_k = {3,4} then it should be merged with cluster_i' as
		// well
		bool isMerged = true;
		while (isMerged) {
			isMerged = false;

			// Check all neighbors of cluster_i
			for (int i = 0; i < (int)cluster_i->getEdges().size(); i++) {
				Cluster* cluster_k = cluster_i->getEdges()[i]->getOtherCluster(cluster_i);

				if (isSubsetVar(cluster_k->getVariables(), cluster_i->getVariables())) {
					isMerged = true;

					cluster_i->merge(cluster_k, true);
					removeCluster(cluster_k);
					delete cluster_k;
					break;
				}
			}
		}
	}

	/**
	 * Uses the mini-bucket schematic to construct a set of clusters
	 * such that each cluster has at most iBound variables
	 */
	void runMiniBucket(vector<int> order, vector<Variable*> variables, int iBound) {

		// Determine the size of the bucket tree
		int numVariables = variables.size();
		assert(variables.size() == order.size());

		//
		// The mini-bucket tree structure
		//
		// The bucket tree has one entry for each variable in the ordering. Each
		// bucket can contain a set of clusters, where each cluster is a set of variables.
		vector<vector<Cluster*> > clustersInBucket = vector<vector<Cluster*> >(numVariables);

		//
		// Assign current set of clusters to buckets
		//
		// Each bucket is associated with a variable and buckets are ordered
		// along the specified elimination order. Place each cluster i
		// in the first bucket j such that variable j occurs in cluster i.
		for (int i = 0; i < (int)m_clusters.size(); i++) {
			bool assigned = false;
			for (int j = 0; j < numVariables; j++) {
				if (containsVar(m_clusters[i]->getVariables(), variables[order[j]])) {
					clustersInBucket[j].push_back(m_clusters[i]);
					assigned = true;
					break;
				}
			}
			assert(assigned);
		}

		// Clear the cluster set, since we will be augmenting it via
		// the mini-bucket schematic
		m_clusters.clear();

		// At each bucket i, partition the original clusters assigned to that
		// bucket and any messages from previous buckets in the ordering into
		// clusters of at most iBound variables.  Perform this partitioning in
		// a greedy fashion based upon number of variables in a message. Then
		// simulate passing of messages from each bucket to some future bucket.
		for (int i = 0; i < numVariables; i++) {

			// Sort the clusters in bucket i largest to smallest by scope size
			sort(clustersInBucket[i].begin(), clustersInBucket[i].end(),compareClustersBySize);

			//
			// Begin the partitioning process. Greedily merge functions
			// while we remain below the specified iBound
			//
			bool foundMerge = true;
			while (foundMerge) {
				foundMerge = false;
				for (int j = 0; j < (int)clustersInBucket[i].size(); j++) {
					for (int k = j + 1; k < (int)clustersInBucket[i].size(); k++) {

						// If merging clusters j and k in bucket i does not violate the iBound
						if (sizeOfMerge(clustersInBucket[i][j], clustersInBucket[i][k]) <= iBound) {
							clustersInBucket[i][j]->merge(clustersInBucket[i][k],false);
							clustersInBucket[i].erase(clustersInBucket[i].begin() + k);
							foundMerge = true;
							break;
						}
					}
					if (foundMerge) {
						break;
					}
				}
			}

			// Now perform the forwarding operation. A message from each cluster is generated
			// by removing the variable in bucket i.  This message is forwarded to the first
			// bucket j that appears in this message scope
			for (int j = 0; j < (int)clustersInBucket[i].size(); j++) {

				// Add each cluster in bucket i to the cluster set
				m_clusters.push_back(clustersInBucket[i][j]);

				// Create the forwarded cluster
				vector<Variable*> clusterScope;
				removeVar(clustersInBucket[i][j]->getVariables(), variables[order[i]], clusterScope);

				// Don't propagate empty messages
				if (!clusterScope.empty()) {
					Cluster* newCluster = new Cluster(clusterScope);

					// Determine which bucket this cluster should be placed in
					bool forwardedCluster = false;
					for (int k = i + 1; k < numVariables; k++) {
						if (containsVar(clusterScope, variables[order[k]])) {
							clustersInBucket[k].push_back(newCluster);
							forwardedCluster = true;
							break;
						}
					}
					assert(forwardedCluster);
				}

			} // for each cluster in bucket i

		} // for each bucket (i = 0; i < numVariables)

		// Now simplify the cluster set by removing redundant and/or subsumed clusters
		sort(m_clusters.begin(), m_clusters.end(), compareClustersBySize);
		for (int i = (int)m_clusters.size() - 1; i >= 0; i--) {
			Cluster* cluster_i = m_clusters[i];
			for (int j = i - 1; j >= 0; j--) {

				// If cluster i is subsumed by cluster j, merge cluster i
				// into cluster j and remove cluster i
				if (isSubsetVar(cluster_i->getVariables(), m_clusters[j]->getVariables())) {
					m_clusters[j]->merge(cluster_i,false);
					removeCluster(cluster_i,false);
					delete cluster_i;
					break;
				}
			}
		}
	}

	/**
	 * Outputs the cluster set
	 */
	string toString() {
		stringstream ss;
		for (int i = 0; i < (int)m_clusters.size(); i++) {
			ss << i << ":" << varsToString(m_clusters[i]->getVariables()) << endl;
		}
		return ss.str();
	}

	/**
	 * Gets the set of variables corresponding to this cluster set
	 */
	vector<vector<Variable*> > getClusterVars() {
		vector<vector<Variable*> > clusterVars;
		for (int i = 0; i < (int)m_clusters.size(); i++) {
			clusterVars.push_back(m_clusters[i]->getVariables());
		}
		return clusterVars;
	}

	/**
	 * Gets the set of variables corresponding to this cluster set,
	 * given that the clusters on this edge have been merged
	 */
	vector<vector<Variable*> > getClusterVars(ClusterEdge* edge) {
		vector<vector<Variable*> > clusterVars;
		Cluster* a = edge->getCluster_i();
		Cluster* b = edge->getCluster_j();
		bool foundA = false;
		bool foundB = false;
		for (int i = 0; i < (int)m_clusters.size(); i++) {
			if (m_clusters[i]->getID() == a->getID()) {
				foundA = true;
			}
			else if (m_clusters[i]->getID() == b->getID()) {
				foundB = true;
			}
			else {
				clusterVars.push_back(m_clusters[i]->getVariables());
			}
		}
		assert(foundA && foundB);

		vector<Variable*> merged;
		getUnionVar(a->getVariables(), b->getVariables(), merged);

		// Remove any clusters subsumed by this merge
		for (int j = (int)clusterVars.size() - 1; j >= 0; j--) {
			if (isSubsetVar(clusterVars[j], merged)) {
				clusterVars.erase(clusterVars.begin() + j);
			}
		}
		clusterVars.push_back(merged);

		return clusterVars;
	}

	/**
	 * Gets the current set of cluster edges
	 */
	const vector<ClusterEdge*>& getClusterEdges() {
		return m_clusterEdges;
	}

	/**
	 * Iterator Functions
	 */
	typedef ClusterList::iterator iterator;
	typedef ClusterList::const_iterator const_iterator;
	iterator begin() {
		return m_clusters.begin();
	}
	iterator end() {
		return m_clusters.end();
	}

};


#endif /* PARTITION_H_ */
