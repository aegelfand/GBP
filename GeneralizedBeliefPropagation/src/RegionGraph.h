/*
 * RegionGraph.h
 *
 *  Created on: Jul 26, 2011
 *      Author: agelfand
 */

#ifndef REGIONGRAPH_H_
#define REGIONGRAPH_H_

#include <algorithm>
#include <sstream>
#include <vector>

#include "Graph.h"
#include "LogFunction.h"
#include "myRandom.h"
#include "Util.h"
#include "Variable.h"

using namespace std;

struct Edge;

struct Region {
protected:
	/**
	 * The set of variables in this region (in ascending order by variable ID)
	 */
	vector<Variable*> m_variables;

	/**
	 * The set of functions assigned to this region.
	 */
	vector<LogFunction*> m_assignedFunctions;

	/**
	 * The base factor of this region
	 */
	LogFunction* m_baseFactor;

	/**
	 * The current belief at this region
	 */
	LogFunction* m_belief;

	/**
	 * The set of edges to/from this region.
	 */
	vector<Edge*> m_edges;

    /**
     * ID of this region
     */
    int m_ID;

    /**
     * Count of regions
     */
    static int s_count;

    /**
     * The counting number for this region.
     */
    int m_countingNumber;

    /**
     * The regions that are ancestors of this region
     */
    vector<Region*> m_ancestors;

    /**
     * Linf distance between current belief and previous belief
     */
    double m_differenceLinf;

public:

    /**
     * Gets the ID that will be assigned to the next Region
     */
    static void resetIDs() {
    	s_count = 0;
    }

    /**
     * Compares regions a and b based upon the number of variables
     * in each region. Used to place in descending order by size.
     */
    static inline bool compareRegionBySize_Descend(Region* a, Region* b) {
    	return (((int)a->getVariables().size()) > ((int)b->getVariables().size()));
    }

    /**
	 * Compares regions a and b based upon the number of variables
	 * in each region. Used to place in descending order by size.
	 */
	static inline bool compareRegionBySize_Ascend(Region* a, Region* b) {
		return (((int)a->getVariables().size()) < ((int)b->getVariables().size()));
	}

    /**
     * Compares regions a and b based upon their IDs.
     */
    static inline bool compareRegionByID(Region* a, Region* b) {
    	return (a->getID() < b->getID());
    }

    /**
	 * Gets the difference of sets A and B and returns the result in C - i.e. C = A - B
	 * Assumes that both <tt>A</tt> and <tt>B</tt> are sorted.
	 */
	static inline void getDifferenceRegion(const vector<Region*>& A, const vector<Region*>& B, vector<Region*>& C) {
		vector<Region*> temp;
		temp.resize(A.size());
		vector<Region*>::iterator curr,last;
		last = temp.end();
		curr = set_difference(A.begin(), A.end(), B.begin(), B.end(), temp.begin(), compareRegionByID);
		temp.erase(curr,last);
		C = temp;
	}

	/**
	 * Gets the intersection of sets A and B and returns the result in C - i.e. C = A cap B
	 * Assumes that both <tt>A</tt> and <tt>B</tt> are sorted.
	 * @return bool -- Returns <tt>true</tt> if the intersection is not empty; <tt>false</tt> otherwise;
	 */
	static inline bool getIntersectionRegion(const vector<Region*>& A, const vector<Region*>& B, vector<Region*>& C) {
		vector<Region*> temp;
		temp.resize(A.size());
		vector<Region*>::iterator curr,last;
		last = temp.end();
		curr = set_intersection(A.begin(), A.end(), B.begin(), B.end(), temp.begin(), Region::compareRegionByID);
		temp.erase(curr,last);
		C = temp;
		return !C.empty();
	}

	/**
	 * Gets the union of sets A and B and returns the result in C - i.e. C = A cup B
	 * Assumes that both <tt>A</tt> and <tt>B</tt> are sorted.
	 */
	static inline void getUnionRegion(const vector<Region*>& A, const vector<Region*>& B, vector<Region*>& C) {
		vector<Region*> temp;
		temp.resize(A.size() + B.size());
		vector<Region*>::iterator curr,last;
		last = temp.end();
		curr = set_union(A.begin(), A.end(), B.begin(), B.end(), temp.begin(), Region::compareRegionByID);
		temp.erase(curr,last);
		C = temp;
	}

	/**
	 * Adds the specified region B to the set A and returns the result in C
	 */
	static inline void addRegion(const vector<Region*>& A, Region* B, vector<Region*>& C) {
		vector<Region*> temp;
		temp.push_back(B);
		getUnionRegion(A,temp,C);
	}

	/**
	 * Removes the specified region B from the set A and returns the result in C
	 */
	static inline void removeRegion(const vector<Region*>& A, Region* B, vector<Region*>& C) {
		vector<Region*> temp;
		temp.push_back(B);
		getDifferenceRegion(A,temp,C);
	}

	/**
	 * Initializes this region by clearing the current region belief.
	 */
	void initialize();

    /**
     * Get ancestors of this region
     */
    const vector<Region*>& getAncestors() {
    	return m_ancestors;
    }

	/**
	 * Default constructor
	 * @param variables - The set of variables in this cluster.
	 *
	 * NOTE:  Assumes that the variables are sorted in ascending order by ID.
	 */
	Region(vector<Variable*> variables) {
		m_variables = variables;
		m_ID = s_count++;

		m_countingNumber = 0;
		m_differenceLinf = 1;
		m_belief = new LogFunction(m_variables);
		m_baseFactor = new LogFunction();
	}

	Region(Variable* var) {
		m_variables.push_back(var);
		m_ID = s_count++;

		m_countingNumber = 0;
		m_differenceLinf = 1;
		m_belief = new LogFunction(m_variables);
		m_baseFactor = new LogFunction();
	}

	~Region() {
		m_variables.clear();
		m_ancestors.clear();
		m_assignedFunctions.clear();
		m_edges.clear();
		delete m_belief;
		delete m_baseFactor;
	}

	/**
	 * Get the counting number of this region
	 */
	int getCountingNumber() {
		return m_countingNumber;
	}

	/**
	 * Sets the counting number for this region
	 */
	void setCountingNumber(int countingNumber) {
		m_countingNumber = countingNumber;
	}

	/**
	 * Gets the Linf difference between the current and previous beliefs
	 */
	double getDifferenceLinf() {
		return m_differenceLinf;
	}

    /**
     * Gets the variables in this region
     */
    const vector<Variable*>& getVariables() {
        return m_variables;
    }

    /**
     * Gets the set of edges to/from this region
     */
    const vector<Edge*>& getEdges() {
    	return m_edges;
    }

	/**
	 * Override of the equals operator. Comparison is based on region ID.
	 */
	bool operator==(const Region& other) { return (m_ID == other.m_ID); }

	/**
	 * Override of the < operator.  Comparison based on region ID.
	 */
	bool operator<(const Region& other) const {
		return (m_ID < other.m_ID);
	};

    /**
     * Assigns a function to this region.  All function additions should be performed
     * prior to message updating.
     */
    void addFunction(LogFunction* function) {
    	// Ensure function scope is a subseteq m_variables
    	vector<Variable*> isect;
    	getIntersectionVar(m_variables, function->getScope(), isect);
    	assert(isect.size() == function->getScope().size());

    	m_assignedFunctions.push_back(function);
    	m_baseFactor->multiplyBy(function);
    }

    /**
     * Gets the current belief at this region
     */
    LogFunction* getBelief() {
    	return m_belief;
    }

    /**
     * Updates the belief at this region
     */
    void updateBelief(double alpha);

    /**
     * Updates the base factor at this region, if it is an outer
     * region
     */
    void updateBaseFactor();

    /**
     * Gets the entropy at this region
     */
    double getEntropy();

    /**
     * Gets the energy at this region
     */
    double getEnergy();

    /**
     * Gets the set of functions assigned to this region.
     */
    const vector<LogFunction*>& getAssignedFunctions() {
    	return m_assignedFunctions;
    }

    /**
     * Adds an edge to/from this region
     */
    void addEdge(Edge* edge);

    /**
     * Removes an edge to/from this region
     */
    void removeEdge(Edge* edge);

    /**
     * Removes the specified region from the set of ancestors
     */
    void removeAncestor(Region* region);

    /**
     * Gets the ID of this region
     */
    int getID() {
    	return m_ID;
    }

    /**
     * Returns a string representation of the region
     */
    string toString() {
    	std::stringstream ss;
    	ss << varsToString(m_variables);
    	return ss.str();
    }

	/**
	 * Initializes a region graph by:
	 *   1) Determining the ancestors of each region; and
	 *   2) Setting the counting number of each region.
	 * NOTE: Assumes that <tt>regions</tt> is sorted, with entries placed in decreasing
	 * order by size (i.e. number of variables in each region).
	 */
	static void determineAncestorsRG(vector<Region*>& regions);

};

struct Edge {
protected:
    /**
     * Parent region
     */
    Region* m_parent;

    /**
     * Child region
     */
    Region* m_child;

	/**
	 * The set of variables eliminated from the clique of the parent
	 */
	vector<Variable*> m_elim;

    /**
     * ID of this edge
     */
    int m_ID;

    /**
     * Count of edges
     */
    static int s_count;

    /**
     * The parent-child msg
     */
    LogFunction* m_parentToChildMsg;

    /**
     * The child-parent msg
     */
    LogFunction* m_childToParentMsg;

public:

    /**
	 * Default constructor
	 * @param Region* parent --
	 * @param Region* child --
	 */
	Edge(Region* parent, Region* child) {
		m_parent = parent;
		m_child = child;
		m_parent->addEdge(this);
		m_child->addEdge(this);
		m_ID = s_count++;

		m_parentToChildMsg = NULL;
		m_childToParentMsg = NULL;

		// Determine the eliminator
		getDifferenceVar(m_parent->getVariables(),m_child->getVariables(),m_elim);
	}

	~Edge() {
		m_elim.clear();
		delete m_parentToChildMsg;
		delete m_childToParentMsg;
	}

	/**
	 * Gets the ID that will be assigned to the next Region
	 */
	static void resetIDs() {
		s_count = 0;
	}

	/**
	 * Creates the parent-child and child-parent messages
	 */
	void initialize() {
		if (m_parentToChildMsg == NULL) {
			m_parentToChildMsg = new LogFunction(m_child->getVariables());
		}
		else {
			m_parentToChildMsg->clear();
		}
		if (m_childToParentMsg == NULL) {
			m_childToParentMsg = new LogFunction(m_child->getVariables());
		}
		else {
			m_childToParentMsg->clear();
		}
	}

	/**
	 * Override of the equals operator. Comparison is based on edge ID.
	 */
	bool operator==(const Edge& other) { return (m_ID == other.m_ID); }
	bool operator !=(const Edge& other) { return (m_ID != other.m_ID); }

    /**
     * Compares edges a and b based upon their IDs.
     */
    static inline bool compareEdgeByID(Edge* a, Edge* b) {
    	return (a->getID() < b->getID());
    }

	/**
	 * Gets the id of this edge
	 */
	int getID() {
		return m_ID;
	}

	/**
	 * Gets the parent region of this edge
	 */
	Region* getParent() {
		return m_parent;
	}

	/**
	 * Gets the child region of this edge
	 */
	Region* getChild() {
		return m_child;
	}

	/**
	 * Gets the eliminator on this edge
	 */
	const vector<Variable*>& getElim() {
		return m_elim;
	}

	/**
	 * Displays a string representation of the edge
	 */
	string toString() {
		std::stringstream ss;
		ss << m_parent->toString().c_str() << "-" << m_child->toString().c_str();
		return ss.str();
	}

	/**
	 * Gets the parent to child message
	 */
	LogFunction* getParentToChildMsg() {
		return m_parentToChildMsg;
	}

	/**
	 * Gets the child to parent message
	 */
	LogFunction* getChildToParentMsg() {
		return m_childToParentMsg;
	}

	/**
	 * Updates the parent-to-child message
	 */
	void updateParentToChildMsg() {
		m_parentToChildMsg->assign(*m_parent->getBelief());
		m_parentToChildMsg->marginalizeOut(m_elim);
		m_parentToChildMsg->divideBy(m_childToParentMsg);
	}

	/**
	 * Updates the child-to-parent message
	 */
	void updateChildToParentMsg() {
		m_childToParentMsg->assign(*m_child->getBelief());
		m_childToParentMsg->divideBy(m_parentToChildMsg);
	}

};

#endif /* REGIONGRAPH_H_ */
