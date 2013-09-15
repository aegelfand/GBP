/*
 * JoinGraph.cpp
 *
 *  Created on: Sep 26, 2011
 *      Author: agelfand
 */

#include "JoinGraph.h"

int JGNode::s_count = 0;
int JGEdge::s_count = 0;

/**
 * Adds an edge to this node
 */
void JGNode::addEdge(JGEdge* edge) {
	// Ensure this edge is either a parent or child of this region
	bool is_i = (edge->getNode_i() == this);
	bool is_j = (edge->getNode_j() == this);
	assert((is_i || is_j) && !(is_i && is_j));
	m_edges.push_back(edge);
}

/**
 * Removes the specified edge from this node
 */
void JGNode::removeEdge(JGEdge* edge) {
	int idx = -1;
	for (int i = 0; i < (int)m_edges.size(); i++) {
		if (m_edges[i] == edge) {
			idx = i;
			break;
		}
	}
	if (idx >= 0) {
		m_edges.erase(m_edges.begin() + idx);
	}
}
