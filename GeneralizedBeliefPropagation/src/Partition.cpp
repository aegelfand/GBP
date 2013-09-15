/*
 * Partition.cpp
 *
 *  Created on: Feb 6, 2013
 *      Author: root
 */

#include "Partition.h"

int Cluster::s_count = 0;
int ClusterEdge::s_count = 0;


string Cluster::toString() {
	stringstream ss;
	ss << m_ID << ":(" << varsToString(m_variables) << "):";
	for (int i = 0; i < (int)m_edges.size(); i++) {
		if (i > 0) {
			ss << ",";
		}
		ss << m_edges[i]->getOtherCluster(this)->getID();
	}

	return ss.str();
}

/**
 * Adds an edge to this cluster
 */
void Cluster::addEdge(ClusterEdge* edge) {
	m_edges.push_back(edge);
}

/**
 * Removes a neighbor to this cluster
 */
void Cluster::removeEdge(ClusterEdge* edge) {
	int idx = -1;
	for (int i = 0; i < (int)m_edges.size(); i++) {
		if (m_edges[i]->getID() == edge->getID()) {
			idx = i;
			break;
		}
	}
	if (idx >= 0) {
		m_edges.erase(m_edges.begin() + idx);
	}
}

/**
 * Gets the neighbors of this cluster whose cluster
 * size is greater than the specified iBound
 */
vector<ClusterEdge*> Cluster::getEdges(int iBound) {
	assert(iBound > 0);
	vector<ClusterEdge*> edges;
	for (int i = 0; i < (int)m_edges.size(); i++) {
		if (((int)m_edges[i]->getOtherCluster(this)->getVariables().size()) > iBound) {
			edges.push_back(m_edges[i]);
		}
	}
	return edges;
}

/**
 * Returns true if cluster u is a neighbor of this cluster
 */
bool Cluster::isNeighbor(Cluster* u) {
	for (int i = 0; i < (int)m_edges.size(); i++) {
		if (m_edges[i]->getOtherCluster(this)->getID() == u->getID()) {
			return true;
		}
	}
	return false;
}

/**
 * Merges the specified cluster into this cluster.
 * If isMergeEdges is <tt>true</tt> the edges of
 * this cluster and other cluster are merged too.
 */
void Cluster::merge(Cluster* other, bool isMergeEdges) {
	getUnionVar(m_variables, other->m_variables, m_variables);
	for (int i = 0; i < (int)other->m_functions.size(); i++) {
		m_functions.push_back(other->m_functions[i]);
	}

	if (!isMergeEdges) {
		return;
	}

	// Remove the edge between this cluster and other
	int idx = -1;
	for (int i = 0; i < (int)m_edges.size(); i++) {
		if (m_edges[i]->containsCluster(other)) {
			idx = i;
			break;
		}
	}
	assert(idx >= 0);
	m_edges.erase(m_edges.begin() + idx);

	// Combine the edges of cluster other with this cluster. For each
	// neighbor u of other:
	//  1) If u is also a neighbor of this cluster, remove the edge
	//     between u and other
	//  2) If u is not a neighbor of this cluster, add an edge between
	//     u and this cluster
	for (int i = (int)other->m_edges.size() - 1; i >= 0; i--) {

		// If u is cluster this
		Cluster* u = other->m_edges[i]->getOtherCluster(other);
		if (u->getID() == this->getID()) {
			continue;
		}
		// If u is a neighbor of this cluster
		else if (isNeighbor(u)) {
			// Remove edge (u,other) from cluster u
			idx = -1;
			for (int j = 0; j < (int)u->m_edges.size(); j++) {
				if (u->m_edges[j]->getOtherCluster(u)->getID() == other->getID()) {
					idx = j;
					break;
				}
			}
			assert(idx >= 0);
//cout << "  Removing edge " << u->m_edges[idx]->toString() << " from common neighbor " << u->toString() << endl;
			u->removeEdge(u->m_edges[idx]);
		}
		// Otherwise u is not a neighbor of this cluster
		else {
			// Locate the edge (u,other)
			idx = -1;
			for (int j = 0; j < (int)u->m_edges.size(); j++) {
				if (u->m_edges[j]->getOtherCluster(u)->getID() == other->getID()) {
					idx = j;
					break;
				}
			}
			assert(idx >= 0);

			// Swap other with this
//cout << "  Swapping cluster " << other->toString() << " on edge " << u->m_edges[idx]->toString() << " for cluster " << this->toString() << endl;
			u->m_edges[idx]->swap(other, this);
		}
	}

	// Mark all edges from this cluster as needing an update
	for (int i = 0; i < (int)m_edges.size(); i++) {
		m_edges[i]->setUpdateStatus(true);
	}

}
