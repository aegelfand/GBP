/*
 * RegionGraph.cpp
 *
 *  Created on: Jul 26, 2011
 *      Author: agelfand
 */

#include <queue>

#include "RegionGraph.h"

int Region::s_count = 0;
int Edge::s_count = 0;

/**
 * Initializes this region by:
 */
void Region::initialize() {

	// Clear the set of beliefs
	m_belief->clear();
	m_belief->multiplyBy(m_baseFactor);
	m_belief->normalize();

	// If this region has no incoming edges it will never be
	// updated, causing message passing to fail to converge
	if (m_edges.empty()) {
		m_differenceLinf = 0;
	}
}

/**
 * Initializes a region graph by:
 *   1) Determining the ancestors of each region; and
 *   2) Setting the counting number of each region;
 * NOTE: Assumes that <tt>regions</tt> is sorted, with entries placed in decreasing
 * order by size (i.e. number of variables in each region).
 */
void Region::determineAncestorsRG(vector<Region*>& regions) {

	// Determine the ancestors of each region
	for (int i = 0; i < (int)regions.size(); i++) {
		int j = 0;
		vector<Region*> ancestors;
		while (j < i) {
			if (isSubsetVar(regions[i]->getVariables(),regions[j]->getVariables())) {
				ancestors.push_back(regions[j]);
			}
			j++;
		}
		std::sort(ancestors.begin(),ancestors.end(),compareRegionByID);
		regions[i]->m_ancestors = ancestors;
	}

	// Determine the counting number and covering number of each region
	for (int i = 0; i < (int)regions.size(); i++) {

		for (int j = 0; j < (int)regions[i]->getAncestors().size(); j++) {

			if (regions[i]->getAncestors()[j]->getID() == regions[i]->getID()) {
				continue;
			}

			regions[i]->m_countingNumber += regions[i]->getAncestors()[j]->getCountingNumber();
		}

		regions[i]->m_countingNumber = 1 - regions[i]->m_countingNumber;
	}
}

/**
 * Adds an edge to this node
 */
void Region::addEdge(Edge* edge) {
	// Ensure this edge is either a parent or child of this region
	bool isChild = (edge->getChild() == this);
	bool isParent = (edge->getParent() == this);
	assert(isChild || isParent);

	// Validate the containment condition
	if (isChild) {
		assert(isSubsetVar(m_variables, edge->getParent()->m_variables));
	}
	else {
		assert(isSubsetVar(edge->getChild()->m_variables, m_variables));
	}

	m_edges.push_back(edge);
}

/**
* Removes an edge to/from this region
*/
void Region::removeEdge(Edge* edge) {
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
 * Removes the specified region from the set of ancestors
 */
void Region::removeAncestor(Region* region) {
	int idx = -1;
	for (int i = 0; i < (int)m_ancestors.size(); i++) {
		if (m_ancestors[i]->getID() == region->getID()) {
			idx = i;
			break;
		}
	}
	if (idx >= 0) {
		m_ancestors.erase(m_ancestors.begin() + idx);
	}
}

/**
 * Gets the entropy at this region
 */
double Region::getEntropy() {
	return m_belief->getEntropy();
}

/**
 * Gets the energy at this region
 */
double Region::getEnergy() {
	double energy = 0;

	if (!m_assignedFunctions.empty()) {
		for (int i = 0; i < m_belief->getNumValues(); i++) {
			Variable::setAddress(m_belief->getScope(),i);

			double partialEnergy = 0;
			for (int k = 0; k < (int)m_assignedFunctions.size(); k++) {
				int idx = Variable::getAddress(m_assignedFunctions[k]->getScope());
				partialEnergy += m_assignedFunctions[k]->getLogValueAt(idx);
			}

			energy += m_belief->getValueAt(i)*partialEnergy;
		}
	}

	return energy;
}

/**
 * Updates the belief at this region
 */
void Region::updateBelief(double alpha) {

	// If this is an outer region
	if (m_ancestors.empty()) {
		vector<LogFunction*> funcs;
		for (int i = 0; i < (int)m_edges.size(); i++) {
			funcs.push_back(m_edges[i]->getChildToParentMsg());
		}

		// Compute the new belief
		LogFunction* temp = new LogFunction(*m_baseFactor);
		vector<Variable*> emptySet;
		temp->multiplyAndMarginalizeOut(funcs,emptySet);
		temp->normalize();

		// Compute the damped update
		temp->weightedAdd(m_belief, 1.0 - alpha);
		//temp->weightedMultiply(m_belief,1.0 - alpha);
		//temp->normalize();

		// Check convergence
		m_differenceLinf = LogFunction::getDistanceLinf(temp,m_belief);

		// Perform the update and clean up
		m_belief->assign(*temp);
		delete temp;
	}
	// This is an inner region
	else {

		vector<LogFunction*> funcs;
		for (int i = 0; i < (int)m_edges.size(); i++) {
			funcs.push_back(m_edges[i]->getParentToChildMsg());
		}

		// Compute the new belief
		LogFunction* temp = new LogFunction();
		vector<Variable*> emptySet;
		temp->multiplyAndMarginalizeOut(funcs,emptySet);
		temp->exponentiateBy(1.0/(m_countingNumber + ((int)m_edges.size())));
		temp->normalize();

		// Compute the damped update
		temp->weightedAdd(m_belief, 1.0 - alpha);
		//temp->weightedMultiply(m_belief,1.0 - alpha);
		//temp->normalize();

		// Check convergence
		m_differenceLinf = LogFunction::getDistanceLinf(temp,m_belief);

		// Perform the update and clean up
		m_belief->assign(*temp);
		delete temp;
	}
}

/**
 * Updates the base factor at this region, if it is an outer
 * region
 */
void Region::updateBaseFactor() {
	// If this is an outer region
	if (m_ancestors.empty()) {
		// Clear the current base factor
		m_baseFactor->clear();

		// Multiply in the original functions
		for (int i = 0; i < (int)m_assignedFunctions.size(); i++) {
			m_baseFactor->multiplyBy(m_assignedFunctions[i]);
		}

		// Multiply in the exponentiated belief of each child region
		for (int i = 0; i < (int)m_edges.size(); i++) {
			Region* child = m_edges[i]->getChild();

			if (child->getCountingNumber() != 0) {
				LogFunction* temp = new LogFunction(*child->getBelief());
				double exponent = ((double)child->getCountingNumber()) / ((int)child->getEdges().size());
				temp->exponentiateBy(-exponent);
				m_baseFactor->multiplyBy(temp);
				delete temp;
			}
		}
	}
}
