/*
 * Variable.h
 *	Implementation of a discrete variable
 *  Created on: Sep 28, 2010
 *      Author: andrew.gelfand
 */

#ifndef VARIABLE_H_
#define VARIABLE_H_

#include <algorithm>
#include <vector>
#include <cassert>
#include <string>
#include <iostream>
#include <sstream>
#include <math.h>
#include <map>

using namespace std;

// Used to indicate that a variable has not been assigned a value
#define INVALID_VALUE -1

/**
 * A data structure for variables.  Variable domains are assumed to run from 0...n-1, where
 * n is the domain size of the variable.
 */
struct Variable {

protected:

	/**
	 * The unique ID of this variable
	 */
	int m_ID;

	/**
	 * The initial domain size.
	 */
	int m_initDomainSize;

	/**
	 * The current domain size.
	 */
	int m_domainSize;

	/**
	 * A dummy value assignment used for addressing purposes
	 */
	int m_addressValue;

	/**
	 * Flag indicating whether the domain has been reduced or not
	 */
	bool m_isReduced;

	/**
	 * Mapping from values in the current domain to values in the original domain. For example,
	 * let A be a ternary variable taking values {0,1,2} and say that the domain of A is to
	 * be reduced to just the values {0,2}. Then this mapping will be a vector of length 2,
	 * with the first element equal '0' and the second element equal '2'.
	 */
	vector<int> m_mapping;

	/**
	 * Gets the vector of valid values for this variable. Returns a boolean vector of length
	 * equal to the initial domain size. If valid[i]=true, then the ith state of this variable's
	 * initial domain is valid.
	 */
	vector<bool> getValid() {
		vector<bool> returnVal = vector<bool>(m_initDomainSize);

		if (m_initDomainSize == m_domainSize) {
			std::fill(returnVal.begin(),returnVal.end(),true);
		}
		else {
			std::fill(returnVal.begin(),returnVal.end(),false);
			for (int i = 0; i < (int)m_mapping.size(); i++) {
				returnVal[m_mapping[i]] = true;
			}
		}

		return returnVal;
	}

	/**
	 * Protected constructor. Creates a copy of the specified
	 * variable, but expands the domain to include all states
	 * in <tt>valid</tt>.
	 */
	Variable(Variable& var, vector<bool> valid) {
		m_ID = var.m_ID;
		m_initDomainSize = var.m_initDomainSize;

		assert(((int)valid.size()) == m_initDomainSize);

		m_addressValue = INVALID_VALUE;

		m_mapping = vector<int>();
		m_domainSize = 0;
		for (int i = 0; i < (int)valid.size(); i++) {
			if (valid[i]) {
				m_mapping.push_back(i);
				m_domainSize++;
			}
		}

		if (m_domainSize == m_initDomainSize) {
			m_isReduced = false;
			m_mapping.clear();
		}
		else {
			m_isReduced = true;
		}
	}

public:

	/**
	 * Creates a discrete variable with domain values [0,domainSize).
	 */
	Variable(const int& id, const int domainSize) {
		assert(domainSize > 0);

		m_ID = id;
		m_initDomainSize = domainSize;
		m_domainSize = domainSize;
		m_addressValue = INVALID_VALUE;
		m_isReduced = false;

		// If this is a unary variable
		if (m_domainSize == 1) {
			m_addressValue = 0;
			m_mapping = vector<int>(1);
			m_mapping[0] = 0;
		}
	}

	/**
	 * Copy constructor
	 */
	Variable(Variable& var) {
		m_ID = var.m_ID;
		m_initDomainSize = var.m_initDomainSize;
		m_domainSize = var.m_domainSize;
		m_addressValue = var.m_addressValue;
		m_isReduced = var.m_isReduced;

		if (m_isReduced) {
			m_mapping = vector<int>(var.m_mapping);
		}
	}

	/**
	 * Default Destructor
	 */
	~Variable(){}

	/**
	 * Override the comparison operators
	 */
	bool operator ==(const Variable& other) { return (m_ID == other.m_ID); }
	bool operator ==(const Variable* other) { return (m_ID == other->m_ID); }
	bool operator !=(const Variable& other) { return (m_ID != other.m_ID); }
	bool operator !=(const Variable* other) { return (m_ID != other->m_ID); }
	bool operator <(const Variable& other) { return (m_ID < other.m_ID);}
	bool operator <(const Variable* other) { return (m_ID < other->m_ID);}
	bool operator <=(const Variable& other) { return (m_ID <= other.m_ID);}
	bool operator <=(const Variable* other) { return (m_ID <= other->m_ID);}
	bool operator >(const Variable& other) { return (m_ID > other.m_ID);}
	bool operator >(const Variable* other) { return (m_ID < other->m_ID);}
	bool operator >=(const Variable& other) { return (m_ID >= other.m_ID);}
	bool operator >=(const Variable* other) { return (m_ID >= other->m_ID);}

	/**
	 * Gets the ID of this variable
	 */
	int getID() {
		return m_ID;
	}

	/**
	 * Gets the current domain size
	 */
	int getDomainSize() const {
		return m_domainSize;
	}

	/**
	 * Gets the initial domain size
	 */
	int getInitDomainSize() const {
		return m_initDomainSize;
	}

	/**
	 * Set the state of this discrete variable. Domain size is reduced to 1 as a result.
	 */
	void setValue(int value) {
		assert(value >= 0);
		assert(value < m_domainSize);

		int mappedValue = value;
		if (m_isReduced) {
			mappedValue = m_mapping[value];

		}
		m_mapping = vector<int>(1);
		m_mapping[0] = mappedValue;

		m_domainSize = 1;
		m_addressValue = 0;
	}

	/**
	 * Returns the state of this discrete variable.  Variable must be unary.
	 */
	int getValue() {
		assert(m_domainSize == 1);
		return m_mapping[0];
	}

	/**
	 * Set the address value of this discrete variable
	 */
	void setAddressValue(int value) {
		assert(value >= 0);
		assert(value < m_domainSize);
		m_addressValue = value;
	}

	/**
	 * Get the address value of this discrete variable
	 */
	int getAddressValue() {
		return m_addressValue;
	}

	/**
	 * Get min domain value (0-based)
	 */
	int getMinDomainValue() {
		if (m_isReduced) {
			return m_mapping[0];
		}
		else {
			return 0;
		}
	}

	/**
	 * Get max domain value (0-based)
	 */
	int getMaxDomainValue() {
		if (m_isReduced) {
			return m_mapping[m_domainSize - 1];
		}
		else {
			return m_domainSize - 1;
		}
	}

	/**
	 * Reduces the domain of this variable. Only the values corresponding to the
	 * <tt>true</tt> elements in the specified vector are retained. This vector
	 * must be the size of the initial domain of this variable.
	 */
/*	void reduce(vector<bool> valid) {
		assert(((int)valid.size()) == getInitDomainSize());

		// Create the mapping object
		m_mapping = vector<int>();
		for (int i = 0; i < getInitDomainSize(); i++) {
			if (valid[i]) {
				m_mapping.push_back(i);
			}
		}
		m_domainSize = m_mapping.size();
		assert(m_domainSize > 0);

		// If the result is a unary variable
		if (m_domainSize == 1) {
			m_addressValue = 0;
			m_value = 0;
		}
		else {
			m_addressValue = INVALID_VALUE;
			m_value = INVALID_VALUE;
		}

		if (m_domainSize < m_initDomainSize) {
			m_isReduced = true;
		}
		else {
			m_mapping.clear();
		}
	}*/

	/**
	 * Prints the id and address of this variable
	 */
	string toString() {
		std::stringstream ss;
		if (m_addressValue != INVALID_VALUE) {
			if (m_isReduced) {
				ss << m_ID << ":" << m_mapping[m_addressValue];
			}
			else {
				ss << m_ID << ":" << m_addressValue;
			}
		}
		else {
			ss << m_ID << ":";
		}
		return ss.str();
	}

	/* **************************************************************
	 *                     STATIC METHODS
	 * **************************************************************/

	/**
	 * Gets the index of the specified address assignment. For example, consider the
	 * ordered set of variables (A,B,C), which have domain sizes 2, 2 and 3,
	 * respectively. If the address assignment is (A=0,B=1,C=2) - i.e. then a
	 * value of: 0*6 + 1*3 + 2*1 = 5 will be returned.
	 */
	inline static int getAddress (const vector<Variable*>& variables) {
		int index = 0;
		int multiplier = 1;
		for (int i = (int)variables.size() - 1; i >= 0 ; i--) {
			assert(variables[i]->getAddressValue() != INVALID_VALUE);
			index += (multiplier * variables[i]->getAddressValue());
			multiplier *= variables[i]->getDomainSize();
		}
		return index;
	}

	/**
	 * Determines the assignment to each of the specified variables given the
	 * specified index.  For example, consider the ordered set of variables
	 * (A,B,C), which have domain sizes 2, 2 and 3, respectively. If given the
	 * index of 9, the variables in assignment will be modified so as to have
	 * the following values (A=1,B=1,C=1) as 1*6 + 1*3 + 0*1 = 9.
	 *
	 * Note: This method does not modify the m_value of each discrete variable,
	 * but uses m_addressValue instead.
	*/
	inline static void setAddress(const vector<Variable*>& variables, const int index) {
		int idx = index;
		assert(idx < getDomainSize(variables));
		for (int i = (int)variables.size() - 1; i >= 0 ; i--) {
			variables[i]->setAddressValue((idx % variables[i]->getDomainSize()));
			idx /= variables[i]->getDomainSize();
		}
	}

	/**
	 * Takes the Cartesian product of the specified variables
	 */
	inline static int getDomainSize(const vector<Variable*>& variables){
		int domain_size = 1;
		for(int i = 0; i < (int)variables.size(); i++) {
			domain_size*=variables[i]->getDomainSize();
		}
		return domain_size;
	}

	/**
	 * Takes the log of the Cartesian product of the specified variable domains
	 */
	inline static double getLogDomainSize(const vector<Variable*>& variables){
		double domain_size = 0;
		for(int i = 0; i < (int)variables.size(); i++) {
			domain_size += log(variables[i]->getDomainSize());
		}
		return domain_size;
	}

	/**
	 * Computes the size in MBs of the table needed to jointly represent the specified set of
	 * variables
	 */
	inline static double getSizeInMBs(const vector<Variable*>& variables) {
		double domain_size = 0;
		for(int i = 0; i < (int)variables.size(); i++) {
			domain_size += log(variables[i]->getDomainSize());
		}
		return exp(domain_size +  log(8) - log(1024) - log(1024));
	}


	/**
	 * Takes the Cartesian product of the specified variables
	 */
	inline static int getInitDomainSize(const vector<Variable*>& variables){
		int domain_size = 1;
		for(int i = 0; i < (int)variables.size(); i++) {
			domain_size*=variables[i]->getInitDomainSize();
		}
		return domain_size;
	}

};

#endif /* VARIABLE_H_ */
