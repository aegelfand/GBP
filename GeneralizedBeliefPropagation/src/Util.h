/*
 * Util.h
 *
 *  Created on: Apr 4, 2011
 *      Author: andrew.gelfand
 */

#ifndef UTIL_H_
#define UTIL_H_

#include <sstream>
#include <vector>
#include "Variable.h"

	const static double DEFAULT_TIME_LIMIT_SECS = 600;// Default time limit (in seconds)
	const double DEFAULT_MEMORY_LIMIT_MBS = 500; // Default memory limit (in MBs)

	static inline bool sortDesendInt (int i, int j) { return (i > j); }

	static inline bool compareByIDVar(Variable* A, Variable* B) {
		return (A->getID() < B->getID());
	}

	/**
	 * Compares two vectors of variables based upon the number of elements in each vector.
	 * Will put in descending order by size.
	 */
	static inline bool compareBySizeVar_Descend(const vector<Variable*>& A, const vector<Variable*>& B) {
		return (A.size() > B.size());
	}

	/**
	 * Compares two vectors of variables based upon the number of elements in each vector.
	 * Will put in asscending order by size.
	 */
	static inline bool compareBySizeVar_Ascend(const vector<Variable*>& A, const vector<Variable*>& B) {
		return (A.size() < B.size());
	}

	/**
	 * Gets the intersection of sets A and B and returns the result in C - i.e. C = A cap B
	 * Assumes that both <tt>A</tt> and <tt>B</tt> are sorted.
	 * @return bool -- Returns <tt>true</tt> if the intersection is not empty; <tt>false</tt> otherwise;
	 */
	inline bool getIntersectionVar(const vector<Variable*>& A, const vector<Variable*>& B, vector<Variable*>& C) {
		vector<Variable*> temp;
		temp.resize(A.size());
		vector<Variable*>::iterator curr,last;
		last = temp.end();
		curr = set_intersection(A.begin(), A.end(), B.begin(), B.end(), temp.begin(), compareByIDVar);
		temp.erase(curr,last);
		C = temp;
		return !C.empty();
	}

	/**
	 * Gets the intersection of sets A and B and returns the result in C - i.e. C = A cap B
	 * Assumes that both <tt>A</tt> and <tt>B</tt> are sorted.
	 * @return bool -- Returns <tt>true</tt> if the intersection is not empty; <tt>false</tt> otherwise;
	 */
	inline bool getIntersectionInt(const vector<int>& A, const vector<int>& B, vector<int>& C) {
		vector<int> temp;
		temp.resize(A.size());
		vector<int>::iterator curr,last;
		last = temp.end();
		curr = set_intersection(A.begin(), A.end(), B.begin(), B.end(), temp.begin());
		temp.erase(curr,last);
		C = temp;
		return !C.empty();
	}

	/**
	 * Gets the union of sets A and B and returns the result in C - i.e. C = A cup B
	 * Assumes that both <tt>A</tt> and <tt>B</tt> are sorted.
	 */
	inline void getUnionVar(const vector<Variable*>& A, const vector<Variable*>& B, vector<Variable*>& C) {
		vector<Variable*> temp;
		temp.resize(A.size() + B.size());
		vector<Variable*>::iterator curr,last;
		last = temp.end();
		curr = set_union(A.begin(), A.end(), B.begin(), B.end(), temp.begin(), compareByIDVar);
		temp.erase(curr,last);
		C = temp;
	}

	/**
	 * Gets the union of sets A and B and returns the result in C - i.e. C = A cup B
	 * Assumes that both <tt>A</tt> and <tt>B</tt> are sorted.
	 */
	inline void getUnionInt(const vector<int>& A, const vector<int>& B, vector<int>& C) {
		vector<int> temp;
		temp.resize(A.size() + B.size());
		vector<int>::iterator curr,last;
		last = temp.end();
		curr = set_union(A.begin(), A.end(), B.begin(), B.end(), temp.begin());
		temp.erase(curr,last);
		C = temp;
	}

	/**
	 * Gets the difference of sets A and B and returns the result in C - i.e. C = A - B
	 * Assumes that both <tt>A</tt> and <tt>B</tt> are sorted.
	 */
	inline void getDifferenceVar(const vector<Variable*>& A, const vector<Variable*>& B, vector<Variable*>& C) {
		vector<Variable*> temp;
		temp.resize(A.size());
		vector<Variable*>::iterator curr,last;
		last = temp.end();
		curr = set_difference(A.begin(), A.end(), B.begin(), B.end(), temp.begin(), compareByIDVar);
		temp.erase(curr,last);
		C = temp;
	}

	/**
	 * Gets the difference of sets A and B and returns the result in C - i.e. C = A - B
	 * Assumes that both <tt>A</tt> and <tt>B</tt> are sorted.
	 */
	inline void getDifferenceInt(const vector<int>& A, const vector<int>& B, vector<int>& C) {
		vector<int> temp;
		temp.resize(A.size());
		vector<int>::iterator curr,last;
		last = temp.end();
		curr = set_difference(A.begin(), A.end(), B.begin(), B.end(), temp.begin());
		temp.erase(curr,last);
		C = temp;
	}

	/**
	 * Removes the specified variable B from the set A and returns the result in C
	 */
	inline void removeVar(const vector<Variable*>& A, Variable* B, vector<Variable*>& C) {
		vector<Variable*> temp;
		temp.push_back(B);
		getDifferenceVar(A,temp,C);
	}

	/**
	 * Removes the specified variable B from the set A and returns the result in C
	 */
	inline void removeInt(const vector<int>& A, int B, vector<int>& C) {
		vector<int> temp;
		temp.push_back(B);
		getDifferenceInt(A,temp,C);
	}

	/**
	 * Adds the specified variable B to the set A and returns the result in C
	 */
	inline void addVar(const vector<Variable*>& A, Variable* B, vector<Variable*>& C) {
		vector<Variable*> temp;
		temp.push_back(B);
		getUnionVar(A,temp,C);
	}

	/**
	 * Adds the specified int B to the set A and returns the result in C
	 */
	inline void addInt(const vector<int>& A, int B, vector<int>& C) {
		vector<int> temp;
		temp.push_back(B);
		getUnionInt(A,temp,C);
	}

	inline string varsToString(const vector<Variable*>& A) {
		std::stringstream ss;
		ss << "(";
		for (int i = 0; i < (int)A.size(); i++) {
			if (i > 0) {
				ss << ",";
			}
			ss << A[i]->getID();
		}
		ss << ")";
		return ss.str();
	}

	inline string intsToString(const vector<int>& A) {
		std::stringstream ss;
		ss << "(";
		for (int i = 0; i < (int)A.size(); i++) {
			if (i > 0) {
				ss << ",";
			}
			ss << A[i];
		}
		ss << ")";
		return ss.str();
	}

	/**
	 * Returns <tt>true</tt> if A is a subseteq of B; returns <tt>false</tt> otherwise.
	 * Assumes that both <tt>A</tt> and <tt>B</tt> are sorted.
	 */
	inline bool isSubsetVar(const vector<Variable*>& A, const vector<Variable*>& B) {
		vector<Variable*> C;
		getIntersectionVar(A, B, C);
		return (A.size() == C.size());
	}

	/**
	 * Returns <tt>true</tt> if A is a subseteq of B; returns <tt>false</tt> otherwise.
	 * Assumes that both <tt>A</tt> and <tt>B</tt> are sorted.
	 */
	inline bool isSubsetInt(const vector<int>& A, const vector<int>& B) {
		vector<int> C;
		getIntersectionInt(A, B, C);
		return (A.size() == C.size());
	}

	/**
	 * Returns <tt>true</tt> if the set A contains variable B.
	 * Assumes that both <tt>A</tt> and <tt>B</tt> are sorted.
	 */
	inline bool containsVar(const vector<Variable*>& A, Variable* B) {
		return binary_search(A.begin(), A.end(), B, compareByIDVar);
	}

	/**
	 * Returns <tt>true</tt> if the set A and B are equivalent; returns <tt>false</tt> otherwise
	 */
	inline bool isEqualVar(const vector<Variable*>& A, const vector<Variable*>& B) {
		if (A.size() != B.size()) {
			return false;
		}
		for (int i = 0; i < (int)A.size(); i++) {
			if (A[i]->getID() != B[i]->getID()) {
				return false;
			}
		}
		return true;
	}

	/**
	 * Gets the symmetric difference of sets A and B and returns the result in C - i.e. C = {A - B} cup {B - A}
	 * Assumes that both <tt>A</tt> and <tt>B</tt> are sorted.
	 */
	inline void getSymmetricDifferenceVar(const vector<Variable*>& A, const vector<Variable*>& B, vector<Variable*>& C) {
		vector<Variable*> temp;
		temp.resize(A.size());
		vector<Variable*>::iterator curr,last;
		last = temp.end();
		curr = set_symmetric_difference(A.begin(), A.end(), B.begin(), B.end(), temp.begin(), compareByIDVar);
		temp.erase(curr,last);
		C = temp;
	}

#endif /* UTIL_H_ */
