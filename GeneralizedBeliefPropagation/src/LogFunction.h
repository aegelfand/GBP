/*
 * LogFunction.h
 *
 *  Created on: Mar 14, 2012
 *      Author: agelfand
 */

#ifndef LOGFUNCTION_H_
#define LOGFUNCTION_H_


#include <vector>
#include <cfloat>
#include <sstream>
#include <cassert>
#include <map>

#include "Util.h"
#include "Variable.h"

using namespace std;

struct LogFunction {

protected:
	/**
	 * The scope of this function, where the order of variables in
	 * this is placed in ascending order by ID.
	 */
	vector<Variable*> m_scopeAscending;

	/**
	 * The original scope of this function. Variables may not be in
	 * ascending order.
	 */
	vector<Variable*> m_scopeOriginal;

	/**
	 * Set of unary variables that have been removed from the scope of
	 * this function
	 */
	vector<Variable*> m_evidenceVariables;

	/**
	 * The values of this function. Values are kept in the log domain.
	 */
	vector<double> m_values;

	/**
	 * Number of values in this function
	 */
	int m_numValues;

	/**
	 * Unique ID for this function
	 */
	int m_ID;

	/**
	 * Count of functions (used to set the function's id)
	 */
	static int s_count;

	/**
	 * Representation of zero in the log domain
	 */
	const static double ZERO = -1000;

	/**
	 * Rearranges the elements of this discrete function according to the variable
	 * ordering in newScope.
	 * @param newScope - New scope
	 */
	void rearrange(vector<Variable*>& newScope) {
		assert(m_scopeAscending.size() == newScope.size());

		int card = Variable::getDomainSize(m_scopeAscending);
		vector<double> values = vector<double>(card);

		for (int i = 0; i < card; i++) {
			Variable::setAddress(m_scopeAscending, i);

			int idxInNewTable = Variable::getAddress(newScope);

			values[idxInNewTable] = m_values[i];
		}

		m_values.clear();
		m_scopeAscending.clear();
		m_values = values;
		m_scopeAscending = newScope;
	}

public:

	/**
	 * Compares functions a and b based upon the number of variables
	 * in each function's scope. Used to place in descending order by size.
	 */
	static inline bool compareBySize(LogFunction* a, LogFunction* b) {
		return (((int)a->getScope().size()) > ((int)b->getScope().size()));
	}

	/**
	 * Creates a trivial function with empty scope and a single
	 * value equal to 1.0.
	 */
	LogFunction() {
		m_values = vector<double>(1);
		m_values[0] = 0.0;
		m_numValues = 1;
		m_ID = s_count++;
	}

    /**
     * Constructs a new discrete valued function.
     *
     * NOTE: The scope vector (and table values) will be rearranged so that the
     * variables are in order from smallest ID to largest ID.
     * @param variables vector<Variable*> -- An ordered set of variables
     * comprising the scope of this discrete function.
     * @param values vector<double> -- A vector of values comprising this discrete
     * function. The values should be a 1-by-d vector, where d is the product of
     * the domain sizes of the variables in the scope of this function.  Tuples
     * are implicitly assumed in ascending order, with the last variable in the
     * scope as the 'least significant'. For example, consider a simple 3 node
     * network:
     *      A   B
     *       \ /
     *        C
     * where both A and C are binary nodes and B is ternary, so that the table
     * size for function C is d = 2*3*2 = 12. Further, assume that the variables
     * in the function are ordered as A, B, C (i.e. with C least significant). Then
     * the mapping to value indices is assumed to be as follows:
     *  | a | b | c || index
     *  | 0 | 0 | 0 ||  0
     *  | 0 | 0 | 1 ||  1
     *  | 0 | 1 | 0 ||  2
     *  | 0 | 1 | 1 ||  3
     *  | 0 | 2 | 0 ||  4
     *  | 0 | 2 | 1 ||  5
     *  | 1 | 0 | 0 ||  6
     *  | 1 | 0 | 1 ||  7
     *  | : | : | : ||  :
     *
     *  NOTE: The variables in the function's scope will be placed in ascending
     *  order by variable ID
     */
    LogFunction(const vector<Variable*> variables, const vector<double> values) {
    	// Validate the arguments
    	assert(Variable::getDomainSize(variables) == (int)values.size());

    	// Create a new function
    	m_scopeAscending = variables;
    	m_scopeOriginal = variables;
    	m_numValues = (int)values.size();

    	m_values = vector<double>(m_numValues);
		for (int i = 0; i < m_numValues; i++) {
			m_values[i] = log(values[i]);
			if (m_values[i] < ZERO) {
				m_values[i] = ZERO;
			}
		}

		m_ID = s_count++;

		// Sort the function scope from smallest ID to largest ID and
		// rearrange the table values accordingly
		vector<Variable*> scope = vector<Variable*>(variables);
		std::sort(scope.begin(), scope.end(), compareByIDVar);
		rearrange(scope);

		// Remove any unary variables from the scope
		vector<Variable*>::iterator it;
		int numEvidence = 0;
		for (int i = (int)m_scopeAscending.size()-1; i >= 0; i--) {
			if (m_scopeAscending[i]->getDomainSize() == 1) {
				m_evidenceVariables.push_back(m_scopeAscending[i]);
				it = m_scopeAscending.begin();
				m_scopeAscending.erase(it+i-numEvidence);
				numEvidence++;
				i--;
			}
		}
    }

    /**
     * Copy constructor
     */
    LogFunction(LogFunction& func) {
    	m_scopeAscending = func.m_scopeAscending;
    	m_scopeOriginal = func.m_scopeOriginal;
    	m_values = vector<double>(func.m_values);
    	m_numValues = func.m_numValues;
    	m_ID = s_count++;
    }

	/**
	 * Creates a new uniform function over the specified scope.
	 */
	LogFunction(const vector<Variable*>& variables) {
		assert(!variables.empty());

		// Create a new function
		m_scopeAscending = variables;
		m_scopeOriginal = variables;
		m_numValues = Variable::getDomainSize(variables);
		m_values = vector<double>(m_numValues);
		m_ID = s_count++;

		// Sort the function scope from smallest ID to largest ID
		std::sort(m_scopeAscending.begin(), m_scopeAscending.end(),compareByIDVar);

		// Initialize uniformly
		std::fill(m_values.begin(),m_values.end(),0.0);
	}

	/**
	 * Creates a new uniform function over the specified variable.
	 */
	LogFunction(Variable* variable) {
		// Create a new function
		m_scopeAscending.push_back(variable);
		m_scopeOriginal.push_back(variable);
		m_numValues = Variable::getDomainSize(m_scopeAscending);
		m_values = vector<double>(m_numValues);
		m_ID = s_count++;

		// Initialize uniformly
		std::fill(m_values.begin(),m_values.end(),0.0);
	}

	/**
	 * Clears this function by making all of its values 1
	 */
	void clear() {
		std::fill(m_values.begin(),m_values.end(),0.0);
	}

    void assign(LogFunction& func) {
    	m_scopeAscending = func.m_scopeAscending;
    	m_scopeOriginal = func.m_scopeOriginal;
		m_values = func.m_values;
		m_numValues = func.m_numValues;
    }

	/**
	 * Gets this function's ID
	 */
	int getID() {
		return m_ID;
	}

    /**
     * Gets the working scope of this function. Note that
     * the variables are placed in ascending order by ID.
     */
    const vector<Variable*>& getScope() {
        return m_scopeAscending;
    }

    /**
     * Gets the original scope of this function. Variables
     * may not be in ascending order by ID.
     */
    const vector<Variable*>& getOriginalScope() {
        return m_scopeOriginal;
    }

    /**
     * Gets the log of the value from this function at the specified index.
     */
    double getLogValueAt(int index) {
    	assert(index >= 0);
    	assert(index < m_numValues);
        return m_values[index];
    }

    /**
	 * Gets the value from this function at the specified index.
	 */
	double getValueAt(int index) {
		assert(index >= 0);
		assert(index < m_numValues);
		return exp(m_values[index]);
	}

	/**
	 * Gets the largest value in this table
	 */
	double getMaxValue() {
		return exp(*max_element(m_values.begin(), m_values.end()));
	}

	/**
	 * Gets the largest log value in this table
	 */
	double getMaxLogValue() {
		return *max_element(m_values.begin(), m_values.end());
	}

	/**
     * Gets the number of values in this table
     */
    int getNumValues() {
    	return m_numValues;
    }

    /**
     * Gets a string representation of this function, for printing
     */
    string toString() {
    	std::stringstream ss;
    	for (unsigned int i = 0; i < m_values.size(); i++) {
			Variable::setAddress(m_scopeAscending, i);
			// Output index
			ss << i << " : ";

			if (m_scopeAscending.empty()) {
				ss << "[]" << m_values[i] << endl;
			}
			else {
				ss << "[";
				for (int j = 0; j < (int)m_scopeAscending.size(); j++) {
					if (j != 0) {
						ss << " ";
					}
					ss << m_scopeAscending[j]->toString();
				}
				ss << "] : " << exp(m_values[i]) << endl;
			}
		}

    	return ss.str();
    }

    /**
     * Reduces the number of values in this discrete function table by removing
     * all table entries inconsistent with the evidence applied to each variable.
     * Assumes that each
     */
    void applyEvidence() {
    	// If the table needs to be updated...
    	int numValues = Variable::getDomainSize(m_scopeAscending);
    	if (m_numValues > numValues) {

    		// Table reduction procedures assumes that table was constructed with
    		// initial domain sizes
    		assert(m_numValues == Variable::getInitDomainSize(m_scopeAscending));

    		// New, reduced table of values
    		vector<double> values = vector<double>(numValues);

    		for (int i = 0; i < numValues; i++) {

    			Variable::setAddress(m_scopeAscending, i);

    			int index = 0;
				int multiplier = 1;
				for (int j = (int)m_scopeAscending.size() - 1; j >= 0 ; j--) {
					if (m_scopeAscending[j]->getDomainSize() == 1) {
						index += (multiplier * m_scopeAscending[j]->getValue());
					}
					else {
						assert(m_scopeAscending[j]->getAddressValue() != INVALID_VALUE);
						index += (multiplier * m_scopeAscending[j]->getAddressValue());
					}

					multiplier *= m_scopeAscending[j]->getInitDomainSize();
				}

				values[i] = m_values[index];
    		}

    		// Update the table
    		m_values.clear();
    		m_values = values;
    		m_numValues = values.size();

    		// Remove any unary variables from the scope
			vector<Variable*>::iterator it;
			int numEvidence = 0;
			for (int i = (int)m_scopeAscending.size()-1; i >= 0; i--) {
				if (m_scopeAscending[i]->getDomainSize() == 1) {
					m_evidenceVariables.push_back(m_scopeAscending[i]);
					it = m_scopeAscending.begin();
					m_scopeAscending.erase(it+i-numEvidence);
					numEvidence++;
					i--;
				}
			}
    	}
    }

    /**
	 * Multiplies this discrete function by the specified set of discrete
	 * functions, while simultaneously marginalizing out the specified set of variables.
	 * @param phis vector<Function*> -
	 * @param variables vector<Variable*> -
	 */
    void multiplyAndMarginalizeOut(const vector<LogFunction*>& phis, const vector<Variable*>& variables) {
    	// Determine the union of the scopes of this factor and the other factors
    	vector<Variable*> unionOfScopes = vector<Variable*>(m_scopeAscending);
    	for (unsigned int i = 0; i < phis.size(); i++) {
    		getUnionVar(unionOfScopes, phis[i]->getScope(), unionOfScopes);
    	}

    	// Cannot have an empty scope and marginalize variables
    	assert(!(unionOfScopes.empty() && !variables.empty()));
    	if (unionOfScopes.empty() && variables.empty()) {
    		m_scopeAscending = vector<Variable*>();
			m_values = vector<double>(1);
			m_values[0] = 0.0;
			m_numValues = 1;
			m_ID = s_count++;
    		return;
    	}

    	// Determine the set of variables remaining after the multiplication
    	// and summation operation
    	vector<Variable*> newScope;
    	getDifferenceVar(unionOfScopes, variables, newScope);

    	// Determine the cardinality of the new function
    	int cardOfNewFunction = Variable::getDomainSize(newScope);
    	vector<double> newValues = vector<double>(cardOfNewFunction);
    	vector<double> maxValues = vector<double>(cardOfNewFunction);
    	fill(maxValues.begin(), maxValues.end(), ZERO);

    	// Determine the cardinality of the intermediate factor due to combination
    	int cardOfUnion = Variable::getDomainSize(unionOfScopes);

    	// Determine the maximum value in the summations at each element
    	for (int i = 0; i < cardOfUnion; i++) {
    		// Set the assignment in union to the new assignment
    		Variable::setAddress(unionOfScopes, i);

    		// Determine the requisite value in this function...
			int idx = Variable::getAddress(m_scopeAscending);
			double temp = m_values[idx];

			// ...and all of the other functions to be combined
			for (unsigned int j = 0; j < phis.size(); j++) {
				idx = Variable::getAddress(phis[j]->getScope());
				assert((idx >= 0) && (idx < phis[j]->m_numValues));
				temp += phis[j]->m_values[idx];
			}

			// Determine the index in the new function
			int idx3 = Variable::getAddress(newScope);

			if (temp > maxValues[idx3]) {
				maxValues[idx3] = temp;
			}
    	}

    	// For each assignment p(i)
    	for (int i = 0; i < cardOfUnion; i++) {
    		// Set the assignment in union to the new assignment
    		Variable::setAddress(unionOfScopes, i);

    		// Determine the requisite value in this function...
    		int idx = Variable::getAddress(m_scopeAscending);
    		double temp = m_values[idx];

    		// ...and all of the other functions to be combined
    		for (unsigned int j = 0; j < phis.size(); j++) {
    			idx = Variable::getAddress(phis[j]->getScope());
    			assert(idx >= 0);
    			assert(idx < phis[j]->getNumValues());
    			temp += phis[j]->m_values[idx];
    		}

    		// Determine the index in the new function
    		int idx3 = Variable::getAddress(newScope);

    		newValues[idx3] += exp(temp - maxValues[idx3]);
    		assert(!isnan(newValues[idx3]));
    	}

    	for (int i = 0; i < cardOfNewFunction; i++) {
    		newValues[i] = maxValues[i] + log(newValues[i]);
    		if (newValues[i] < ZERO) {
    			newValues[i] = ZERO;
			}
    	}

    	m_scopeAscending.clear();
    	m_values.clear();
    	m_scopeAscending = newScope;
    	m_values = newValues;
    	m_numValues = m_values.size();
    }

    inline void multiplyAndMarginalizeOut(const vector<LogFunction*>& phis, Variable* var) {
    	vector<Variable*> vars;
    	vars.push_back(var);
    	multiplyAndMarginalizeOut(phis,vars);
    }

    /**
	 * Exponentiates the tabular function. Raises each element to the specified power.
	 */
	inline void exponentiateBy(double power) {
		for (int i = 0; i < m_numValues; i++) {
			m_values[i] *= power;
		}
	}

	/**
	 * Exponentiates the tabular function. Raises each element to the specified power.
	 * Note:  This function is not modified.
	 */
	inline LogFunction* exponentiate(double power) {
		if (power == 1.0) {
			return new LogFunction(*this);
		}
		// Create the new function
		LogFunction* newFunc = new LogFunction();
		newFunc->m_scopeAscending = m_scopeAscending;
		newFunc->m_values = m_values;
		newFunc->m_evidenceVariables = m_evidenceVariables;
		newFunc->m_numValues = m_numValues;

		for (int i = 0; i < m_numValues; i++) {
			newFunc->m_values[i] *= power;
		}
		return newFunc;
	}

	/**
	 * Multiples by the specified function in place
	 */
	inline void multiplyBy(double c) {
		for(int i = 0; i < m_numValues; i++) {
			m_values[i] += log(c);
		}
	}

    /**
     * Multiples by the specified function in place
     */
    inline void multiplyBy(LogFunction* phi) {
    	// Determine the union of the scopes of this factor and the other factors
		vector<Variable*> unionOfScopes = vector<Variable*>(m_scopeAscending);
		getUnionVar(unionOfScopes, phi->getScope(), unionOfScopes);

		assert(!unionOfScopes.empty());

		// Determine the cardinality of the intermediate factor due to combination
		int cardOfUnion = Variable::getDomainSize(unionOfScopes);
		vector<double> newValues = vector<double>(cardOfUnion);

		// For each assignment in the union
		for (int i = 0; i < cardOfUnion; i++) {
			// Set the assignment in union to the new assignment
			Variable::setAddress(unionOfScopes, i);

			// Determine the requisite value in this function...
			int idx = Variable::getAddress(m_scopeAscending);
			assert(idx >= 0);
			assert(idx < m_numValues);
			newValues[i] = m_values[idx];

			// ...and in the other function
			idx = Variable::getAddress(phi->getScope());
			assert(idx >= 0);
			assert(idx < phi->getNumValues());
			newValues[i] += phi->m_values[idx];

			assert(!isnan(newValues[i]));
			assert(!isinf(newValues[i]));
		}

		m_scopeAscending.clear();
		m_values.clear();
		m_scopeAscending = unionOfScopes;
		m_values = newValues;
		m_numValues = cardOfUnion;
    }

    /**
     * Marginalizes the specified set of variables from this function.
     * Assumes that the variables in <tt>variables</tt> are placed in
     * ascending order by ID.
     * @param variables --
     */
    inline void marginalizeOut(vector<Variable*> variables) {
    	// Nothing to be marginalized
		if (variables.empty()) {
			return;
		}

		// Determine the set of variables remaining after summation operation
		vector<Variable*> newScope;
		getDifferenceVar(m_scopeAscending, variables, newScope);

		// Determine the cardinality of the new function
		int cardOfNewFunction = Variable::getDomainSize(newScope);
		vector<double> newValues = vector<double>(cardOfNewFunction);
		vector<double> maxValues = vector<double>(cardOfNewFunction);
		fill(maxValues.begin(), maxValues.end(), ZERO);

		// Determine the maximum value in the summations at each element
		for (int i = 0; i < m_numValues; i++) {
			Variable::setAddress(m_scopeAscending,i);

			// Determine the index in the new function
			int idx = Variable::getAddress(newScope);

			if (m_values[i] > maxValues[idx]) {
				maxValues[idx] = m_values[i];
			}
		}

		// For each assignment in the union
		for (int i = 0; i < m_numValues; i++) {
			Variable::setAddress(m_scopeAscending,i);

			// Determine the index in the new function
			int idx = Variable::getAddress(newScope);

			newValues[idx] += exp(m_values[i] - maxValues[idx]);
			assert(!isnan(newValues[idx]));
		}

		for (int i = 0; i < cardOfNewFunction; i++) {
			newValues[i] = maxValues[i] + log(newValues[i]);
			if (newValues[i] < ZERO) {
				newValues[i] = ZERO;
			}
		}

		// Create the new function
		m_scopeAscending = newScope;
		m_values = newValues;
		m_numValues = cardOfNewFunction;
    }

    /**
      * Marginalizes the specified set of variables from this function. This
      * function is not modified.  Assumes that the variables in <tt>variables</tt>
      * are placed in ascending order by ID.
      * @param variables --
      */
     inline LogFunction* marginalize(vector<Variable*> variables) {
     	// Nothing to be marginalized
     	if (variables.empty()) {
     		LogFunction* newFunc = new LogFunction();
 			newFunc->m_scopeAscending = m_scopeAscending;
 			newFunc->m_values = m_values;
 			newFunc->m_evidenceVariables = m_evidenceVariables;
 			newFunc->m_numValues = m_numValues;

 			return newFunc;
     	}

     	// Determine the set of variables remaining after summation operation
 		vector<Variable*> newScope;
 		getDifferenceVar(m_scopeAscending, variables, newScope);

     	// Determine the cardinality of the new function
     	int cardOfNewFunction = Variable::getDomainSize(newScope);
     	vector<double> newValues = vector<double>(cardOfNewFunction);
     	vector<double> maxValues = vector<double>(cardOfNewFunction);
     	fill(maxValues.begin(), maxValues.end(), ZERO);

		// Determine the maximum value in the summations at each element
		for (int i = 0; i < m_numValues; i++) {
			Variable::setAddress(m_scopeAscending,i);

			// Determine the index in the new function
			int idx = Variable::getAddress(newScope);

			if (m_values[i] > maxValues[idx]) {
				maxValues[idx] = m_values[i];
			}
		}

		// For each assignment in the union
		for (int i = 0; i < m_numValues; i++) {
			Variable::setAddress(m_scopeAscending,i);

			// Determine the index in the new function
			int idx = Variable::getAddress(newScope);

			newValues[idx] += exp(m_values[i] - maxValues[idx]);
			assert(!isnan(newValues[idx]));
		}

		for (int i = 0; i < cardOfNewFunction; i++) {
			newValues[i] = maxValues[i] + log(newValues[i]);
			if (newValues[i] < ZERO) {
				newValues[i] = ZERO;
			}
		}

     	// Create the new function
     	LogFunction* newFunc = new LogFunction();
     	newFunc->m_scopeAscending = newScope;
     	newFunc->m_values = newValues;
     	newFunc->m_evidenceVariables = m_evidenceVariables;
     	newFunc->m_numValues = cardOfNewFunction;

     	return newFunc;
     }

     /**
      * Computes the entropy of this tabular function.
      *   H(X) = SUM{p(x)log(p(x))}
      * where the sum is over all configurations of the variables X
      * in this function.
      */
     double getEntropy() {
    	 double entropy = 0;
    	 for (int i = 0; i < m_numValues; i++) {
    		 entropy -= exp(m_values[i])*m_values[i];
    	 }
    	 return entropy;
     }

     /**
      * Normalizes this function so that the sum of its elements equals 1.
      */
     void normalize() {
    	 // Get the maximum value in this function
    	 double maxValue = *max_element(m_values.begin(), m_values.end());
    	 if (maxValue < ZERO) {
    		 cout << "WARNING: maxValue < 1e-500!!!" << endl;
    	 }

    	 // Compute the sum: exp(log( (a + b + c + ....) / max))
    	 double total = 0.0;
    	 for (int i = 0; i < m_numValues; i++) {
    		 total += exp(m_values[i] - maxValue);
    	 }

    	 // Compute: log(a + b + c + ...)
    	 total = log(total) + maxValue;
    	 assert(total >= ZERO);

    	 // Now normalize in the log domain
    	 for (int i = 0; i < m_numValues; i++) {
    		 m_values[i] = m_values[i] - total;
    	 }
     }

     /**
      * Gets the log normalization constant for this function
      */
     double getlogPR() {
    	 // Get the maximum value in this function
		 double maxValue = *max_element(m_values.begin(), m_values.end());

		 // Compute the sum: exp(log( (a + b + c + ....) / max))
		 double total = 0.0;
		 for (int i = 0; i < m_numValues; i++) {
			 total += exp(m_values[i] - maxValue);
		 }

		 // Compute: log(a + b + c + ...)
		 total = log(total) + maxValue;

		 return total;
     }

     /**
	  * Divides this factor by the specified constant.
	  */
	 inline void divideBy(double c) {
		 for (int i = 0; i < m_numValues; i++) {
			 m_values[i] -= log(c);
		 }
	 }

	/**
	 * Divides this factor by the specified factor.  The scope of the denominator
	 * must be a subseteq of the variables in this factor's scope.
	 * @param phi -
	 */
	inline void divideBy(LogFunction* phi) {
		// Validate the divisor
		assert(isSubsetVar(phi->getScope(), m_scopeAscending));

		// Determine the cardinality of the new function
		vector<double> newValues = vector<double>(m_numValues);

		// For each assignment
		for (int i = 0; i < m_numValues; i++) {
			Variable::setAddress(m_scopeAscending,i);

			// Determine the index in the new function
			int idx = Variable::getAddress(phi->getScope());

			newValues[i] = m_values[i] - phi->m_values[idx];
		}

		// Update this function
		m_values = newValues;
	}

    /**
     * Adds this function to the specified function, using the specified
     * weight.  In particular, if this function is phi1 and the specified
     * function is phi2, each element is added as:
     * 		phi1[i] = (1.0 - alpha)*phi1[i] + alpha*phi2[i]
     * where lambda is a weighting factor assumed to be in [0,1].
     * @param phi2 -
     * @param alpha
     */
    void weightedAdd(LogFunction* phi2, double alpha) {
    	// Validate the args the
    	assert((alpha >= 0) && (alpha <= 1.0));
    	assert(m_scopeAscending.size() == phi2->getScope().size());
    	vector<Variable*> isect;
    	getIntersectionVar(m_scopeAscending, phi2->getScope(), isect);
    	assert(isect.size() == m_scopeAscending.size());

    	// Catch the boundary conditions
    	if (alpha == 0) {
    		return;
    	}
    	else if (alpha == 1) {
    		m_values = phi2->m_values;
    		return;
    	}

    	vector<double> newValues = vector<double>(m_numValues);
    	for (int i = 0; i < m_numValues; i++) {

    		if (m_values[i] > phi2->m_values[i]) {
    			newValues[i] = log((1.0 - alpha)*exp(m_values[i] - m_values[i]) + alpha*exp(phi2->m_values[i] - m_values[i])) + m_values[i];
    		}
    		else {
    			newValues[i] = log((1.0 - alpha)*exp(m_values[i] - phi2->m_values[i]) + alpha*exp(phi2->m_values[i] - phi2->m_values[i])) + phi2->m_values[i];
    		}
    	}

    	// Update the function
    	m_values = newValues;
    }

    /**
     * Multiplies this function by the specified function, using the specified
     * weight to exponentiate.  In particular, if this function is phi1 and the
     * specified function is phi2, each element is computed as:
     * 		phi1[i] = phi1[i]^(1.0 - alpha) * phi2[i]^alpha
     * where lambda is a weighting factor assumed to be in [0,1].
     * @param phi2 -
     * @param alpha
     */
    void weightedMultiply(LogFunction* phi2, double alpha) {
    	// Validate the args
    	assert((alpha >= 0) && (alpha <= 1.0));
    	assert(m_scopeAscending.size() == phi2->getScope().size());
    	vector<Variable*> isect;
    	getIntersectionVar(m_scopeAscending, phi2->getScope(), isect);
    	assert(isect.size() == m_scopeAscending.size());

    	// Catch the boundary conditions
    	if (alpha == 0) {
    		return;
    	}
    	else if (alpha == 1) {
    		m_values = phi2->m_values;
    		return;
    	}

    	vector<double> newValues = vector<double>(m_numValues);
    	for (int i = 0; i < m_numValues; i++) {
    		newValues[i] = (1.0 - alpha)*m_values[i] + alpha*phi2->m_values[i];
    	}

    	// Update the function
    	m_values = newValues;
    }

    /**
     * Computes the maximum of the absolute difference between the values in
     * the specified discrete functions.  The two functions must have the same
     * number of elements and are assumed to be over the same scope.
     */
    static double getDistanceLinf(LogFunction* phi1, LogFunction* phi2) {
    	assert(phi1->getNumValues() == phi2->getNumValues());
    	double maxDiff = 0.0;

    	for (int i = 0; i < phi1->getNumValues(); i++) {
			double diff = fabs(exp(phi1->m_values[i]) - exp(phi2->m_values[i]));
			if (diff > maxDiff) {
				maxDiff = diff;
			}
		}

    	return maxDiff;
    }

    /**
	* Computes the L1 distance between tabular function phi1 and tabular
	* function phi2.
	*/
    static double getDistanceL1(LogFunction* phi1, LogFunction* phi2) {
    	assert(phi1->getNumValues() == phi2->getNumValues());
    	double distance = 0;

		for (int i = 0; i < phi1->getNumValues(); i++) {
			distance += fabs(exp(phi1->m_values[i]) - exp(phi2->m_values[i]));
		}

		return distance;
    }

};

#endif /* LOGFUNCTION_H_ */
