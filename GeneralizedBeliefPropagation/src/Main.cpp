/*
 * Main.cpp
 *
 *  Created on: Dec 13, 2011
 *      Author: agelfand
 */

#include <bitset>
#include <fstream>
#include <iostream>
#include <map>
#include <math.h>
#include <queue>
#include <set>
#include <sstream>
#include <string>
#include <time.h>
#include <utility>
#include <vector>

#include "binaryheap.h"
#include "LogFunction.h"
#include "Graph.h"
#include "JoinGraph.h"
#include "myRandom.h"
#include "RegionGraph.h"
#include "Partition.h"
#include "Variable.h"

// Region Graph Design Type
typedef enum {
	BETHE, IJGP
} REGION_DESIGN_TYPE;

//
// Graphical Model
//
vector<LogFunction*> m_functions;
vector<Variable*> m_variables;
int m_numVariables;
int m_numFunctions;
Graph* m_graph; // The undirected graph representing the markov network

//
// The Region graph
//
vector<Region*> m_regions;
vector<Edge*> m_edges;

// Any contribution to the log-partition function coming from isolated regions
double m_logPR_Partial;

//
// Run-time parameters
//
int m_maxIters = 1000;
double m_convergenceThreshold = 1e-8;
double m_alpha = 0.5; // Default step size
bool m_isDoubleLoop = false;
int m_innerLoopMaxIters = 10;
bool m_isVerbose = false;

myRandom* m_random;
int m_seed;

string m_uaiFileName;
string m_rgFileName;
REGION_DESIGN_TYPE m_rgType = BETHE;

bool m_isFromFile = false;
bool m_isCVM = false;
int m_iBound = 1;

/**
 * Displays usage
 */
void help() {
	cout << "Usage:" << endl;
	cout << "  GeneralizedBeliefPropagation <uaifilename> <seed>" << endl << endl;
	cout << "  Runs Genearlized Belief Propagation (GBP) on the network specified" << endl;
	cout << "  by <uaifilename>. The network must be in .uai format." << endl;
	cout << "  Computes: " << endl;
	cout << "      PR - The log10 value of the partition function (probability of evidence)" << endl;
	cout << "           and outputs to a file named uaifilename.PR" << endl;
	cout << "      BEL - Computes the belief approximations for all cliques in the model" << endl;
	cout << "           and outputs to a file named uaifilename.BEL" << endl;
	cout << " Optional Arguments:" << endl;
	cout << "    -iters <value>  maximum number of iterations to run message passing for." << endl;
	cout << "    -alpha <value>  step size used in message passing." << endl;
	cout << "    -thresh <value>  convergence threshold." << endl;
	cout << "    -dbl  runs a convergent, double-loop form of message passing." << endl;
	cout << "    -iters_inner <value>  maximum number of iterations to run inner loop" << endl;
	cout << "          of double-loop algorithm for." << endl;
	cout << "    -cvm  runs GBP on a region graph (RG) constructed using the " << endl;
	cout << "		   cluster variation method (cvm)." << endl;
	cout << "    -ijgp <value>  Runs the mini-bucket schematic, where <value>" << endl;
	cout << "          specifies the iBound." << endl;
	cout << "    -rg <rgfilename>" << endl;
	cout << "      Constructs the a region graph using the outer regions specified in" << endl;
	cout << "      the <rgfilename>. Inner regions are filled in using cvm." << endl;
	cout << "    -verbose  Runs in verbose mode." << endl;
}

void readParameters(int argc, char* argv[]) {
	// Must include model file and seed
	if (argc < 3) {
		help();
		exit(-1);
	}

	int idx = 1;
	string next;
	while (idx < (argc - 2)) {
		next = argv[idx];
		if (next == "-cvm") {
			m_isCVM = true;
		}
		else if (next == "-iters") {
			idx++;
			m_maxIters = atoi(argv[idx]);
			assert(m_maxIters > 0);
		}
		else if (next == "-iters_inner") {
			idx++;
			m_innerLoopMaxIters = atoi(argv[idx]);
			assert(m_innerLoopMaxIters > 0);
		}
		else if (next == "-alpha") {
			idx++;
			m_alpha = atof(argv[idx]);
			assert(m_alpha >= 0);
			assert(m_alpha <= 1.0);
		}
		else if (next == "-thresh") {
			idx++;
			m_convergenceThreshold = atof(argv[idx]);
			assert(m_convergenceThreshold >= 0);
			assert(m_convergenceThreshold < 1.0);
		}
		else if (next == "-dbl") {
			m_isDoubleLoop = true;
		}
		else if (next == "-ijgp") {
			idx++;
			m_rgType = IJGP;
			m_iBound = atoi(argv[idx]);
			assert(m_iBound > 0);
			assert(!m_isFromFile);
		}
		else if (next == "-rg") {
			idx++;
			m_rgFileName = argv[idx];
			m_isFromFile = true;
		}
		else if (next == "-verbose") {
			m_isVerbose = true;
		}
		else {
			// Unexpected option
			help();
			exit(-1);
		}
		idx++;
	}

	// Get the model file and seed
	m_uaiFileName = argv[idx]; idx++;
	m_random = new myRandom(atoi(argv[idx])); idx++;
}

/**
 * Reads a network file from disk.
 * @param uaiFileName - Name of the file containing the network
 */
void readUAI(const char* uaiFileName) {
	ifstream infile(uaiFileName);
	cout << "Reading model from file: " << uaiFileName << endl;

	if (!infile.is_open()) {
		cout << "Failed to open file:" << uaiFileName << endl;
		exit(-1);
	}

	// Validate the network type
	string tempString;
	infile >> tempString;
	if ((tempString.compare("BAYES") != 0) && (tempString.compare("MARKOV") != 0))  {
		cerr << "Unexpected Network Type: " << tempString.c_str() << endl;
		exit(-1);
	}

	// Determine the number of variables
	int numVariables;
	infile >> numVariables;
	m_numVariables = numVariables;

	// Read the domain sizes
	m_variables = vector<Variable*> (numVariables);
	for (int i = 0; i < numVariables; i++) {
		int domainSize;
		infile >> domainSize;

		// And create a corresponding variable
		m_variables[i] = new Variable(i, domainSize);
	}

	// Get the number of functions
	infile >> m_numFunctions;

	// Get the scope of each function
	vector<vector<Variable*> > scope(m_numFunctions);
	for (int i = 0; i < m_numFunctions; i++) {
		int numVarsInScope;
		infile >> numVarsInScope;

		scope[i] = vector<Variable*> (numVarsInScope);
		for (int j = 0; j < numVarsInScope; j++) {
			int varID;
			infile >> varID;
			scope[i][j] = m_variables[varID];
		}
	}

	cout << "Reading Network File...";

	// Create each discrete-valued function
	m_functions = vector<LogFunction*> (m_numFunctions);
	for (int i = 0; i < m_numFunctions; i++) {
		// Get number of table entries
		int numEntries;
		infile >> numEntries;

		// Ensure that the number of values agrees with the product
		// of the domains from the preamble
		int card = Variable::getDomainSize(scope[i]);
		assert(numEntries == card);

		// Create an appropriately sized table and
		// populate with the probabilities/values
		vector<double> values = vector<double>(card);
		for (int j = 0; j < card; j++) {
			double value;
			infile >> value;

			values[j] = value;
		}

		// Create a new function and add it to the model
		m_functions[i] = new LogFunction(scope[i], values);
	}

	cout << "done!" << endl;
	infile.close();
}

/**
 * Computes an estimate of the log partition function value given the current
 * region graph.
 */
double getLogPR() {
	double avgEnergy = 0.0;
	double entropy = 0.0;

	for (int i = 0; i < (int)m_regions.size(); i++) {
		avgEnergy += m_regions[i]->getCountingNumber()*m_regions[i]->getEnergy();
		entropy += m_regions[i]->getCountingNumber()*m_regions[i]->getBelief()->getEntropy();
	}

	return (avgEnergy + entropy + m_logPR_Partial);
}

void buildGraph() {
	// Create the graph structure
	m_graph = new Graph(m_numVariables);
	for (int i = 0; i < (int)m_functions.size(); i++) {
		for (int j = 0; j < (int)m_functions[i]->getScope().size(); j++) {
			for (int k = j + 1; k < (int)m_functions[i]->getScope().size(); k++) {
				m_graph->addEdge(m_functions[i]->getScope()[j]->getID(), m_functions[i]->getScope()[k]->getID());
			}
		}
	}
}

/**
 * Get an elimination order
 */
vector<int> getEliminationOrder(vector<vector<int> >& cliques) {
	vector<int> order;
	cliques.clear();
	int minWidth = m_graph->getNumberOfVertices();
	int numVariables = m_graph->getNumberOfVertices();

	// Get an elimination ordering using the min-fill heuristic
	for (int i = 0; i < 100; i++) {
		Graph* g = new Graph(*m_graph);

		vector<int> order_i;
		vector<vector<int> > cliques_i;

		int width = g->getMinFillOrdering(order_i,cliques_i,m_random,minWidth);

		// Keep track of the min-width ordering
		if ((width > 0) && (width < minWidth)) {
			order = order_i;
			cliques = cliques_i;
			assert(numVariables == (int)order.size());
			assert(numVariables == (int)cliques.size());
			minWidth = width;
		}

		if (width == 1) {
			break;
		}

		delete g;
	}
	cout << "Max cluster size of order: " << minWidth + 1 << endl;

	return order;
}

vector<int> getEliminationOrder() {
	vector<int> order;
	vector<vector<int> > cliques;
	return getEliminationOrder(cliques);
}

/**
 * Gets a sequence of maximal cliques C_1,...,C_n consistent with
 * some elimination order.
 */
vector<vector<Variable*> > getMaximalCliqueSequence() {
	vector<vector<int> > cliques;
	vector<int> order = getEliminationOrder(cliques);

	// Sort the vertices in each clique
	for (int i = 0; i < (int)cliques.size(); i++) {
		sort(cliques[i].begin(), cliques[i].end());
	}

	// Remove non-maximal clusters from the cluster sequence. If C_1,...,C_n
	// is the cluster sequence and C_j is a subset of C_i, where i < j and i
	// is the largest such index in the sequence, then C_j is a non-maximal
	// cluster and can be replaced by cluster C_i in the sequence
	bool mergedClusters = true;
	int idx_j = cliques.size() - 1;
	while(mergedClusters) {
		mergedClusters = false;
		vector<vector<int> >::iterator it;

		// For j = n...1, find cluster i, such that C_j subseteq C_i
		for (int j = idx_j; j >= 0; j--) {
			for (int i = j - 1; i >= 0; i--) {
				// If C_j subseseq C_i
				if (isSubsetInt(cliques[j], cliques[i])) {

					// Replace C_j with C_i
					it = cliques.begin();
					vector<int> clique_i = cliques[i];
					cliques.erase(it + i);
					it = cliques.begin();
					cliques.erase(it + j - 1);
					it = cliques.begin();
					cliques.insert(it + j - 1, clique_i);

					mergedClusters = true;
					idx_j = j - 1;
					break;
				}
			}
		}
	}

	// Create a variable cluster for each maximal vertex clique
	vector<vector<Variable*> > clusterSequence;
	for (int i = 0; i < (int)cliques.size(); i++) {
		vector<Variable*> clusterVars;
		for (int j = 0; j < (int)cliques[i].size(); j++) {
			clusterVars.push_back(m_variables[cliques[i][j]]);
		}
		clusterSequence.push_back(clusterVars);
	}

	return clusterSequence;
}

/**
 * Display the region graph
 */
void outputRG() {
	for (int i = 0; i < (int)m_regions.size(); i++) {
		cout << "Region[" << i << "]:" << m_regions[i]->getCountingNumber() << ":" << m_regions[i]->toString() << endl;
	}
	cout << endl;
	for (int i = 0; i < (int)m_edges.size(); i++) {
		cout << "Edge[" << i << "]:" << m_edges[i]->toString() << endl;
	}
}


/**
 * Creates a region graph corresponding to the Bethe Approximation.
 */
void createRG_bethe() {
	m_regions.clear();
	m_edges.clear();

	// Sort the functions from largest to smallest
	vector<LogFunction*> functions_sorted = vector<LogFunction*>(m_functions);
	sort(functions_sorted.begin(), functions_sorted.end(), LogFunction::compareBySize);

	// Initialize the factor graph
	m_logPR_Partial = 0;

	// Keep track of the factors mentioning each variable
	vector<vector<int> > varIDToFactorIdxMap = vector<vector<int> >(m_numVariables);

	// Create a region for each non-trivial function
	for	(int i = 0; i < (int)functions_sorted.size(); i++) {
		// Don't create regions for empty factors
		if (functions_sorted[i]->getNumValues() > 1) {

			// Don't create a region for any functions that are subsumed by other functions
			bool addedToExistingFactor = false;
			for (int j = 0; j < (int)m_regions.size(); j++) {
				if (isSubsetVar(functions_sorted[i]->getScope(),m_regions[j]->getVariables())) {
					addedToExistingFactor = true;
					break;
				}
			}
			if (addedToExistingFactor) {
				continue;
			}

			// Create a new region
			Region* region = new Region(functions_sorted[i]->getScope());

			// Keep track of the variables in this factor
			for (int j = 0; j < (int)region->getVariables().size(); j++) {
				varIDToFactorIdxMap[region->getVariables()[j]->getID()].push_back((int)m_regions.size());
			}

			// Add region to the list of regions
			m_regions.push_back(region);
		}
	}

	// Create an inner region for each variable
	for (int i = 0; i < m_numVariables; i++) {

		// Ignore evidence variables and variables with a single parent
		if ((m_variables[i]->getDomainSize() > 1) && (((int)varIDToFactorIdxMap[i].size()) > 1)) {

			Region* region = new Region(m_variables[i]);
			m_regions.push_back(region);

			// Add edges from the 'covering' outer regions
			for (int j = 0; j < (int)varIDToFactorIdxMap[i].size(); j++) {
				Edge* edge = new Edge(m_regions[varIDToFactorIdxMap[i][j]], region);
				m_edges.push_back(edge);
			}
		}
	}
}

/**
 * Given a set of clusters <tt>A</tt>, produces a set of inner clusters <tt>B</tt>
 * using the cluster variation method. In particular every element in <tt>B</tt>,
 * is formed from the intersection of two clusters in <tt>A</tt> and all elements of
 * <tt>B</tt> that are subsets (or equal) of elements in <tt>B</tt> are removed.
 * @param clusters
 * @return int -- The number of clusters added
 */
int getInnerClusters_CVM(vector<vector<Variable*> >& clusters) {
	int totalAdded = 0;
	int numAdded = clusters.size();

	// Don't allow outer regions that are subsumed by outer regions
	for (int i = 0; i < (int)clusters.size(); i++) {
		for (int j = 0; j < (int)clusters.size(); j++) {
			if (i != j) {
				assert (!isSubsetVar(clusters[i], clusters[j]));
			}
		}
	}
	vector<vector<Variable*> > addedClusters;
	vector<vector<Variable*> > newClusters = clusters;

	while (numAdded > 0) {
		numAdded = 0;

		vector<vector<Variable*> > intersections;

		for (int i = 0; i < (int)newClusters.size(); i++) {
			for (int j = i + 1; j < (int)newClusters.size(); j++) {
				vector<Variable*> isect;
				getIntersectionVar(newClusters[i], newClusters[j], isect);

				if (!isect.empty()) {

					bool isRedundant = false;
					for (int k = 0; k < (int)addedClusters.size(); k++) {
						if (isEqualVar(addedClusters[k], isect)) {
							isRedundant = true;
							break;
						}
					}

					if (!isRedundant) {
						addedClusters.push_back(isect);
						numAdded++;
						intersections.push_back(isect);
					}
				}
			}
		}

		newClusters.clear();
		newClusters = intersections;

		totalAdded += numAdded;
	} // While numAdded > 0

	cout << "Added " << totalAdded << " inner regions!" << endl;

	sort(addedClusters.begin(), addedClusters.end(), compareBySizeVar_Descend);
	for (int j = 0; j < (int)addedClusters.size(); j++) {
		clusters.push_back(addedClusters[j]);
	}

	return totalAdded;
}

/**
 * Creates an RG where the outer regions are as specified and the
 * inner regions are constructed using the cluster variation method.
 */
void createRG_cvm(vector<vector<Variable*> > outerRegionVars) {
	m_regions.clear();
	m_edges.clear();

	// The set of variables defining each region
	vector<vector<Variable*> > regionsVar;

	// Create a region for each set of variables
	int numOuterRegions = outerRegionVars.size();
	for (int i = 0; i < numOuterRegions; i++) {
		Region* region = new Region(outerRegionVars[i]);
		m_regions.push_back(region);
		regionsVar.push_back(outerRegionVars[i]);
	}

	// Determine the complete set of regions
	int numInnerRegions = getInnerClusters_CVM(regionsVar);
	int numRegions = numInnerRegions + numOuterRegions;
	assert(numRegions == (int)regionsVar.size());

	// Create regions for all clusters of variables found by the cvm
	for (int i = numOuterRegions; i < numRegions; i++) {
		Region* region = new Region(regionsVar[i]);
		m_regions.push_back(region);

		for (int j = 0; j < numOuterRegions; j++) {
			if (isSubsetVar(region->getVariables(), m_regions[j]->getVariables())) {
				Edge* edge = new Edge(m_regions[j], region);
				m_edges.push_back(edge);
			}
		}
	}
}

/**
 * Reads a region graph specification file from disk.
 *
 * NOTE:  Assumes that the uai network file has already been read/processed.
 *
 * @param rgFileName - Name of the file containing the network
 * @param outerRegions - The outer regions specified in the file
 */
void createRGFromFile(const char* rgFileName) {
	cout << "Reading region graph from file: " << rgFileName << endl;

	ifstream infile(rgFileName);

	if (!infile.is_open()) {
		cout << "Failed to open file:" << rgFileName << endl;
		exit(-1);
	}

	// Get the number of outer regions
	int numOuterRegions;
	infile >> numOuterRegions;

	// Get the variables in each outer region
	vector<vector<Variable*> > outerRegionVars(numOuterRegions);
	for (int i = 0; i < numOuterRegions; i++) {
		int numVarsInRegion;
		infile >> numVarsInRegion;
		outerRegionVars[i] = vector<Variable*> (numVarsInRegion);
		for (int j = 0; j < numVarsInRegion; j++) {
			int varID;
			infile >> varID;
			assert(varID >= 0);
			assert(varID < m_numVariables);
			outerRegionVars[i][j] = m_variables[varID];
		}
	}

	// Add outer regions for each function not covered by some region in rg
	vector<LogFunction*> functions_sorted = vector<LogFunction*>(m_functions);
	sort(functions_sorted.begin(), functions_sorted.end(), LogFunction::compareBySize);
	for (int i = 0; i < (int)functions_sorted.size(); i++) {
		bool isCovered = false;
		for (int j = 0; j < (int)outerRegionVars.size(); j++) {
			if (isSubsetVar(functions_sorted[i]->getScope(), outerRegionVars[j])) {
				isCovered = true;
				break;
			}
		}

		if (!isCovered) {
			vector<Variable*> outerRegion = vector<Variable*>(functions_sorted[i]->getScope());
			outerRegionVars.push_back(outerRegion);
		}
	}

	createRG_cvm(outerRegionVars);
}

/**
 * Constructs a RG using the Junction Graph method
 */
void buildJoinGraph() {
	// Build the join graph using the mini-bucket heuristic
	vector<int> order = getEliminationOrder();

	JoinGraph* jg = JoinGraph::runMiniBucketSchematic(order,m_variables, m_functions, m_iBound);

	// And simplify the join graph to remove any unnecessary nodes/edges
	jg->simplify();

	//
	// Create a region graph from the join graph as follows:
	// 	1) The outer regions in the RG are all the JGNodes; and
	//  2) The inner regions are all the JG Edges.
	//
	map<int,int> nodeIDToIndexMap;
	int numOuterRegions = jg->getNodes().size();
	for (int i = 0; i < numOuterRegions; i++) {
		JGNode* node_i = jg->getNodes()[i];

		// Update the map from node id to index
		nodeIDToIndexMap[node_i->getID()] = i;

		// Create an outer region for node i and add necessary original functions
		Region* newOuterRegion = new Region(node_i->getVariables());
		for (int j = 0; j < (int)node_i->getOriginalFunctions().size(); j++) {
			newOuterRegion->addFunction(node_i->getOriginalFunctions()[j]);
		}
		m_regions.push_back(newOuterRegion);
	}

	// Each inner region corresponds to an edge in the join graph. However,
	// since multiple edges in the join graph may be defined over the same
	// separators, each inner region may have more than two parents.
	for (int i = 0; i < (int)jg->getEdges().size(); i++) {
		JGEdge* edge_i = jg->getEdges()[i];

		// Determine if this inner region is redundant
		bool isRedundant = false;
		for (int j = numOuterRegions; j < (int)m_regions.size(); j++) {
			if (isEqualVar(edge_i->getVariables(),m_regions[j]->getVariables())) {
				isRedundant = true;

				// Determine if node i or node j, is the new parent outer region
				bool found_node_i = false;
				bool found_node_j = false;
				for (int k = 0; k < (int)m_regions[j]->getEdges().size(); k++) {
					// If an edge already exists from the region corresponding to node i to
					// the redundant inner region
					if (isEqualVar(edge_i->getNode_i()->getVariables(), m_regions[j]->getEdges()[k]->getParent()->getVariables())) {
						found_node_i = true;
						break;
					}
					// If an edge already exists from the region corresponding to node j to
					// the redundant inner region
					if (isEqualVar(edge_i->getNode_j()->getVariables(), m_regions[j]->getEdges()[k]->getParent()->getVariables())) {
						found_node_j = true;
						break;
					}
				}

				assert(found_node_i || found_node_j);

				// If the region corresponding to node i was redundant, the region
				// corresponding to node j is new
				if (found_node_i) {
					Edge* newEdge = new Edge(m_regions[nodeIDToIndexMap[edge_i->getNode_j()->getID()]], m_regions[j]);
					m_edges.push_back(newEdge);
				}
				// Otherwise the region corresponding to node j was redundant and the region
				// corresponding to node i is new
				else {
					Edge* newEdge = new Edge(m_regions[nodeIDToIndexMap[edge_i->getNode_i()->getID()]], m_regions[j]);
					m_edges.push_back(newEdge);
				}
				break;
			}
		}

		if (isRedundant) {
			continue;
		}

		//
		// Create an inner region for edge i.
		//
		Region* newInnerRegion = new Region(edge_i->getVariables());
		m_regions.push_back(newInnerRegion);

		// Add the necessary edges to this region
		Edge* newRGEdge_i = new Edge(m_regions[nodeIDToIndexMap[edge_i->getNode_i()->getID()]], newInnerRegion);
		m_edges.push_back(newRGEdge_i);

		Edge* newRGEdge_j = new Edge(m_regions[nodeIDToIndexMap[edge_i->getNode_j()->getID()]], newInnerRegion);
		m_edges.push_back(newRGEdge_j);
	}

	// Place the regions in descending order by size
	sort(m_regions.begin(), m_regions.end(), Region::compareRegionBySize_Descend);

	// Specify the ancestors of each region...but set the counting numbers using the join
	// graph structure.
	Region::determineAncestorsRG(m_regions);

	// Set the counting numbers
	for (int i = 0; i < (int)m_regions.size(); i++) {
		// If this is an isolated region, give it counting number 1
		if (m_regions[i]->getEdges().empty()) {
			m_regions[i]->setCountingNumber(1);
		}
		else {
			// If this region is a parent...
			if (m_regions[i]->getEdges()[0]->getParent()->getID() == m_regions[i]->getID()) {
				m_regions[i]->setCountingNumber(1);
			}
			else {
				m_regions[i]->setCountingNumber(1 - ((int)m_regions[i]->getEdges().size()));
			}
		}
	}
}

/**
 * Verifies that the region graph is 1-connected and 1-balanced
 */
void verifyRG() {

	for (int i = 0; i < m_numVariables; i++) {

		if (m_variables[i]->getDomainSize() == 1) {
			continue;
		}

		// Get the set of regions containing this variable
		vector<Region*> involvedRegions;
		map<int,int> regionIDToIndexMap;
		int sumCountingNums = 0;
		for (int j = 0; j < (int)m_regions.size(); j++) {
			if (containsVar(m_regions[j]->getVariables(), m_variables[i])) {
				regionIDToIndexMap[m_regions[j]->getID()] = involvedRegions.size();
				involvedRegions.push_back(m_regions[j]);
				sumCountingNums += m_regions[j]->getCountingNumber();
			}
		}

		// Every variable must appear in some region
		assert(((int)involvedRegions.size()) > 0);

		if (sumCountingNums != 1) {
			cout << "RG is not 1-balanced for variable " << m_variables[i]->getID() << endl;

			cout << "Variable " << m_variables[i]->getID() << " in:" << endl;
			for (int j = 0; j < (int)involvedRegions.size(); j++) {
				cout << "  " << involvedRegions[j]->toString() << " c:" << involvedRegions[j]->getCountingNumber() << endl;
			}

		}
		assert(sumCountingNums == 1);

		// Traverse the set of involved regions
		queue<Region*> fringe;
		fringe.push(involvedRegions[0]);
		vector<bool> visited = vector<bool>(involvedRegions.size());
		fill(visited.begin(),visited.end(),false);
		while(!fringe.empty()) {
			Region* next = fringe.front();
			fringe.pop();

			if (visited[regionIDToIndexMap[next->getID()]]) {
				continue;
			}

			visited[regionIDToIndexMap[next->getID()]] = true;

			for (int j = 0; j < (int)next->getEdges().size(); j++) {
				if (containsVar(next->getEdges()[j]->getChild()->getVariables(), m_variables[i])) {

					if (next->getEdges()[j]->getParent() == next) {
						if (!visited[regionIDToIndexMap[next->getEdges()[j]->getChild()->getID()]]) {
							fringe.push(next->getEdges()[j]->getChild());
						}
					}
					else {
						assert(next->getEdges()[j]->getChild() == next);

						if (!visited[regionIDToIndexMap[next->getEdges()[j]->getParent()->getID()]]) {
							fringe.push(next->getEdges()[j]->getParent());
						}
					}
				}
			}
		}

		if (((int)involvedRegions.size()) == 1) {
			continue;
		}

		// Ensure that all regions have been visited
		for (int j = 0; j < (int)visited.size(); j++) {
			if (!visited[j]) {
				cout << "Region[" << involvedRegions[j]->getID() << "] " << varsToString(involvedRegions[j]->getVariables()) << endl;
				cout << "not visited for variable " << i << endl;
			}
			assert(visited[j]);
		}
	}
}

/**
 * Outputs the current estimate of logPR given the beliefs at each region in the
 * current region graph. This method assumes that the beliefs at each region have
 * been computed, so that the region contains valid energy and entropies.
 */
void outputPR(string outFileName, double log10PR) {
	cout << "log10PR: " << log10PR << endl;
	FILE *pFile = fopen(outFileName.c_str(), "w");
	if (pFile == NULL) {
		cout << "Failed to create handle for file " << outFileName.c_str() << endl;
		exit(-1);
	}
	fprintf(pFile,"%0.15e\n",log10PR);
	fclose(pFile);
}

/**
 * Outputs the current estimate of the marginals.
 */
void outputBEL(string outFileName) {
	FILE *pFile = fopen(outFileName.c_str(), "w");
	if (pFile == NULL) {
		cout << "Failed to create handle for file " << outFileName.c_str() << endl;
		exit(-1);
	}

	// Output the number of beliefs
	fprintf(pFile,"%d\n",m_numFunctions);
	for (int i = 0; i < m_numFunctions; i++) {
		// Output the number of elements in this belief
		fprintf(pFile,"%d ",m_functions[i]->getNumValues());

		// Find smallest region containing this function, keeping in mind
		// that regions are placed in descending order by size.
		for (int j = (int)m_regions.size() - 1; j >= 0; j--) {
			if (isSubsetVar(m_functions[i]->getScope(), m_regions[j]->getVariables())) {
				// Determine the set of variables to be eliminated
				vector<Variable*> elim;
				getDifferenceVar(m_regions[j]->getVariables(), m_functions[i]->getScope(), elim);

				// Compute the marginal belief
				LogFunction* temp = m_regions[j]->getBelief()->marginalize(elim);

				// Output the marginal in an order consistent with the variable
				// ordering specified in the input UAI file
				assert(temp->getNumValues() == m_functions[i]->getNumValues());
				for (int k = 0; k < temp->getNumValues(); k++) {

					Variable::setAddress(m_functions[i]->getOriginalScope(), k);
					int idx = Variable::getAddress(temp->getScope());

					fprintf(pFile,"%0.10e ", temp->getValueAt(idx));
				}

				// Clean up
				delete temp;

				// Compute the next marginal
				break;
			}
		}
		fprintf(pFile,"\n");
	}

	fclose(pFile);
}

/**
 * The main computational loop.
 */
double runSingleLoop(bool &converged, int & iter, int maxIters, double convergenceThreshold, bool verbose = true, bool initialize = true) {

	// Initialize each of the regions and determine the set of inner regions
	vector<Region*> innerRegions;
	for (int i = 0; i < (int)m_regions.size(); i++) {
		if (initialize) {
			m_regions[i]->initialize();
		}

		if (!m_regions[i]->getAncestors().empty()) {
			innerRegions.push_back(m_regions[i]);
		}
	}

	// Initialize each of the edges
	if (initialize) {
		for (int i = 0; i < (int)m_edges.size(); i++) {
			m_edges[i]->initialize();
		}
	}

	//
	// Main computational loop
	//
	iter = 0;
	converged = false;
	if (verbose) {
		cout << "maxIters:" << maxIters << " alpha:" << m_alpha << endl;
	}
	int numInnerRegions = (int)innerRegions.size();
	vector<int> order = m_random->getRandomPerm(numInnerRegions);
	double maxDiff = 0;
	while (!converged && (iter < maxIters)) {

		m_random->shuffle(order);

		// For each inner region
		for (int i = 0; i < numInnerRegions; i++) {

			// Update each parent-to-child message
			for (int j = 0; j < (int)innerRegions[order[i]]->getEdges().size(); j++) {
				innerRegions[order[i]]->getEdges()[j]->updateParentToChildMsg();
			}

			// Update the child belief
			innerRegions[order[i]]->updateBelief(m_alpha);

			// Update each parent of this child
			for (int j = 0; j < (int)innerRegions[order[i]]->getEdges().size(); j++) {
				innerRegions[order[i]]->getEdges()[j]->updateChildToParentMsg();
				innerRegions[order[i]]->getEdges()[j]->getParent()->updateBelief(m_alpha);
			}
		}

		// Check convergence
		maxDiff = 0;
		for (int i = 0; i < (int)m_regions.size(); i++) {
			maxDiff = fmax(maxDiff, m_regions[i]->getDifferenceLinf());
		}
		converged = false;
		if (maxDiff < convergenceThreshold) {
			converged = true;
		}
		if (verbose) {
			cout << "  iter " << iter << " maxDiff:" << maxDiff << endl;
		}
		iter++;

	} // end while !converged

	if (converged && verbose) {
		cout << "Converged after " << iter << " iterations!" << endl;
	}

	// Return the max difference
	return maxDiff;
}

void runDoubleLoop(bool &converged, int &iter, bool verbose = true) {

	// Number of regions
	int numRegions = (int)m_regions.size();

	// Save the original counting numbers of the inner regions
	vector<int> countingNumbers;
	vector<Region*> innerRegions;
	for (int i = 0; i < numRegions; i++) {
		// If this is an inner region
		if (!m_regions[i]->getAncestors().empty()) {
			innerRegions.push_back(m_regions[i]);
			countingNumbers.push_back(m_regions[i]->getCountingNumber());
		}
	}

	// Keep region beliefs to check convergence
	vector<LogFunction*> oldBeliefs = vector<LogFunction*>(numRegions);
	vector<Region*> outerRegions;
	for (int i = 0; i < numRegions; i++) {
		m_regions[i]->initialize();
		oldBeliefs[i] = new LogFunction(*m_regions[i]->getBelief());

		if (m_regions[i]->getAncestors().empty()) {
			outerRegions.push_back(m_regions[i]);
		}
	}

	// Initialize the edges
	for (int i = 0; i < (int)m_edges.size(); i++) {
		m_edges[i]->initialize();
	}

	//
	// Main computational loop
	//
	iter = 0;
	converged = false;
	int numOuterRegions = (int)outerRegions.size();
	m_alpha = 1.0; // Don't damp in double-loop and run to convergence...
	while (!converged && (iter < m_maxIters)) {

		// Restore the original counting numbers
		for (int i = 0; i < (int)innerRegions.size(); i++) {
			innerRegions[i]->setCountingNumber(countingNumbers[i]);
		}

		// Update the belief at each outer region
		if (iter > 0) {
			for (int i = 0; i < numOuterRegions; i++) {
				outerRegions[i]->updateBaseFactor();
			}
		}

		// Bound the inner regions
		for (int i = 0; i < (int)innerRegions.size(); i++) {
			innerRegions[i]->setCountingNumber(0);
		}

		// Run the inner loop for a few iterations
		bool innerLoopConverged = false;
		int innerLoopIters = 0;
		runSingleLoop(innerLoopConverged, innerLoopIters, m_innerLoopMaxIters, 0.01*m_convergenceThreshold, false, false);

		// Check convergence
		double maxDiff = 0;
		for (int i = 0; i < numRegions; i++) {
			maxDiff = fmax(maxDiff, LogFunction::getDistanceLinf(oldBeliefs[i], m_regions[i]->getBelief()));
			oldBeliefs[i]->assign(*m_regions[i]->getBelief());
		}

		converged = false;
		if (maxDiff < m_convergenceThreshold) {
			converged = true;
		}

		if (verbose) {
			cout << "iter " << iter << " maxDiff:" << maxDiff << endl;
		}
		iter++;
	}

	if (converged) {
		cout << "Double-Loop converged after " << iter << " iterations" << endl;
	}

	// Restore the original counting numbers
	for (int i = 0; i < (int)innerRegions.size(); i++) {
		innerRegions[i]->setCountingNumber(countingNumbers[i]);
	}

	// Clean-up after
	for (int i = 0; i < numRegions; i++) {
		delete oldBeliefs[i];
	}
	oldBeliefs.clear();
}
/**
 * Assign non-empty functions to each region in the region graph
 */
void assignFunctions() {

	m_logPR_Partial = 0;

	for (int i = 0; i < (int)m_functions.size(); i++) {
		if (m_functions[i]->getNumValues() == 1) {
			m_logPR_Partial += m_functions[i]->getLogValueAt(0);
		}
		else {
			bool isFunctionAssigned = false;
			for (int j = 0; j < (int)m_regions.size(); j++) {
				// If this is an outer region that subsumes this function
				if (m_regions[j]->getAncestors().empty() &&
					(m_functions[i]->getScope().size() <= m_regions[j]->getVariables().size()) &&
					isSubsetVar(m_functions[i]->getScope(), m_regions[j]->getVariables())) {
					m_regions[j]->addFunction(m_functions[i]);
					isFunctionAssigned = true;
					break;
				}
			}
			if (!isFunctionAssigned) {
				cout << "Failed to assign function " << varsToString(m_functions[i]->getScope()) << endl;
			}
			assert(isFunctionAssigned);
		}
	}
}

/**
 * Removes the specified edge from the current region graph
 */
bool removeEdge(Edge* edge) {
	edge->getChild()->removeEdge(edge);
	edge->getParent()->removeEdge(edge);

	int idx = -1;
	for (int i = 0; i < (int)m_edges.size(); i++) {
		if (m_edges[i]->getID() == edge->getID()) {
			idx = i;
			break;
		}
	}
	if (idx >= 0) {
		m_edges.erase(m_edges.begin() + idx);
		return true;
	}

	return false;
}

/**
 * Simplifies the region graph.
 */
void simplifyRG() {

	int numEdgesRemoved = 0;
	int numRegionsRemoved = 0;

	sort(m_regions.begin(), m_regions.end(), Region::compareRegionBySize_Descend);

	Region::determineAncestorsRG(m_regions);

	// Search for redundant constraints as follows:
	//  For each region \beta, look for ancestors \beta' of
	//  \beta such that \beta and \beta' have two or more
	//  outer regions in common
	for (int i = (int)m_regions.size() - 1; i >= 0; i--) {

		// For each region \beta' that is ancestor of region \beta
		for (int j = 0; j < (int)m_regions[i]->getAncestors().size(); j++) {

			// Determine the ancestral regions that \beta and \beta' have in common
			vector<Region*> isect;
			Region::getIntersectionRegion(m_regions[i]->getAncestors(), m_regions[i]->getAncestors()[j]->getAncestors(), isect);

			// Determine which of the common ancestral regions are outer regions
			vector<Region*> outerCoveringRegions;
			for (int k = 0; k < (int)isect.size(); k++) {
				// This is an outer region
				if (isect[k]->getAncestors().empty()) {
					outerCoveringRegions.push_back(isect[k]);
				}
			}

			// If the set of common ancestral outer regions is of size 2 or more,
			// then some of the edges are redundant
			for (int k = 1; k < (int)outerCoveringRegions.size(); k++) {

				// Remove edge from common outer region to the inner region
				for (int l = 0; l < (int)outerCoveringRegions[k]->getEdges().size(); l++) {
					Edge* edge = outerCoveringRegions[k]->getEdges()[l];
					if (edge->getChild()->getID() == m_regions[i]->getID()) {

						if (removeEdge(edge)) {
							if (m_isVerbose) {
								cout << "Removing edge " << edge->toString() << endl;
							}
							numEdgesRemoved++;
							delete edge;
						}
						break;
					}
				}
			}
		}
	}

	// Remove unnecessary inner regions. An unnecessary inner region has only
	// 1 parent and counting number zero
	for (int i = (int)m_regions.size() - 1; i >= 0; i--) {

		// If this is an inner region, with counting number zero and
		// only 1 parent it can be removed
		Region* region_i = m_regions[i];
		if (!region_i->getAncestors().empty() &&
			(region_i->getCountingNumber() == 0) &&
			(((int)region_i->getEdges().size()) == 1)) {

			numRegionsRemoved++;

			// Remove this region and its lone edge
			Edge* edge = region_i->getEdges()[0];
			removeEdge(edge);

			if (m_isVerbose) {
				cout << "Removing region " << region_i->toString() << " and edge " << edge->toString() << endl;
			}

			delete edge;

			// Remove this region from the ancestor lists its descendants
			for (int j = 0; j < (int)m_regions.size(); j++) {
				if ((j != i) && isSubsetVar(m_regions[j]->getVariables(), m_regions[i]->getVariables())) {
					m_regions[j]->removeAncestor(m_regions[i]);
				}
			}

			delete m_regions[i];
			m_regions.erase(m_regions.begin()+i);
		}
	}

	// Look for inner regions \beta such that c_{\beta} + n_{\beta} = 0,
	// where c_{\beta} is the counting number of region \beta and
	// n_{\beta} is the number of outer regions covering region \beta
	for (int i = 0; i < (int)m_regions.size(); i++) {
		if (!m_regions[i]->getAncestors().empty() &&
			(m_regions[i]->getCountingNumber() + ((int)m_regions[i]->getEdges().size()) == 0)) {

			// If we find such a region, add an edge to a 'new' covering outer region
			bool isAddEdge = false;
			vector<Region*> ancestors_i = vector<Region*>(m_regions[i]->getAncestors());
			sort(ancestors_i.begin(), ancestors_i.end(), Region::compareRegionBySize_Ascend);
			for (int j = 0; j < (int)ancestors_i.size(); j++) {

				// If region j is an outer region and there isn't already an
				// edge between region i and region j...
				isAddEdge = true;
				Region* region_j = ancestors_i[j];
				if (region_j->getAncestors().empty() && (!region_j->getVariables().empty())) {

					for (int k = 0; k < (int)m_regions[i]->getEdges().size(); k++) {
						if (region_j->getID() == m_regions[i]->getEdges()[k]->getParent()->getID()) {
							isAddEdge = false;
							break;
						}
					}
				}
				else {
					isAddEdge = false;
				}

				// Add an edge between region i and outer region j
				if (isAddEdge) {
					if (m_isVerbose) {
						cout << "Adding edge to fix coveringNumber:" << region_j->toString() << "-" << m_regions[i]->toString() << endl;
					}
					Edge* edge = new Edge(region_j, m_regions[i]);
					m_edges.push_back(edge);
					numEdgesRemoved--;
					break;
				}
			} // end for each ancestor of region i
			assert(isAddEdge);
		} // If inner region such that c_{\beta} + n_{\beta} = 0
	}// For each region i

	if ((numEdgesRemoved > 0) || (numRegionsRemoved > 0)) {
		cout << "Removed " << numEdgesRemoved << " Edges and " << numRegionsRemoved << " Regions!" << endl;
	}
}

/**
 * Clears the current region graph
 */
void clearRG() {
	for (int i = 0; i < (int)m_regions.size(); i++) {
		delete m_regions[i];
	}
	m_regions.clear();
	for (int i = 0; i < (int)m_edges.size(); i++) {
		delete m_edges[i];
	}
	m_edges.clear();

	Region::resetIDs();
	Edge::resetIDs();
}

/**
 * The main method
 */
int main(int argc, char* argv[]) {

	// Record the start time
	clock_t startTime = clock();

	// Read the parameters
	readParameters(argc, argv);

	// Read the uai file
	readUAI(m_uaiFileName.c_str());

	// Build the undirected graph corresponding to this MN
	buildGraph();

	// Create the Region Graph
	if (m_isFromFile) {
		createRGFromFile(m_rgFileName.c_str());
		simplifyRG();
		assignFunctions();
	}
	else {
		if (m_rgType == BETHE) {
			if (m_isCVM) {
				cout << "Running GBP w/ cluster variation on original functions." << endl;

				// Get a cluster collection of the original functions
				ClusterSet* clusterSet = new ClusterSet(m_functions);

				// Create the region graph given this clustering
				vector<vector<Variable*> > outerRegionVars = clusterSet->getClusterVars();
				createRG_cvm(outerRegionVars);
			}
			else {
				cout << "Running BP" << endl;
				createRG_bethe();
			}

			// Simplify the region graph before running message passing
			simplifyRG();

			// Assign functions to the regions
			assignFunctions();
		}
		else if (m_rgType == IJGP) {

			if (m_isCVM) {

				cout << "Running GBP w/ CVM using mini-bucket schema on iBound=" << m_iBound << endl;

				// Get a cluster collection of the original functions
				ClusterSet* clusterSet = new ClusterSet(m_functions);

				// Get an elimination order
				vector<int> order = getEliminationOrder();

				// Get a set of clusters by running the mini-bucket schematic
				clusterSet->runMiniBucket(order, m_variables, m_iBound);

				// Create the region graph given this clustering
				vector<vector<Variable*> > outerRegionVars = clusterSet->getClusterVars();
				createRG_cvm(outerRegionVars);

				// Simplify the region graph before running message passing
				simplifyRG();

				// Assign functions to the regions
				assignFunctions();
			}
			else {
				cout << "Running IJGP w/ iBound=" << m_iBound << endl;

				// Build the join graph with the specified iBound
				buildJoinGraph();
			}
		}

	} // if create rg from file

	// Validate the region graph before running message passing
	verifyRG();

	if (m_isVerbose) {
		outputRG();
		int sumCnt = 0;
		for (int i = 0; i < (int)m_regions.size(); i++) {
			sumCnt += m_regions[i]->getCountingNumber();
		}
		cout << "Sum of Counting Numbers:" << sumCnt << endl;
	}

	//
	// Run message passing on the specified region graph
	//
	bool converged = false;
	int iters = 0;
	if (m_isDoubleLoop) {
		runDoubleLoop(converged, iters, m_isVerbose);
	}
	else {
		runSingleLoop(converged, iters, m_maxIters, m_convergenceThreshold, m_isVerbose);
	}

	//
	// Compute and output both BEL and PR
	//
	string PRFileName = m_uaiFileName + ".PR";
	string BELFileName = m_uaiFileName + ".BEL";
	outputPR(PRFileName, getLogPR()/log(10));
	outputBEL(BELFileName);

	//
	// Output summary information
	//
	string summaryFileName = m_uaiFileName + ".summary";
	double totalExecTime = (((double)clock() - startTime)/CLOCKS_PER_SEC);
	FILE *pFile = fopen(summaryFileName.c_str(), "w");
	if (pFile == NULL) {
		cout << "Failed to create handle for file " << summaryFileName.c_str() << endl;
		exit(-1);
	}
	fprintf(pFile,"%d %d %0.3f\n",converged,iters,totalExecTime);
	fclose(pFile);

	// Normal exit
	cout << "Total Execution Time(sec): " << totalExecTime << endl;
	return 0;
}
