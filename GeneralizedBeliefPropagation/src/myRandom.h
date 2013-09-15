#ifndef MY_RANDOM_H_
#define MY_RANDOM_H_

#include "randomc.h"

using namespace std;

struct myRandom {
	int seed;
	CRandomMersenne RanGen;

	/**
	 * Default constructor.  Uses a default seed set to current time()
	 */
	myRandom() {
		seed=time(NULL);
		RanGen = CRandomMersenne(seed);
	}

	/**
	 * Default constructor.  Uses specified int as seed
	 */
	myRandom(int seed_) {
		seed = seed_;
		RanGen = CRandomMersenne(seed);
	}

	/**
	 * Sets the see used by the random number generator
	 */
	void setSeed(int seed_) {
		seed = seed_;
		RanGen = CRandomMersenne(seed);
	}

	/**
	 * Draws a random floating point number in [double.MIN_VALUE,double.MAX_VALUE)
	 */
	double getDouble() {
		return RanGen.Random();
	}

	/**
	 * Draws a random integer in [int.MIN_VALUE,int.MAX_VALUE)
	 */
	int getInt() {
		return RanGen.BRandom();
	}

	/**
	 * Draws a value in [0,max_value)
	 */
	int getInt(int max_value) {
		return RanGen.IRandomX(0,max_value-1);
	}

	/**
	 * Returns a random permutation of integers in [0,...,N).
	 */
	vector<int> getRandomPerm(int N) {
		vector<int> perm;
		for (int i = 0; i < N; i++) {
			perm.push_back(i);
		}

		for (int i = 0; i < N; i++) {
			int j = getInt(N-i);

			// Swap i and j
			int k = perm[i];
			perm[i] = perm[i+j];
			perm[i+j] = k;
		}

		return perm;
	}

	/**
	 * Randomly permutes the specified set of integers
	 */
	void shuffle(vector<int>& perm) {
		int N = (int)perm.size();
		for (int i = 0; i < N; i++) {
			int j = getInt(N-i);

			// Swap i and j
			int k = perm[i];
			perm[i] = perm[i+j];
			perm[i+j] = k;
		}
	}
};
#endif
