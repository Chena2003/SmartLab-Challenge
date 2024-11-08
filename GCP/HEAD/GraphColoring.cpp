#include "GraphColoring.h"
#include "doublevector.h"
#include "staticvector.h"
#include "head.h"

#include <random>
#include <iostream>
#include <vector>


using namespace std;


namespace szx {

	class Solver {
		// random number generator.
		mt19937 pseudoRandNumGen;
		void initRand(int seed) { pseudoRandNumGen = mt19937(seed); }
		int fastRand(int lb, int ub) { return (pseudoRandNumGen() % (ub - lb)) + lb; }
		int fastRand(int ub) { return pseudoRandNumGen() % ub; }
		int rand(int lb, int ub) { return uniform_int_distribution<int>(lb, ub - 1)(pseudoRandNumGen); }
		int rand(int ub) { return uniform_int_distribution<int>(0, ub - 1)(pseudoRandNumGen); }

		const int TabuIters[12] = { 3000, 3000, 5000, 5000, 6000, 6000, 8000, 8000, 25000, 7000, 40000, 30000 };
		int initTabuIter(int nodeNum, int edgeNum) {
			int i, j;
			switch (nodeNum) {
			case 125: i = 0; break;
			case 250: i = 1; break;
			case 500: i = 2; break;
			case 1000: i = 3; break;
			}
			int p1 = 0.1 * nodeNum * nodeNum;
			int p2 = 3 * p1;

			if (edgeNum < p1) j = 0;
			else if (edgeNum > p2) j = 2;
			else j = 1;

			return TabuIters[3 * i + j];
		};

	public:
		void solve(NodeColors& output, GraphColoring& input, function<long long()> restMilliSec, int seed) {

			int TabuIter = initTabuIter(input.nodeNum, input.edgeNum); // ��ʼ����������
			head solver(input.nodeNum, input.edgeNum, input.colorNum, TabuIter, input.edges, seed);
			solver.HybridEvolution();
			solver.WriteSolution(output);
		};
	};
		// solver.
	void solveGraphColoring(NodeColors& output, GraphColoring& input, function<long long()> restMilliSec, int seed) {
		Solver().solve(output, input, restMilliSec, seed);
	};
}
