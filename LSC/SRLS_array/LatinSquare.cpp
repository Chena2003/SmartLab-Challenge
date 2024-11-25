#include "LatinSquare.h"

#include <random>
#include <iostream>
#include "work.h"
#include <chrono>

using namespace std;
using namespace std::chrono;


namespace szx {

class Solver {
	// random number generator.
	mt19937 pseudoRandNumGen;
	void initRand(int seed) { pseudoRandNumGen = mt19937(seed); }
	int fastRand(int lb, int ub) { return (pseudoRandNumGen() % (ub - lb)) + lb; }
	int fastRand(int ub) { return pseudoRandNumGen() % ub; }
	int rand(int lb, int ub) { return uniform_int_distribution<int>(lb, ub - 1)(pseudoRandNumGen); }
	int rand(int ub) { return uniform_int_distribution<int>(0, ub - 1)(pseudoRandNumGen); }

public:
	void solve(ostream& outputStream, LatinSquare& input, function<long long()> restMilliSec, int seed) {
		// Debug
		// steady_clock::time_point StartTime = steady_clock::now();
		Work solver(input.n * input.n, input.n, input.fixedNums, seed);

		int iter = solver.solve(restMilliSec);
		solver.savesol(outputStream);
		// steady_clock::time_point EndTime = steady_clock::now();
		// uint64_t duration = duration_cast<milliseconds>(EndTime - StartTime).count();

		// cerr << "Run Time: " << duration << " Frequency: " << iter << endl;
		// solver.check();
	}
};

// solver.
void solveLatinSquare(ostream& outputStream, LatinSquare& input, function<long long()> restMilliSec, int seed) {
	Solver().solve(outputStream, input, restMilliSec, seed);
}

}
