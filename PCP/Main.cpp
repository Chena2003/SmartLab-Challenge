#include <iostream>
#include <string>
#include <chrono>
#include <fstream>

#include "PCenter.h"

// #define __ENABLE_DEBUG__ 

using namespace std;
using namespace std::chrono;
using namespace szx;


void loadInput(istream& is, PCenter& pc) {
	is >> pc.nodeNum >> pc.centerNum;

	pc.coverages.resize(pc.nodeNum);
	pc.coveredNodeNums.resize(pc.nodeNum);
	pc.serives.assign(pc.nodeNum, vector<bool>(pc.nodeNum, 0));
	pc.fixNodes.assign(pc.nodeNum, 0);

	int t = 0;
	for (int i = 0; i < pc.nodeNum; i++) {
		NodeId coveredNodeNum;
		is >> coveredNodeNum;
		pc.coveredNodeNums[i] = coveredNodeNum;
		pc.coverages[i].resize(coveredNodeNum);
		for (int j = 0; j < coveredNodeNum; j++) {
			is >> pc.coverages[i][j];
			pc.serives[0][0] = true;
			pc.serives[pc.coverages[i][j]][i] = true;
		}

		if (coveredNodeNum == 1 && pc.coverages[i][0] == i) {
			pc.fixNodes[i] = true;
			t++;
		}
	}
	pc.fixNum = t;

	EdgeId minEdgeLenRank;
	EdgeId maxEdgeLenRank;
	is >> maxEdgeLenRank >> minEdgeLenRank;
	pc.nodesWithDrops.resize(maxEdgeLenRank - minEdgeLenRank);
	for (auto r = pc.nodesWithDrops.begin(); r != pc.nodesWithDrops.end(); ++r) {
		NodeId nodeNumToDrop;
		is >> nodeNumToDrop;
		r->resize(nodeNumToDrop);
		for (auto node = r->begin(); node != r->end(); ++node) { is >> *node; }
	}
}

void saveOutput(ostream& os, Centers& centers) {
	for (auto center = centers.begin(); center != centers.end(); ++center) { os << *center << endl; }
}

void test(istream& inputStream, ostream& outputStream, long long secTimeout, int randSeed) {
	cerr << "load input." << endl;
	PCenter pc;
	loadInput(inputStream, pc);

	cerr << "solve." << endl;
	steady_clock::time_point start = steady_clock::now();
	steady_clock::time_point endTime = steady_clock::now() + seconds(secTimeout);
	Centers centers(pc.centerNum);
	solvePCenter(centers, pc, [&]() { return duration_cast<milliseconds>(endTime - steady_clock::now()).count(); }, randSeed);
	steady_clock::time_point end = steady_clock::now();

	cerr << "save output." << endl;
	saveOutput(outputStream, centers);

#ifdef __ENABLE_DEBUG__
	// 输出运行时间
	auto duration = duration_cast<milliseconds>(end - start);
	cerr << "Execution time: " << duration.count() << "ms" << endl;

	// 输出结果节点
	cerr << "Solution centers: " << endl;
	for (NodeId i : centers) {
		cerr << i << ' ';
	}
	cerr << endl;

	// 输出未覆盖节点
	cerr << "Uncovered centers: " << endl;
	vector<bool> flag(pc.nodeNum, 0);
	for (NodeId i : centers) 
		for (NodeId j : pc.coverages[i]) 
			flag[j] = true;

	for (NodeId i = 0; i < pc.nodeNum; i++) {
		if (!flag[i]) {
			cerr << i << endl;
		}
	}
#endif

}

void test(istream& inputStream, ostream& outputStream, long long secTimeout) {
	return test(inputStream, outputStream, secTimeout, static_cast<int>(time(nullptr) + clock()));
}

int main(int argc, char* argv[]) {
	cerr << "load environment." << endl;
	if (argc > 2) {
		long long secTimeout = atoll(argv[1]);
		int randSeed = atoi(argv[2]);
		test(cin, cout, secTimeout, randSeed);
	}
	else {
		ifstream ifs("instance/input.txt");
		ofstream ofs("instance/solution.txt");
		test(ifs, ofs, 1000000); // for self-test.
	}
	return 0;
}
