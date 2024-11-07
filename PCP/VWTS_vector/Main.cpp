#include <iostream>
#include <string>
#include <chrono>
#include <fstream>

#include "PCenter.h"


using namespace std;
using namespace std::chrono;
using namespace szx;


static const int ObjTable[15] = { 0, 108039, 50812, 32694, 24102, 18847, 9076, 5738, 4153, 3154, 2774, 2295, 1945, 1604, 1449 };
int find(int n) {
	int i;
	if (n <= 50) {
		i = n / 10;
	}
	else {
		i = n / 50 + 4;
	}

	return ObjTable[i];
}

void loadInput(istream& is, PCenter& pc) {
	ios::sync_with_stdio(false); // 禁用同步
	cin.tie(0);
	cout.tie(0);

	is >> pc.nodeNum >> pc.centerNum;

	// 直接构造具有所需大小的向量
	pc.coverages = vector<vector<NodeId>>(pc.nodeNum);
	pc.coveredNodeNums = vector<NodeId>(pc.nodeNum);
	pc.fixNodes.Nodes = vector<bool>(pc.nodeNum, false);

	int t = 0;
	for (int i = 0; i < pc.nodeNum; i++) {
		NodeId coveredNodeNum;
		is >> coveredNodeNum;
		pc.coveredNodeNums[i] = coveredNodeNum;
		pc.coverages[i] = vector<NodeId>(coveredNodeNum);
		for (int j = 0; j < coveredNodeNum; j++) {
			is >> pc.coverages[i][j];
		}

		if (coveredNodeNum == 1 && pc.coverages[i][0] == i) {
			pc.fixNodes.Nodes[i] = true;
			t++;
		}
	}
	pc.fixNodes.Num = t;

	EdgeId minEdgeLenRank;
	EdgeId maxEdgeLenRank;
	is >> maxEdgeLenRank >> minEdgeLenRank;
	if (maxEdgeLenRank == 0) {
		pc.obj = -1;
		return;
	}
	int obj = find(pc.centerNum);
	pc.L = minEdgeLenRank;
	pc.U = maxEdgeLenRank;
	pc.obj = obj;
	pc.nodesWithDrops = vector<vector<NodeId>>(maxEdgeLenRank - minEdgeLenRank);
	for (auto& r : pc.nodesWithDrops) {
		NodeId nodeNumToDrop;
		is >> nodeNumToDrop;
		r = vector<NodeId>(nodeNumToDrop);
		--maxEdgeLenRank;
		for (auto& node : r) {
			is >> node;

			// 缩小半径到最小半径
			if (maxEdgeLenRank >= obj) {
				//cerr << maxEdgeLenRank << endl;
				pc.coverages[node].pop_back();
				pc.coveredNodeNums[node]--;
				if (pc.coveredNodeNums[node] == 0) {
					pc.fixNodes.Nodes[node] = true;
					++pc.fixNodes.Num;
				}
			}
		}
	}

	cerr << maxEdgeLenRank << endl;
}

void saveOutput(ostream& os, Centers& centers) {
	for (auto center = centers.begin(); center != centers.end(); ++center) { os << *center << endl; }
}

void test(istream& inputStream, ostream& outputStream, long long secTimeout, int randSeed) {
	cerr << "load input." << endl;
	PCenter pc;
	loadInput(inputStream, pc);

	vector<Nodes> coverages;
	int nodeNum;
	//coverages = pc.coverages;
	//nodeNum = pc.nodeNum;

	cerr << "solve." << endl;
	steady_clock::time_point start = steady_clock::now();
	steady_clock::time_point endTime = steady_clock::now() + seconds(secTimeout);
	Centers centers(pc.centerNum);
	solvePCenter(centers, pc, [&]() { return duration_cast<milliseconds>(endTime - steady_clock::now()).count(); }, randSeed);
	steady_clock::time_point end = steady_clock::now();

	cerr << "save output." << endl;
	saveOutput(outputStream, centers);

	// DEBUG
	//  // 输出执行时间
	//  auto duration = duration_cast<milliseconds>(end - start);
	//  cerr << "Execution time: " << duration.count() << "ms" << endl;

	//  // 输出解
	//  cerr << "Solution centers: " << endl;
	//  for (NodeId i : centers) {
	//  	cerr << i << ' ';
	//  }
	//  cerr << endl;

	//  // 输出未覆盖的节点
	//  cerr << "Uncovered centers: " << endl;
	//  vector<bool> flag(pc.nodeNum, 0);
	//  for (NodeId i : centers)
	//  	for (NodeId j : pc.coverages[i])
	//  		flag[j] = true;

	//  int tot = 0;
	//  for (NodeId i = 0; i < pc.nodeNum; i++) {
	//  	if (!flag[i]) {
	//  		 cerr << i << endl;
	//  		tot++;
	//  	}
	//  }

	//  cerr << "uncovered center num: " << tot << endl;
}

void test(istream& inputStream, ostream& outputStream, long long secTimeout) {
	return test(inputStream, outputStream, secTimeout, static_cast<int>(time(nullptr) + clock()));
}

void check(istream& ansStream, istream& refStream) {
	vector<int> ans, ref;

	int t;
	while (ansStream >> t)
		ans.push_back(t);
	while (refStream >> t)
		ref.push_back(t);

	for (int i : ref) {
		if (find(ans.begin(), ans.end(), i) == ans.end()) {
			cerr << "Error: " << i << " not found in answer." << endl;
		}
	}
};


int main(int argc, char* argv[]) {
	cerr << "load environment." << endl;
	if (argc > 2) {
		long long secTimeout = atoll(argv[1]);
		int randSeed = atoi(argv[2]);
		test(cin, cout, secTimeout, randSeed);
	}
	else {
		ifstream ifs("instance/pcb3038p200r141.txt");
		ofstream ofs("instance/solution.txt");
		test(ifs, ofs, 30000); // for self-test.

		ifs.close();
		ofs.close();
		//ifstream ref("instance/input1_ref.txt");
		//ifstream ans("instance/solution.txt");
		//check(ans, ref);
	}
	return 0;
}
