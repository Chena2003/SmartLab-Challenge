#include "PCenter.h"

#include <random>
#include <iostream>
#include <unordered_set>
#include <algorithm>


using namespace std;


namespace szx {

class Solver {
private:
	// random number generator.
	mt19937 pseudoRandNumGen;
	void initRand(int seed) { pseudoRandNumGen = mt19937(seed); }
	int fastRand(int lb, int ub) { return (pseudoRandNumGen() % (ub - lb)) + lb; }
	int fastRand(int ub) { return pseudoRandNumGen() % ub; }
	int rand(int lb, int ub) { return uniform_int_distribution<int>(lb, ub - 1)(pseudoRandNumGen); }
	int rand(int ub) { return uniform_int_distribution<int>(0, ub - 1)(pseudoRandNumGen); }

	using UCenters = unordered_set<int>;
	#define MODNUM 10e8
	vector<bool> inCenters;

	// 贪婪搜索，得到初始解
	void greedySearch(Centers& output, PCenter& input, function<long long()> restMilliSec) {
		// 每次加入覆盖未覆盖节点最多的节点
		PCenter node(input);
		UCenters coveredCenter;

		for(NodeId k = 0; (restMilliSec() > 0) && (k < node.centerNum); k ++) {
			vector<pair<int, int>> t;
			vector<int> cnts(0);
			
			for(NodeId i = 0; (restMilliSec() > 0) && (i < node.nodeNum); i ++) {
				cnts[node.coveredNodeNums[i]] ++;
				t.push_back({node.coveredNodeNums[i], i});
			}

			sort(t.begin(), t.end(),  [](const pair<int, int>& a, const pair<int, int>& b) { return a.first > b.first;});

			int choice = rand(cnts[t[0].first]);
			output[k] = t[choice].second;

			// 更新
			for(auto i : node.coverages[output[k]]) 
				coveredCenter.insert(i);

			for(NodeId i = 0; i < node.nodeNum; i ++) {
				for(auto j : node.coverages[i]) {
					if(coveredCenter.find(j) != coveredCenter.end())
						node.coveredNodeNums[i] --;
				}
			}
		}
	}

	// 初始化delta
	void initDelta(Centers& solution, PCenter& input, Deltas& delta) {
		for(NodeId i = 0; i < input.nodeNum; i ++) {
			Centers t(solution);
			if(inCenters[i]) {
				t.erase(remove(t.begin(), t.end(), i), t.end());
				UCenters ucenter = U(t, input);
				for(auto j : input.coverages[i]) {
					if(ucenter.find(j) == ucenter.end())
						delta[i] += input.weight[j];
				}
			}
			else {
				UCenters ucenter = U(t, input);
				for(auto j : input.coverages[i])
					if(ucenter.find(j) == ucenter.end())
					delta[i] += input.weight[j];
			}
		}
	}	

	// U 操作
	UCenters U(Centers& solution, PCenter& input) {
		UCenters ucenter;
		for(auto i : solution) {
			for(auto j : input.coverages[i]) {
				ucenter.insert(j);
			}
		}

		return move(ucenter);
	}

	void findPair() {


	}


	void makeMove() {


	}


public:
	void solve(Centers& output, PCenter& input, function<long long()> restMilliSec, int seed) {
		initRand(seed);

		// initFixCenterx(); 寻找只覆盖器自己的节点

		coverAllNodesUnderFixedRadius(output, input, restMilliSec, seed);
		for (auto r = input.nodesWithDrops.begin(); (restMilliSec() > 0) && (r != input.nodesWithDrops.end()); ++r) {
			reduceRadius(input, *r);
			if(!coverAllNodesUnderFixedRadius(output, input, restMilliSec, seed)) {
				break;
			}
		}
	}

	bool coverAllNodesUnderFixedRadius(Centers& output, PCenter& input, function<long long()> restMilliSec, int seed) {
		vector<bool> SL(MODNUM, 0);       // 初始化SL
		vector<int> AL(input.nodeNum, 0); // 初始化AL
		vector<int> age(input.nodeNum, 0);// 初始化节点年龄

		Centers solution(input.centerNum);
		inCenters.assign(input.nodeNum, 0);// 初始化节点是否被选中

		greedySearch(solution, input, restMilliSec);// 求解初始解

		Deltas delta(input.nodeNum, 0); // 初始节点delta
		initDelta(solution, input, delta);

		Centers solution_(solution), output(solution);

		bool findSolution = false;
		for(int iter = 1; (restMilliSec() > 0) && U(solution, input).size() == 0; iter ++) {
			findPair();

			makeMove();

			// 比较

			// 更新
		}

		if(findSolution) {
			output = move(output);
			return true;
		}

		return false;

		// TODO: the following code in this function is for illustration only and can be deleted.
		// print some information for debugging.
		cerr << input.nodeNum << '\t' << input.centerNum << endl;
		for (NodeId n = 0; (restMilliSec() > 0) && (n < input.centerNum); ++n) { cerr << n << '\t' << output[n] << endl; }
	}

	void reduceRadius(PCenter& input, Nodes nodesWithDrop) {
		for (auto n = nodesWithDrop.begin(); n != nodesWithDrop.end(); ++n) {
			input.coverages[*n].pop_back();
		}
	}
};

// solver.
void solvePCenter(Centers& output, PCenter& input, function<long long()> restMilliSec, int seed) {
	Solver().solve(output, input, restMilliSec, seed);
}

}
