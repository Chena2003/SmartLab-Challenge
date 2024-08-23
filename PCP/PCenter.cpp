#include "PCenter.h"

#include <random>
#include <iostream>
#include <unordered_set>
#include <algorithm>
#include <cmath>

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

	int Gamma[3] = {1.3, 1.7, 2.3};

	// 贪婪搜索，得到初始解
	void greedySearch(Centers& solution, PCenter& input, function<long long()> restMilliSec) {
		// 每次加入覆盖未覆盖节点最多的节点
		PCenter nodes(input);
		UCenters coveredCenter;

		for(NodeId k = 0; (restMilliSec() > 0) && (k < nodes.centerNum); k ++) {
			vector<pair<int, int>> t;
			vector<int> cnts(0);
			
			for(NodeId i = 0; (restMilliSec() > 0) && (i < nodes.nodeNum); i ++) {
				cnts[nodes.coveredNodeNums[i]] ++;
				t.push_back({nodes.coveredNodeNums[i], i});
			}

			sort(t.begin(), t.end(),  [](const pair<int, int>& a, const pair<int, int>& b) { return a.first > b.first;});

			int choice = rand(cnts[t[0].first]);
			solution[k] = t[choice].second;

			// 更新
			for(auto i : nodes.coverages[solution[k]]) 
				coveredCenter.insert(i);

			for(NodeId i = 0; i < nodes.nodeNum; i ++) {
				for(auto j : nodes.coverages[i]) {
					if(coveredCenter.find(j) != coveredCenter.end())
						nodes.coveredNodeNums[i] --;
				}
			}
		}
	}

	// 初始化delta
	void initDelta(Centers& solution, solverNodes& solver) {
		PCenter nodes(solver.Nodes);
		for(NodeId i = 0; i < nodes.nodeNum; i ++) {
			Centers t(solution);
			if(solver.inCenter[i]) {
				t.erase(remove(t.begin(), t.end(), i), t.end());
				UCenters ucenter = U(t, nodes);
				for(auto j : nodes.coverages[i]) {
					if(ucenter.find(j) == ucenter.end())
						solver.delta[i] += solver.weight[j];
				}
			}
			else {
				UCenters ucenter = U(t, nodes);
				for(auto j : nodes.coverages[i])
					if(ucenter.find(j) == ucenter.end())
					solver.delta[i] += solver.weight[j];
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

	// 更新delta
	void updateDelta(Centers& solution, solverNodes& solver, NodeId i) {
		for(auto v : solver.Nodes.coverages[i]) {
			for(auto x : solution) {
				if(solver.Nodes.serives[v][x] == true) {
					solver.delta[x] -= solver.weight[v];
				}
			}
		}
	}

	uint32_t hash(Centers& solution, int k) {
		
		uint32_t tot = 0;
		for(auto i : solution) 
			tot = (((uint32_t)pow(i, Gamma[k]) % (uint32_t)MODNUM) + tot) % (uint32_t)MODNUM;

		return tot;
	}

	// 判断swap(i, j)是否禁忌
	// iter > AL_j
	// 判断solution是否出现过
	bool Tabu(Centers& solution, vector<bool>& SL, vector<int>& AL, int iter, int i, int j) {
		if(iter > AL[j])
			return false;

		for(int k = 0; k < 3; k ++) {
			if(!SL[hash(solution, k)]) {
				return false;
			}		
		}

		return true;
	}

	pair<int, int> findPair(Centers& solution, solverNodes& solver, vector<bool>& SL, vector<int>& AL, int iter) {
		pair<int, int> bestpair;
		int obj = INTMAX;

		PCenter nodes = solver.Nodes;
		UCenters ucenters = U(solution, nodes);
		NodeId k = rand(ucenters.size());
		
		Deltas delta_(solver.delta);
		for(NodeId i = 0; i < nodes.nodeNum; i ++) {
			if(nodes.serives[k][i] == false)
				continue; 
			
			updateDelta(solution, solver, i);

			for(auto j : solution) {
				if(Tabu(solution, SL, AL, iter, i, j))
					continue;


				
			}

		}
	}


	void makeMove(Centers& solution, int i, int j) {


	}


public:
	void solve(Centers& output, PCenter& input, function<long long()> restMilliSec, int seed) {
		initRand(seed);

		coverAllNodesUnderFixedRadius(output, input, restMilliSec, seed);
		for (auto r = input.nodesWithDrops.begin(); (restMilliSec() > 0) && (r != input.nodesWithDrops.end()); ++r) {
			reduceRadius(input, *r);
			if(!coverAllNodesUnderFixedRadius(output, input, restMilliSec, seed)) {
				break;
			}
		}
	}

	bool coverAllNodesUnderFixedRadius(Centers& output, PCenter& input, function<long long()> restMilliSec, int seed) {
		solverNodes solver; // 初始化解结构
		solver.Nodes = input;
		solver.weight.assign(solver.Nodes.nodeNum, 1); // 初始化节点权重为1
		solver.inCenter.assign(solver.Nodes.nodeNum, 0);  // 开始所有点都没有加入中心点中
		solver.age.assign(solver.Nodes.nodeNum, 0); // 初始化节点年龄为0
		solver.delta.assign(solver.Nodes.nodeNum, 0); // 初始化节点delta为0
 		
		// initFixCenterx(); 寻找只覆盖器自己的节点
		
		vector<bool> SL(MODNUM, 0);       // 初始化SL
		vector<int> AL(solver.Nodes.nodeNum, 0); // 初始化AL

		Centers solution(solver.Nodes.centerNum);
		greedySearch(solution, solver.Nodes, restMilliSec);// 求解初始解

		initDelta(solution, solver); // 初始化节点delta

		bool findSolution = false;
		Centers solution_(solution), output(solution);
		for(int iter = 1; (restMilliSec() > 0) && U(solution, solver.Nodes).size() == 0; iter ++) {
			pair<int, int> t = findPair(solution_, solver, SL, AL, iter);

			makeMove(solution, t.first, t.second);

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
			int j = input.coverages[*n].back();
			input.coverages[*n].pop_back();
			input.serives[j][*n] = false;
 		}
	}
};

// solver.
void solvePCenter(Centers& output, PCenter& input, function<long long()> restMilliSec, int seed) {
	Solver().solve(output, input, restMilliSec, seed);
}

}
