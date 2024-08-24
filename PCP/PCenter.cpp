#include "PCenter.h"

#include <random>
#include <iostream>
#include <unordered_set>
#include <algorithm>
#include <cmath>

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

		const double Gamma[3] = { 1.3, 1.7, 2.3 };
		const int RESSERVE = 1;

		// 贪婪搜索，得到初始解
		void greedySearch(Centers& solution, solverNodes& solver, function<long long()> restMilliSec) {
			// 每次加入覆盖未覆盖节点最多的节点
			PCenter nodes(solver.Nodes);
			UCenters coveredCenters;

			for (NodeId k = 0; (restMilliSec() > RESSERVE) && (k < nodes.centerNum); k ++) {
				vector<pair<int, int>> t;
				vector<int> cnts(nodes.nodeNum, 0);
				UCenters coveredCenter;

				for (NodeId i = 0; i < nodes.nodeNum; i ++) {
					cnts[nodes.coveredNodeNums[i]] ++;
					t.push_back({nodes.coveredNodeNums[i], i });
				}

				sort(t.begin(), t.end(), [](const pair<int, int>& a, const pair<int, int>& b) { return a.first > b.first; });

				int choice = rand(cnts[t[0].first]);
				solution[k] = t[choice].second;
				solver.inCenter[t[choice].second] = true;

				// 更新
				for (auto i : nodes.coverages[solution[k]]) {
					if (coveredCenters.find(i) == coveredCenters.end())
						coveredCenters.insert(i);

					coveredCenters.insert(i);
				}

				for (NodeId i = 0; i < nodes.nodeNum; i++) {
					for (auto j : nodes.coverages[i]) {
						if (coveredCenter.find(j) != coveredCenter.end())
							nodes.coveredNodeNums[i] --;
					}
				}
			}
		}

		// 初始化delta
		void initDelta(Centers& solution, solverNodes& solver) {
			PCenter nodes(solver.Nodes);
			for (NodeId i = 0; i < nodes.nodeNum; i++) {
				if (solver.inCenter[i]) {
					for (NodeId j : nodes.coverages[i]) {
						if (solver.seriveNodes[j].size() == 1)
							solver.delta[i] += solver.weight[j];
					}
				}
				else {
					for (NodeId j : nodes.coverages[i])
						if (solver.seriveNodes[j].size() == 0)
							solver.delta[i] += solver.weight[j];
				}
			}
		}

		void initseriveNodeNums(Centers& solution, solverNodes& solver) {
			for (NodeId i : solution) {
				for (NodeId j : solver.Nodes.coverages[i]) {
					solver.seriveNodes[j].insert(i);
				}
			}
		}

		void initUCenters(solverNodes& solver) {
			for (NodeId i = 0; i < solver.Nodes.nodeNum; i++) {
				if (solver.seriveNodes[i].empty())
					solver.ucenters.insert(i);
			}
		}

		// U 操作
		UCenters U(solverNodes& solver) {
			UCenters ucenter;
			for (NodeId i = 0; i < solver.Nodes.nodeNum; i++) {
				if (solver.seriveNodes[i].empty())
					ucenter.insert(i);
			}

			return move(ucenter);
		}

		// 更新delta
		void updateDelta(Centers& solution, solverNodes& solver, NodeId i) {
			for (auto v : solver.Nodes.coverages[i]) {
				if (solver.seriveNodes[v].size() == 1) {
					solver.delta[*solver.seriveNodes[v].begin()] -= solver.weight[v];
					solver.seriveNodes[v].insert(i);
				}
			}
		}

		// 哈希函数
		uint32_t hash(Centers& solution, int k) {
			uint32_t tot = 0;
			for (auto i : solution)
				tot = (((uint32_t)pow(i, Gamma[k]) % (uint32_t)MODNUM) + tot) % (uint32_t)MODNUM;

			return tot;
		}

		// 判断swap(i, j)是否禁忌
		// iter > AL_j
		// 判断solution是否出现过
		bool Tabu(Centers& solution, vector<bool>& SL, vector<int>& AL, int iter, int i, int j) {
			if (iter > AL[j])
				return false;

			for (int k = 0; k < 3; k++) {
				if (!SL[hash(solution, k)]) {
					return false;
				}
			}

			return true;
		}

		pair<int, int> findPair(Centers& solution, solverNodes& solver, vector<bool>& SL, vector<int>& AL, int iter) {
			pair<int, int> bestpair = { -1, -1 };
			uint32_t obj = INTMAX;

			solverNodes solver_(solver);
			PCenter nodes = solver.Nodes;
			UCenters ucenters(solver.ucenters);
			auto i = ucenters.begin();
			advance(i, rand(ucenters.size()));
			NodeId k = *i;

			Deltas delta_(solver.delta);
			for (NodeId i = 0; i < nodes.nodeNum; i++) {
				if (nodes.serives[k][i] == false)
					continue;

				updateDelta(solution, solver_, i);

				for (auto j : solution) {
					if (Tabu(solution, SL, AL, iter, i, j))
						continue;

					uint32_t u = obj - solver_.delta[i] + solver_.delta[j];

					if (u < obj) {
						obj = u;
						bestpair = { i, j };
					}
					else if (u == obj && solver_.age[j] < solver_.age[bestpair.second]) {
						bestpair = { i, j };
					}
				}

				solver_.delta = delta_;
			}

			return bestpair;
		}


		void makeMove(Centers& solution, solverNodes& solver, int i, int j) {
			// 加入i
			for (NodeId v : solver.Nodes.coverages[i]) {
				if (solver.seriveNodes[v].size() == 1) {
					solver.delta[*solver.seriveNodes[v].begin()] -= solver.weight[v];
				}
				else if (solver.seriveNodes[v].size() == 0) {
					for (NodeId l = 0; l < solver.Nodes.nodeNum; l++)
						if (l != i && solver.Nodes.serives[v][l] == true)
							solver.delta[l] -= solver.weight[v];

					solver.ucenters.erase(v);
				}

				solver.seriveNodes[v].insert(i);
			}

			// swap(i, j)
			for (auto t = solution.begin(); t != solution.end(); t++) {
				if (*t == j) {
					*t = i;
					solver.inCenter[j] = false;
					solver.inCenter[i] = true;
					break;
				}
			}

			// 关闭j
			for (NodeId v : solver.Nodes.coverages[j]) {
				solver.seriveNodes[v].erase(j);

				if (solver.seriveNodes[v].size() == 2) {
					solver.delta[*solver.seriveNodes[v].begin()] += solver.weight[v];
				}
				else if (solver.seriveNodes[v].size() == 1) {
					for (NodeId l = 0; l < solver.Nodes.nodeNum; l++)
						if (l != j && solver.Nodes.serives[v][l] == true)
							solver.delta[l] += solver.weight[v];

					solver.ucenters.insert(v);
				}
			}
		}

	public:
		void solve(Centers& output, PCenter& input, function<long long()> restMilliSec, int seed) {
			initRand(seed);

			coverAllNodesUnderFixedRadius(output, input, restMilliSec, seed);
			for (auto r = input.nodesWithDrops.begin(); (restMilliSec() > RESSERVE) && (r != input.nodesWithDrops.end()); ++r) {
				reduceRadius(input, *r);
				if (!coverAllNodesUnderFixedRadius(output, input, restMilliSec, seed)) {
					break;
				}
			}
		}

		bool coverAllNodesUnderFixedRadius(Centers& output, PCenter& input, function<long long()> restMilliSec, int seed) {
			solverNodes solver; // 初始化解结构
			solver.Nodes = input;
			solver.weight.assign(solver.Nodes.nodeNum, 1); // 初始化节点权重为1
			solver.inCenter.assign(solver.Nodes.nodeNum, 0);  // 开始所有点都没有加入中心点
			solver.age.assign(solver.Nodes.nodeNum, 0); // 初始化节点年龄为0
			solver.delta.assign(solver.Nodes.nodeNum, 0); // 初始化节点delta为0
			solver.seriveNodes.resize(solver.Nodes.nodeNum);

			// initFixCenterx(); 寻找只覆盖器自己的节点

			vector<bool> SL(MODNUM, 0);       // 初始化SL
			vector<int> AL(solver.Nodes.nodeNum, 0); // 初始化AL

			Centers solution(solver.Nodes.centerNum);
			greedySearch(solution, solver, restMilliSec);// 求解初始解

			initseriveNodeNums(solution, solver); // 初始化节点被服务中心点数
			initUCenters(solver); // 初始化ucenter

			if (solver.ucenters.size() == 0) {
				output = solution;
				return true;
			}

			initDelta(solution, solver); // 初始化节点delta

			bool findSolution = false;
			Centers solution_(solution), bestsolution(solution);
			UCenters bestUCenters(solver.ucenters);
			for (int iter = 1; restMilliSec() > RESSERVE; iter++) {
				pair<int, int> t = findPair(solution_, solver, SL, AL, iter);

				if (t.first == -1 && t.second == -1)
					break;

				makeMove(solution, solver, t.first, t.second);

				// 比较
				if (solver.ucenters.size() < bestUCenters.size()) {
					bestsolution = solution;
					bestUCenters = solver.ucenters;

					if (bestUCenters.size() == 0) {
						output = bestsolution;
						findSolution = true;
						break;
					}
				}
				else if (solver.ucenters.size() > bestUCenters.size()) {
					for (NodeId i : solver.ucenters) {
						solver.weight[i] ++;

						for (NodeId j : solver.Nodes.coverages[i])
							solver.delta[j] ++;
					}
				}

				// 更新
				for (int k = 0; k < 3; k++)
					SL[hash(solution, k)] = true;
				AL[t.second] = iter + 2;
				solver.age[t.second] = iter;
			}

			// TODO: the following code in this function is for illustration only and can be deleted.
			// print some information for debugging.
//			cerr << input.nodeNum << '\t' << input.centerNum << endl;
//			for (NodeId n = 0; (restMilliSec() > 0) && (n < input.centerNum); ++n) { cerr << n << '\t' << output[n] << endl; }

			return findSolution;
		}

		void reduceRadius(PCenter& input, Nodes nodesWithDrop) {
			for (auto n = nodesWithDrop.begin(); n != nodesWithDrop.end(); ++n) {
				int j = input.coverages[*n].back();
				input.coverages[*n].pop_back();
				input.serives[j][*n] = false;
				input.coveredNodeNums[*n] --;
			}
		}
	};

	// solver.
	void solvePCenter(Centers& output, PCenter& input, function<long long()> restMilliSec, int seed) {
		Solver().solve(output, input, restMilliSec, seed);
	}

}
