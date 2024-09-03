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

		const double Gamma[3] = { 1.3, 1.7, 2.3 }; // 哈希参数
		const int RESSERVE = 1;  // 预留1ms退出时间

		// 寻找只覆盖自己的节点
		void initOutput(Centers& output, PCenter& input) {
			for (NodeId i = 0, k = 0; i < input.nodeNum; i++) {
				if (input.fixNodes[i]) {
					output[k++] = i;
				}
			}
		}

		// 贪婪搜索，得到初始解
		void greedySearch(Centers& solution, solverNodes& solver, function<long long()> restMilliSec) {
			// 每次加入覆盖未覆盖节点最多的节点
			PCenter nodes(solver.Nodes);
			Flags coveredCenters(nodes.nodeNum, 0);

			NodeId k;
			for (k = 0; (restMilliSec() > RESSERVE) && (k < nodes.centerNum); k++) {
				if (k < solver.Nodes.fixNum) {
					solver.inCenter[solution[k]] = true;
					coveredCenters[solution[k]] = true;
					nodes.coveredNodeNums[solution[k]]--;
					continue;
				}

				vector<pair<int, int>> t;
				vector<int> cnts(nodes.nodeNum, 0);

				for (NodeId i = 0; i < nodes.nodeNum; i++) {
					if (nodes.coveredNodeNums[i] == 0)
						continue;

					cnts[nodes.coveredNodeNums[i]] ++;
					t.push_back({ nodes.coveredNodeNums[i], i });
				}

				// 覆盖全部节点，直接退出
				if (t.empty())
					break;

				sort(t.begin(), t.end(), [](const pair<int, int>& a, const pair<int, int>& b) { return a.first > b.first; });

				int choice = rand(cnts[t[0].first]);
				solution[k] = t[choice].second;
				solver.inCenter[t[choice].second] = true;

				// 更新
				for (auto i : nodes.coverages[solution[k]]) {
					if (coveredCenters[i] == true)
						continue;

					coveredCenters[i] = true;

					for (NodeId j = 0; j < nodes.nodeNum; j++)
						if (nodes.serives[i][j] == true)
							nodes.coveredNodeNums[j] --;
				}
			}

			// 有多余空中心点情况
			while ((restMilliSec() > RESSERVE) && k < nodes.centerNum) {
				int choice = rand(nodes.nodeNum);

				if (!solver.inCenter[choice]) {
					solution[k] = choice;
					solver.inCenter[choice] = true;

					k++;
				}
			}
		}

		// 初始化节点被服务的中心点集合
		void initseriveNodes(Centers& solution, solverNodes& solver) {
			for (NodeId i : solution)
				for (NodeId j : solver.Nodes.coverages[i])
					solver.seriveNodes[j].insert(i);
		}

		// 初始化未覆盖节点
		void initUCenters(solverNodes& solver) {
			for (NodeId i = 0; i < solver.Nodes.nodeNum; i++) {
				if (solver.seriveNodes[i].empty())
					solver.ucenters.insert(i);
			}
		}

		// 初始化delta
		void initDelta(Centers& solution, solverNodes& solver) {
			PCenter nodes(solver.Nodes);
			for (NodeId i = 0; i < nodes.nodeNum; i++) {
				for (NodeId j : nodes.coverages[i]) {
					if (solver.seriveNodes[j].size() == solver.inCenter[i])
						solver.delta[i] += solver.weight[j];
				}
			}
		}

		// 更新delta
		void updateDelta(solverNodes& solver, NodeId i) {
			for (auto v : solver.Nodes.coverages[i]) {
				if (solver.seriveNodes[v].size() == 1) {
					solver.delta[*solver.seriveNodes[v].begin()] -= solver.weight[v];
					solver.seriveNodes[v].insert(i);
				}

				solver.ucenters.erase(v);
			}
		}

		// 哈希函数
		vector<uint32_t> hash(Centers& solution) {
			vector<uint32_t> tot(3, 0);
			for (auto i : solution) {
				tot[0] = (((uint32_t)pow(i, Gamma[0]) % (uint32_t)MODNUM) + tot[0]) % (uint32_t)MODNUM;
				tot[1] = (((uint32_t)pow(i, Gamma[1]) % (uint32_t)MODNUM) + tot[1]) % (uint32_t)MODNUM;
				tot[2] = (((uint32_t)pow(i, Gamma[2]) % (uint32_t)MODNUM) + tot[2]) % (uint32_t)MODNUM;
			}
			return tot;
		}

		// 判断swap(i, j)是否禁忌
		// 1. 判断iter > AL_j
		// 2. 判断solution是否出现过
		bool Tabu(Centers& solution, vector<bool>& SL, vector<int>& AL, int iter, int i, int j) {
			if (iter <= AL[j])
				return true;

			vector<uint32_t> t = hash(solution);
			return SL[t[0]] & SL[t[1]] & SL[t[2]];
		}

		// 寻找最佳交换swap(i, j)
		pair<int, int> findPair(Centers& solution, solverNodes& solver, vector<bool>& SL, vector<int>& AL, int iter) {
			pair<int, int> bestpair = { -1, -1 };
			uint32_t obj = INTMAX;

			PCenter nodes = solver.Nodes;
			UCenters ucenters(solver.ucenters);

			// 随机选择未覆盖节点
			auto t = ucenters.begin();
			advance(t, rand(ucenters.size()));
			NodeId k = *t;

			// 计算f(x)
			uint32_t fx = 0;
			for (auto i : ucenters)
				fx += solver.weight[i];

			solverNodes solver_(solver);
			for (NodeId i = 0; i < nodes.nodeNum; i++) {
				if (nodes.serives[k][i] == false)
					continue;

				// 将节点i加入中心点
				updateDelta(solver_, i);

				for (auto t = solution.begin(); t != solution.end(); t++) {
					if (solver.Nodes.fixNodes[*t]) // 固定节点不能换出
						continue;

					NodeId j = *t; // 交换i, j
					*t = i;

					// 判断swap(i, j)是否为禁忌
					if (Tabu(solution, SL, AL, iter, i, j)) {
						*t = j;
						continue;
					}

					// 更新
					uint32_t u = fx - solver_.delta[i] + solver_.delta[j];
					if (u < obj) {
						obj = u;
						bestpair = { i, j };
					}
					else if (u == obj && solver_.age[j] < solver_.age[bestpair.second]) {
						bestpair = { i, j };
					}

					*t = j; // 恢复
				}

				solver_ = solver;
			}

			return bestpair;
		}


		void makeMove(Centers& solution, solverNodes& solver, int i, int j) {
			// 加入i
			for (NodeId v : solver.Nodes.coverages[i]) {
				if (solver.seriveNodes[v].size() == 1) { // 节点v只有一个中心点
					solver.delta[*solver.seriveNodes[v].begin()] -= solver.weight[v];
				}
				else if (solver.seriveNodes[v].size() == 0) { // 节点v没有中心点
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

				if (solver.seriveNodes[v].size() == 1) { // 节点v有两个中心点
					solver.delta[*solver.seriveNodes[v].begin()] += solver.weight[v];
				}
				else if (solver.seriveNodes[v].size() == 0) { // 节点v只有j中心点
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
			initOutput(output, input); // 初始化固定节点

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
			solver.seriveNodes.resize(solver.Nodes.nodeNum); // 初始化服务节点

			vector<bool> SL(MODNUM, 0);       // 初始化SL
			vector<int> AL(solver.Nodes.nodeNum, 0); // 初始化AL

			Centers solution(output);
			greedySearch(solution, solver, restMilliSec); // 求解初始解

			initseriveNodes(solution, solver); // 初始化节点被服务中心点数
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
				// 寻找最佳交换
				pair<int, int> t = findPair(solution_, solver, SL, AL, iter);

				if (t.first == -1 && t.second == -1)
					break;

				// 交换
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
				else {
					for (NodeId i : solver.ucenters) {
						solver.weight[i] ++;

						for (NodeId j : solver.Nodes.coverages[i])
							solver.delta[j] ++;
					}
				}

				// 更新
				vector<uint32_t> hashNum = hash(solution);
				for (int k = 0; k < 3; k++)
					SL[hashNum[k]] = true;
				AL[t.first] = iter + 2;
				solver.age[t.second] = iter;
				solution_ = solution;
			}

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
