#include "PCenter.h"

#include <random>
#include <iostream>
#include <unordered_set>
#include <algorithm>
#include <cmath>
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

		const double Gamma[3] = { 1.3, 1.7, 2.3 }; // 哈希参数
		const int RESSERVE = 1;  // 预留1ms退出时间

		// 寻找只覆盖自己的节点
		void initOutput(Centers& output, PCenter& input) {
			for (NodeId i = 0, k = 0; i < input.nodeNum; i++) {
				if (input.fixNodes.Nodes[i]) {
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
				if (k < solver.Nodes.fixNodes.Num) {
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

				int choice = fastRand(cnts[t[0].first]);
				solution[k] = t[choice].second;
				solver.inCenter[t[choice].second] = true;

				// 更新
				for (auto i : nodes.coverages[solution[k]]) {
					if (coveredCenters[i] == true)
						continue;

					coveredCenters[i] = true;

					for (NodeId j : nodes.coverages[i])
						nodes.coveredNodeNums[j] --;
				}
			}

			// 有多余空中心点情况
			while ((restMilliSec() > RESSERVE) && k < nodes.centerNum) {
				int choice = fastRand(nodes.nodeNum);

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
				for (NodeId j : solver.Nodes.coverages[i]) {
					solver.seriveNodes[j].insert(i);
					++ solver.seriveNodeNums[j];
				}
		}

		// 初始化未覆盖节点
		void initUCenters(solverNodes& solver) {
			for (NodeId i = 0; i < solver.Nodes.nodeNum; i++) {
				if (solver.seriveNodeNums[i] == 0) {
					solver.ucenters.insert(i);
				}
			}
		}

		// 初始化delta
		void initDelta(Centers& solution, solverNodes& solver) {
			PCenter nodes(solver.Nodes);
			for (NodeId i = 0; i < nodes.nodeNum; i++) {
				for (NodeId j : nodes.coverages[i]) {
					if (solver.seriveNodeNums[j] == solver.inCenter[i])
						solver.delta[i] += solver.weight[j];
				}
			}
		}

		// 更新delta
		void updateDelta(solverNodes& solver, NodeId i) {
			for (auto v : solver.Nodes.coverages[i]) {
				if (solver.seriveNodeNums[v] == 1) {
					solver.delta[solver.seriveNodes[v][0]] -= solver.weight[v];
				}
			}
		}

		// 哈希函数
		vector<uint32_t> hash(Centers& solution) {
			vector<uint32_t> tot(3);
			tot[0] = tot[1] = tot[2] = 0;
			for (auto i : solution) {
				//tot[0] = (((uint32_t)pow(i, Gamma[0]) % (uint32_t)MODNUM) + tot[0]) % (uint32_t)MODNUM;
				//tot[1] = (((uint32_t)pow(i, Gamma[1]) % (uint32_t)MODNUM) + tot[1]) % (uint32_t)MODNUM;
				//tot[2] = (((uint32_t)pow(i, Gamma[2]) % (uint32_t)MODNUM) + tot[2]) % (uint32_t)MODNUM;
			
				tot[0] = ((uint32_t)pow(i, Gamma[0]) + tot[0]) % (uint32_t)MODNUM;
				tot[1] = ((uint32_t)pow(i, Gamma[1]) + tot[1]) % (uint32_t)MODNUM;
				tot[2] = ((uint32_t)pow(i, Gamma[2]) + tot[2]) % (uint32_t)MODNUM;
			}
			return tot;
		}

		// 判断swap(i, j)是否禁忌
		// 1. 判断iter > AL_j
		// 2. 判断solution是否出现过
		bool Tabu(Centers& solution, vector<bool>& SL, vector<int>& AL, int iter, int i, int j) {
			if (iter <= AL[j] || iter <= AL[i])
			//if (iter <= AL[j])
				return true;

			//vector<uint32_t> t = hash(solution);
			vector<uint32_t> tot(3);
			tot[0] = tot[1] = tot[2] = 0;
			for (auto i : solution) {
				tot[0] = ((uint32_t)pow(i, Gamma[0]) + tot[0]) % (uint32_t)MODNUM;
				tot[1] = ((uint32_t)pow(i, Gamma[1]) + tot[1]) % (uint32_t)MODNUM;
				tot[2] = ((uint32_t)pow(i, Gamma[2]) + tot[2]) % (uint32_t)MODNUM;
			}

			return SL[tot[0]] && SL[tot[1]] && SL[tot[2]];
		}

		// 寻找最佳交换swap(i, j)
		pair<int, int> findPair(Centers& solution, solverNodes& solver, vector<int>& AL, int iter) {
			pair<int, int> bestpair = { -1, -1 };
			uint32_t obj = INTMAX;

			// 随机选择未覆盖节点
			NodeId k = solver.ucenters[fastRand(solver.ucenters.size())];

			// 计算f(x)
			uint32_t fx = 100000;
			// 计算f(x)
			//for (auto i : solver.ucenters)
			//	fx += solver.weight[i];

			for (NodeId i : solver.Nodes.coverages[k]) {
				if (iter <= AL[i])
					continue;

				//vector<int> delta_(solver.delta);
				// 将节点i加入中心点
				for (auto v : solver.Nodes.coverages[i]) {
					if (solver.seriveNodeNums[v] == 1) {
						solver.delta[solver.seriveNodes[v][0]] -= solver.weight[v];
					}
				}

				for (int t = solver.Nodes.fixNodes.Num; t != solver.Nodes.centerNum; ++ t) {

					NodeId j = solution[t]; // 交换i, j
					solution[t] = i;

					// 判断swap(i, j)是否为禁忌
					if (iter <= AL[j]) {
						solution[t] = j;
						continue;
					}

					// 更新
					uint32_t u = fx - solver.delta[i] + solver.delta[j];
					if (u < obj) {
						obj = u;
						bestpair = { i, j };
					}
					else if (u == obj && solver.age[j] < solver.age[bestpair.second]) {
						bestpair = { i, j };
					}

					solution[t] = j; // 恢复
				}

				/// 恢复delta
				for (auto v : solver.Nodes.coverages[i]) {
					if (solver.seriveNodeNums[v] == 1) {
						solver.delta[solver.seriveNodes[v][0]] += solver.weight[v];
					}
				}
			}

			return bestpair;
		}


		void makeMove(Centers& solution, solverNodes& solver, int i, int j) {
			// 加入i
			for (NodeId v : solver.Nodes.coverages[i]) {
				if (solver.seriveNodeNums[v] == 1) { // 节点v只有一个中心点
					solver.delta[solver.seriveNodes[v][0]] -= solver.weight[v];
				}
				else if (solver.seriveNodeNums[v] == 0) { // 节点v没有中心点
					for (NodeId l : solver.Nodes.coverages[v])
						solver.delta[l] -= solver.weight[v];
					solver.delta[i] += solver.weight[v];

					solver.ucenters.erase(v);
				}

				solver.seriveNodes[v].insert(i);
				++ solver.seriveNodeNums[v];
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
				-- solver.seriveNodeNums[v];

				if (solver.seriveNodeNums[v] == 1) { // 节点v有两个中心点
					solver.delta[solver.seriveNodes[v][0]] += solver.weight[v];
				}
				else if (solver.seriveNodeNums[v] == 0) { // 节点v只有j中心点
					for (NodeId l : solver.Nodes.coverages[v])
						solver.delta[l] += solver.weight[v];
					solver.delta[j] -= solver.weight[v];

					solver.ucenters.insert(v);
				}
			}
		}

	public:
		void solve(Centers& output, PCenter& input, function<long long()> restMilliSec, int seed) {
			initRand(seed);
			
			uint64_t start = restMilliSec();
			coverAllNodesUnderFixedRadius(output, input, restMilliSec, seed);
			uint64_t end = restMilliSec();
			uint64_t duration = end - start; // 统计每次求解时间

			int n = input.nodesWithDrops.size();
			for (int t = 0; (restMilliSec() > RESSERVE) && (t < n); ++ t) {
				reduceRadius(input, input.nodesWithDrops[t]);
				
				start = restMilliSec();
				steady_clock::time_point startTime = steady_clock::now();
				steady_clock::time_point endTime = steady_clock::now() + milliseconds(5 * duration); // 单次求解时间不能超过上次求解时间的5倍
				bool flag = coverAllNodesUnderFixedRadius(output, input, [&]() { return duration_cast<milliseconds>(endTime - steady_clock::now()).count(); }, seed);
				end = restMilliSec();
				duration = end - start;
				
				if (!flag) {
					break;
				}
			}
		}

		bool coverAllNodesUnderFixedRadius(Centers& output, PCenter& input, function<long long()> restMilliSec, int seed) {
			solverNodes solver; // 初始化解结构
			solver.Nodes = input;

			vector<int> AL(solver.Nodes.nodeNum, 0); // 初始化AL
			solver.inCenter.assign(solver.Nodes.nodeNum, 0);  // 开始所有点都没有加入中心点
			solver.seriveNodes.resize(solver.Nodes.nodeNum); // 初始化节点被服务的中心点集合
			solver.ucenters.init(solver.Nodes.nodeNum); // 初始化未覆盖节点

			// 统一初始化，减少时间开销
			for (auto i = 0; i < solver.Nodes.nodeNum; ++i) {
				solver.weight.push_back(1);
				solver.age.push_back(0);
				solver.delta.push_back(0);
				solver.seriveNodeNums.push_back(0);
				solver.seriveNodes[i].init(solver.Nodes.nodeNum);
			}								 

			Centers solution(input.centerNum);
			initOutput(solution, input); // 初始化固定节点
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
			uint32_t bestUCenterNum = solver.ucenters.size();
			uint32_t lastUCenterNum = bestUCenterNum;
			for (int iter = 1; restMilliSec() > RESSERVE; iter++) {
				// 寻找最佳交换
				pair<int, int> t = findPair(solution_, solver, AL, iter);

				if (t.first == -1 && t.second == -1)
					break;

				// 交换
				makeMove(solution, solver, t.first, t.second);

				// 比较
				if (solver.ucenters.size() < bestUCenterNum) {
					bestsolution = solution;
					bestUCenterNum = solver.ucenters.size();

					if (bestUCenterNum == 0) {
						output = bestsolution;
						findSolution = true;
						break;
					}
				}
				else if (solver.ucenters.size() >= lastUCenterNum ) {
				//else {
					for (int i = 0; i < solver.ucenters.size(); ++ i) {
						NodeId t = solver.ucenters[i];	
						solver.weight[t] ++;

						for (NodeId j : solver.Nodes.coverages[t])
							solver.delta[j] ++;
					}
				}

				// 更新
				AL[t.first] = iter + 2;
				AL[t.second] = iter + 2;
				solver.age[t.second] = iter;
				solution_ = solution;
				lastUCenterNum = solver.ucenters.size();
			}

			return findSolution;
		}

		void reduceRadius(PCenter& input, Nodes nodesWithDrop) {
			for(auto n : nodesWithDrop) {
				input.coverages[n].pop_back();
				input.coveredNodeNums[n] --;

				if (input.coveredNodeNums[n] == 0) {
					input.fixNodes.Nodes[n] = true;
					++ input.fixNodes.Num;
				}
			}
		}
	};

	// solver.
	void solvePCenter(Centers& output, PCenter& input, function<long long()> restMilliSec, int seed) {
		Solver().solve(output, input, restMilliSec, seed);
	}

}
