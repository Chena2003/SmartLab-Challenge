#include "PCenter.h"

#include <random>
#include <iostream>
#include <unordered_set>
#include <algorithm>
#include <cmath>
#include <chrono>
#include <cstring>

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

		void greedySearch(Centers& solution, solverNodes& solver, function<long long()> restMilliSec) {
			PCenter& nodes = solver.Nodes;
			Flags coveredCenters(nodes.nodeNum, false);
			vector<int> uncoveredNodeCounts(nodes.nodeNum, 0);

			for (NodeId i = 0; i < nodes.nodeNum; ++i) {
				uncoveredNodeCounts[i] = nodes.coveredNodeNums[i];
			}

			NodeId k = 0;
			NodeId* bestNode = new int[nodes.nodeNum];
			int len = 0;
			while (restMilliSec() > RESSERVE && k < nodes.centerNum) {
				// 直接选择固定节点
				if (k < solver.Nodes.fixNodes.Num) {
					solver.inCenter[solution[k]] = true;
					coveredCenters[solution[k]] = true;
					--uncoveredNodeCounts[solution[k]];
					++k;
					continue;
				}

				// 找到最大覆盖节点
				int maxCover = 0;
				for (NodeId i = 0; i < nodes.nodeNum; ++i) {
					if (uncoveredNodeCounts[i] > maxCover) {
						maxCover = uncoveredNodeCounts[i];
						len = 0;
						bestNode[len++] = i;
					}
					else if(maxCover && uncoveredNodeCounts[i] == maxCover) {
						bestNode[len++] = i;
					}
				}

				// 覆盖完全部节点退出
				if (len == 0) break;

				// 随机选择一个节点
				int choice = fastRand(len);
				solution[k] = bestNode[choice];
				solver.inCenter[solution[k]] = true;

				// 更新
				for (NodeId coveredNode : nodes.coverages[solution[k]]) {
					if (!coveredCenters[coveredNode]) {
						coveredCenters[coveredNode] = true;
						for (NodeId neighbor : nodes.coverages[coveredNode]) {
							--uncoveredNodeCounts[neighbor];
						}
					}
				}
				++k;
			}

			delete[]bestNode;
			// 有多余空中心点情况
			while (restMilliSec() > RESSERVE && k < nodes.centerNum) {
				int choice = fastRand(nodes.nodeNum);
				if (!solver.inCenter[choice]) {
					solution[k] = choice;
					solver.inCenter[choice] = true;
					++k;
				}
			}
		}

		// 初始化节点被服务的中心点集合
		void initseriveNodes(Centers& solution, solverNodes& solver) {
			for (NodeId i : solution)
				for (NodeId j : solver.Nodes.coverages[i]) {
					solver.seriveNodes[j].insert(i);
					solver.tseriveNodes[j] = i;
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

		// 初始化SLdata
		vector<vector<uint32_t>> SLdata;
		void initSLdata(solverNodes& solver) {
			for (int i = 0; i < solver.Nodes.nodeNum; ++i) {
				SLdata[i][0] = (uint32_t)pow(i, Gamma[0]);
				SLdata[i][1] = (uint32_t)pow(i, Gamma[1]);
				SLdata[i][2] = (uint32_t)pow(i, Gamma[2]);
			}
		}

		// 哈希函数
		bool SLtabu(Centers& solution, vector<bool>& SL) {
			vector<uint32_t> tot(3);
			tot[0] = tot[1] = tot[2] = 0;
			for (auto i : solution) {
				tot[0] = (SLdata[i][0] + tot[0]) % (uint32_t)MODNUM;
				tot[1] = (SLdata[i][1] + tot[1]) % (uint32_t)MODNUM;
				tot[2] = (SLdata[i][2] + tot[2]) % (uint32_t)MODNUM;
			}
			
			return SL[tot[0]] && SL[tot[1]] && SL[tot[2]];
		}

		// 判断swap(i, j)是否禁忌
		// 1. 判断iter > AL_j
		// 2. 判断solution是否出现过
		bool Tabu(Centers& solution, vector<bool>& SL, vector<int>& AL, int iter, int i, int j) {
			if (iter <= AL[j] || iter <= AL[i])
			//if (iter <= AL[j])
				return true;

			vector<uint32_t> tot(3, 0);
			for (auto i : solution) {
				tot[0] = ((uint32_t)pow(i, Gamma[0]) + tot[0]) % (uint32_t)MODNUM;
				tot[1] = ((uint32_t)pow(i, Gamma[1]) + tot[1]) % (uint32_t)MODNUM;
				tot[2] = ((uint32_t)pow(i, Gamma[2]) + tot[2]) % (uint32_t)MODNUM;
			}

			return SL[tot[0]] && SL[tot[1]] && SL[tot[2]];
		}

		// 寻找最佳交换swap(i, j)
		int Tabu_open, Tabu_close;
		int *delta_;
		vector<pair<int, int>> Pair{ 2000 };
		pair<int, int> findPair(Centers& solution, solverNodes& solver, vector<vector<int>> &coverages) {
			pair<int, int> bestpair = { -1, -1 };
			int obj = INT32_MAX;

			// 随机选择未覆盖节点
			NodeId k = solver.ucenters[fastRand(solver.ucenters.size())];

			int len = 0;
			for (NodeId i : coverages[k]) {
				//if (iter < AL[i])
				if(i == Tabu_open || i == Tabu_close)
					continue;

				memcpy(delta_, solver.delta, solver.Nodes.nodeNum * sizeof(int));

				// 将节点i加入中心点
				for (auto v : coverages[i]) {
					if (solver.seriveNodeNums[v] == 1) {
						//delta_[solver.seriveNodes[v][0]] -= solver.weight[v];
						delta_[solver.tseriveNodes[v]] -= solver.weight[v];
					}
				}

				for (int t = solver.Nodes.fixNodes.Num; t != solver.Nodes.centerNum; ++ t) {
					NodeId j = solution[t]; // 交换i, j

					// 判断swap(i, j)是否为禁忌
					if (j == Tabu_open || j == Tabu_close)
						continue;

					// 更新
					int u = delta_[j] - delta_[i];

					if (u < obj) {
						len = 0;
						Pair[len++] = {i, j};
						obj = u;
					}
					else if (u == obj) {
						Pair[len++] = {i, j};
					}
				}
			}

			if (len == 0)
				return { -1, -1 };

			int choice = fastRand(len);
			return Pair[choice];
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

				solver.tseriveNodes[v] = i;
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
					solver.tseriveNodes[v] = solver.seriveNodes[v][0];
				}
				else if (solver.seriveNodeNums[v] == 0) { // 节点v只有j中心点
					for (NodeId l : solver.Nodes.coverages[v])
						solver.delta[l] += solver.weight[v];
					solver.delta[j] -= solver.weight[v];

					solver.ucenters.insert(v);
				}
			}
		}

		// 判断当前解是否满足下一次半径
		bool check(Centers& output, PCenter& input) {
			Flags flag(input.nodeNum, 0);
			for (NodeId i : output)
				for (NodeId j : input.coverages[i])
					flag[j] = true;

			int tot = 0;
			for (NodeId i = 0; i < input.nodeNum; ++i) 
				if (!flag[i]) 
					++tot;
			return tot == 0;
		}

		int totiter = 0; // 总迭代次数
	public:
		void solve(Centers& output, PCenter& input, function<long long()> restMilliSec, int seed) {
			initRand(seed);
			
			coverAllNodesUnderFixedRadius(output, input, restMilliSec, seed);
			cerr << "Iteration : " << input.U-- << endl;

			// 缩小半径半径
			int n = input.U - input.obj;
			for (int t = 0; (restMilliSec() > RESSERVE) && (t < n); ++ t) {
				reduceRadius(input, input.nodesWithDrops[t]);

				--input.U;
				//if (--input.U > input.obj)
				//	continue;

				//if (check(output, input)) {
				//	cerr << "Iteration : " << input.U-- << " skip. " << endl;
				//	//continue;
				//}
				//else {
				//	if (!coverAllNodesUnderFixedRadius(output, input, restMilliSec, seed)) {
				//		break;
				//	}
				//	cerr << "Iteration : " << input.U-- << endl;
				//}
			}

			if (check(output, input))
				return ;

			cerr << input.U << ' ' << input.obj << endl;
			coverAllNodesUnderFixedRadius(output, input, restMilliSec, seed);
			
			cerr << "Tototal Iters: " << totiter << endl;
		}

		bool coverAllNodesUnderFixedRadius(Centers& output, PCenter& input, function<long long()> restMilliSec, int seed) {
			solverNodes solver(input); // 初始化解结构

			Tabu_open = Tabu_close = -1;
			
			solver.inCenter.assign(solver.Nodes.nodeNum, 0);  // 开始所有点都没有加入中心点
			solver.seriveNodes.resize(solver.Nodes.nodeNum); // 初始化节点被服务的中心点集合
			solver.ucenters.init(solver.Nodes.nodeNum); // 初始化未覆盖节点
			solver.delta = new int[solver.Nodes.nodeNum];
			memset(solver.delta, 0, solver.Nodes.nodeNum * sizeof(int));
			delta_ = new int[solver.Nodes.nodeNum];
			vector<vector<int>> coverages(input.coverages);

			// 统一初始化，减少时间开销
			for (auto i = 0; i < solver.Nodes.nodeNum; ++i) {
				solver.weight.push_back(1);
				//solver.age.push_back(0);
				//solver.delta.push_back(0);
				solver.seriveNodeNums.push_back(0);
				solver.tseriveNodes.push_back(-1);
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
				pair<int, int> t = findPair(solution_, solver, coverages);

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

						totiter += iter;
						break;
					}
				}
				else if (solver.ucenters.size() >= lastUCenterNum ) {
				//else {
					for (int i = 0; i < solver.ucenters.size(); ++ i) {
						NodeId t = solver.ucenters[i];	
						++ solver.weight[t];

						for (NodeId j : solver.Nodes.coverages[t])
							++ solver.delta[j];
					}
				}

				// 更新
				Tabu_open = t.first;
				Tabu_close = t.second;
				solution_ = solution;
				lastUCenterNum = solver.ucenters.size();
			}

			delete[] solver.delta;
			delete[] delta_;
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
