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
		const int INTSIZE = sizeof(int);

		// 全局变量
		NodeId nodeNum, centerNum, fixNum;
		int* coveredNodeNums;
		int** coverages;
		int* weight;
		bool* inCenter;
		int* delta;
		doublevector ucenters;
		doublevector* seriveNodes;
		int* seriveNodeNums;
		int* tseriveNodes;


		// 寻找只覆盖自己的节点
		void initOutput(Centers& output, PCenter& input) {
			for (NodeId i = 0, k = 0; i < input.nodeNum; i++) {
				if (input.fixNodes.Nodes[i]) {
					output[k++] = i;
				}
			}
		}

		void greedySearch(Centers& solution, function<long long()> restMilliSec) {
			Flags coveredCenters(nodeNum, false);
			vector<int> uncoveredNodeCounts(nodeNum, 0);

			for (NodeId i = 0; i < nodeNum; ++i) {
				uncoveredNodeCounts[i] = coveredNodeNums[i];
			}

			NodeId k = 0;
			NodeId* bestNode = new int[nodeNum];
			int len = 0;
			while (restMilliSec() > RESSERVE && k < centerNum) {
				// 直接选择固定节点
				if (k < fixNum) {
					inCenter[solution[k]] = true;
					coveredCenters[solution[k]] = true;
					--uncoveredNodeCounts[solution[k]];
					++k;
					continue;
				}

				// 找到最大覆盖节点
				int maxCover = 0;
				for (NodeId i = 0; i < nodeNum; ++i) {
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
				inCenter[solution[k]] = true;

				// 更新
				for(int i = 0; i < coveredNodeNums[solution[k]]; ++i) {
					NodeId coveredNode = coverages[solution[k]][i];
					if (!coveredCenters[coveredNode]) {
						coveredCenters[coveredNode] = true;
						for (int j = 0; j < coveredNodeNums[coveredNode]; ++j) {
							NodeId neighbor = coverages[coveredNode][j];
							--uncoveredNodeCounts[neighbor];
						}
					}
				}
				++k;
			}

			delete[]bestNode;
			// 有多余空中心点情况
			while (restMilliSec() > RESSERVE && k < centerNum) {
				int choice = fastRand(nodeNum);
				if (!inCenter[choice]) {
					solution[k] = choice;
					inCenter[choice] = true;
					++k;
				}
			}
		}

		// 初始化节点被服务的中心点集合
		void initseriveNodes(Centers& solution) {
			for (NodeId i : solution)
				for (int t = 0; t < coveredNodeNums[i]; ++t) {
					NodeId j = coverages[i][t];
					seriveNodes[j].insert(i);
					tseriveNodes[j] = i;
					++seriveNodeNums[j];
				}
		}

		// 初始化未覆盖节点
		void initUCenters() {
			for (NodeId i = 0; i < nodeNum; ++i) {
				if (seriveNodeNums[i] == 0) {
					ucenters.insert(i);
				}
			}
		}

		// 初始化delta
		void initDelta(Centers& solution) {
			//PCenter nodes(Nodes);
			for (NodeId i = 0; i < nodeNum; ++i) {
				for (int t = 0; t < coveredNodeNums[i]; ++t) {
					NodeId j = coverages[i][t];
					if (seriveNodeNums[j] == inCenter[i])
						delta[i] += weight[j];
				}
			}
		}

		vector<vector<uint32_t>> SLdata;
		void initSLdata(solverNodes& solver) {
			for (int i = 0; i < nodeNum; ++i) {
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
		int* delta_;
		vector<pair<int, int>> Pair{ 2000 };
		pair<int, int> findPair(Centers& solution) {
			pair<int, int> bestpair = { -1, -1 };
			int obj = INT32_MAX;

			// 随机选择未覆盖节点
			NodeId k = ucenters[fastRand(ucenters.size())];

			int len = 0;
			for(int ii = 0; ii < coveredNodeNums[k]; ++ii) {
				NodeId i = coverages[k][ii];
			//for (NodeId i : coverages[k]) {
				//if (iter < AL[i])
				if(i == Tabu_open || i == Tabu_close)
					continue;

				memcpy(delta_, delta, nodeNum * sizeof(int));

				// 将节点i加入中心点
				for(int vv = 0; vv < coveredNodeNums[i]; ++vv) {
					NodeId v = coverages[i][vv];
					if (seriveNodeNums[v] == 1) {
						delta_[tseriveNodes[v]] -= weight[v];
					}
				}

				for (int t = fixNum; t != centerNum; ++ t) {
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


		void makeMove(Centers& solution, int i, int j) {
			// 加入i
			for(int vv = 0; vv < coveredNodeNums[i]; ++vv) {
				NodeId v = coverages[i][vv];
				if (seriveNodeNums[v] == 1) { // 节点v只有一个中心点
					delta[seriveNodes[v][0]] -= weight[v];
				}
				else if (seriveNodeNums[v] == 0) { // 节点v没有中心点
					for (int ll = 0; ll < coveredNodeNums[v]; ++ll) {
						NodeId l = coverages[v][ll];
						delta[l] -= weight[v];
					}
					delta[i] += weight[v];

					ucenters.erase(v);

				}

				tseriveNodes[v] = i;
				seriveNodes[v].insert(i);
				++ seriveNodeNums[v];
			}

			// swap(i, j)
			for (auto t = solution.begin(); t != solution.end(); t++) {
				if (*t == j) {
					*t = i;
					inCenter[j] = false;
					inCenter[i] = true;
					break;
				}
			}

			// 关闭j
			for (int vv = 0; vv < coveredNodeNums[j]; ++vv) {
				NodeId v = coverages[j][vv];

				seriveNodes[v].erase(j);
				-- seriveNodeNums[v];

				if (seriveNodeNums[v] == 1) { // 节点v有两个中心点
					delta[seriveNodes[v][0]] += weight[v];
					tseriveNodes[v] = seriveNodes[v][0];
				}
				else if (seriveNodeNums[v] == 0) { // 节点v只有j中心点
					for (int ll = 0; ll < coveredNodeNums[v]; ++ll) {
						NodeId l = coverages[v][ll];
						delta[l] += weight[v];
					}
					delta[j] -= weight[v];

					ucenters.insert(v);
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

		// 释放申请内存空间
		void free() {
			delete[]coveredNodeNums;
			delete[]weight;
			delete[]inCenter;
			delete[]delta;
			delete[]delta_;
			delete[]seriveNodeNums;
			delete[]tseriveNodes;
			delete[]seriveNodes;
			for (int i = 0; i < nodeNum; ++i)
				delete[]coverages[i];
		}

		int totiter = 0; // 总迭代次数
	public:
		int solve(Centers& output, PCenter& input, function<long long()> restMilliSec, int seed) {
			initRand(seed);
			
			coverAllNodesUnderFixedRadius(output, input, restMilliSec, seed);
			cerr << "Iteration : " << input.U-- << endl;

			// 缩小半径
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
				return totiter;

			cerr << input.U << ' ' << input.obj << endl;
			coverAllNodesUnderFixedRadius(output, input, restMilliSec, seed);
			
			cerr << "Tototal Iters: " << totiter << endl;
			return totiter;
		}

		bool coverAllNodesUnderFixedRadius(Centers& output, PCenter& input, function<long long()> restMilliSec, int seed) {
			nodeNum = input.nodeNum;
			centerNum = input.centerNum;
			fixNum = input.fixNodes.Num;

			Tabu_open = Tabu_close = -1;
			inCenter = new bool[nodeNum];
			seriveNodes = new doublevector[nodeNum];
			coveredNodeNums = new int[nodeNum];
			seriveNodeNums = new int[nodeNum];
			tseriveNodes = new int[nodeNum];
			delta = new int[nodeNum];
			delta_ = new int[nodeNum];
			weight = new int[nodeNum];

			coverages = new int*[nodeNum];
			for (int i = 0; i < nodeNum; ++i) {
				weight[i] = 1;
				seriveNodes[i].init(nodeNum);
				coveredNodeNums[i] = input.coveredNodeNums[i];
				coverages[i] = new int[input.coveredNodeNums[i]];

				for (int j = 0; j < coveredNodeNums[i]; ++j) 
					coverages[i][j] = input.coverages[i][j];
			}

			int SIZE = nodeNum * INTSIZE;
			ucenters.init(nodeNum);
			memset(inCenter, 0, nodeNum);
			memset(delta, 0, SIZE);
			memset(seriveNodeNums, 0, SIZE);
			memset(tseriveNodes, -1, SIZE);

			Centers solution(input.centerNum);
			initOutput(solution, input); // 初始化固定节点
			greedySearch(solution, restMilliSec); // 求解初始解

			initseriveNodes(solution); // 初始化节点被服务中心点数
			initUCenters(); // 初始化ucenter

			if (ucenters.size() == 0) {
				output = solution;
				return true;
			}

			initDelta(solution); // 初始化节点delta

			bool findSolution = false;
			Centers solution_(solution), bestsolution(solution);
			uint32_t bestUCenterNum = ucenters.size();
			uint32_t lastUCenterNum = bestUCenterNum;
			for (int iter = 1; restMilliSec() > RESSERVE; iter++) {
				// 寻找最佳交换
				pair<int, int> t = findPair(solution_);

				if (t.first == -1 && t.second == -1)
					break;

				// 交换
				makeMove(solution, t.first, t.second);

				// 比较
				if (ucenters.size() < bestUCenterNum) {
					bestsolution = solution;
					bestUCenterNum = ucenters.size();

					if (bestUCenterNum == 0) {
						output = bestsolution;
						findSolution = true;

						totiter += iter;
						break;
					}
				}
				else if (ucenters.size() >= lastUCenterNum ) {
					for (int i = 0; i < ucenters.size(); ++ i) {
						NodeId t = ucenters[i];	
						++ weight[t];

						for (int jj = 0; jj < coveredNodeNums[t]; ++jj) {
							int j = coverages[t][jj];
							++delta[j];
						}
					}
				}

				// 更新
				Tabu_open = t.first;
				Tabu_close = t.second;
				solution_ = solution;
				lastUCenterNum = ucenters.size();
			}

			free();

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

	// 
	int solvePCenter(Centers& output, PCenter& input, function<long long()> restMilliSec, int seed) {
		return Solver().solve(output, input, restMilliSec, seed);
	}

}
