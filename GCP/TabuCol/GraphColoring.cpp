#include "GraphColoring.h"
//#include "doublevector.h"

#include <random>
#include <iostream>
#include <vector>


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
		double frand(double ub) { return uniform_real_distribution<double>(0.0, ub)(pseudoRandNumGen); };

	public:
		void solve(NodeColors& output, GraphColoring& input, function<long long()> restMilliSec, int seed) {
			initRand(seed);

			uint32_t f = 0;
			vector<vector<int>> TabuTable(input.nodeNum, vector<int>(input.colorNum));
			vector<vector<int>> ACTable(input.nodeNum, vector<int>(input.colorNum));
			vector<vector<int>> NeighborTable(input.nodeNum);

			vector<int> conflictNodes; // 冲突节点
			vector<int> pos(input.nodeNum + 1, -1);

			// 初始化节点颜色
			for (int i = 0; i < input.nodeNum; ++i) {
				output[i] = rand(input.colorNum);
			}
			
			// 初始化ACTable和f
			for (auto edge : input.edges) {
				int i = edge[0], j = edge[1];

				NeighborTable[i].push_back(j);
				NeighborTable[j].push_back(i);

				++ ACTable[i][output[j]];
				++ ACTable[j][output[i]];

				if (output[i] == output[j]) {
					++ f; // 计算评价函数

					if (pos[i] == -1) {
						conflictNodes.push_back(i);
						pos[i] = conflictNodes.size() - 1;
					}

					if (pos[j] == -1) {
						conflictNodes.push_back(j);
						pos[j] = conflictNodes.size() - 1;
					}
				}
			}

			if (f == 0)
				return;

			const int RESERVE = 1; // 预留输出时间
			int bestf = f;
			int cnt = 0;
			int iter;
			for (iter = 0; (restMilliSec() > RESERVE); ++iter) {
				// 寻找最佳关键动作
				vector<int> bestu(conflictNodes.size() * input.colorNum), bestc(conflictNodes.size() * input.colorNum);
				vector<int> tbestu(conflictNodes.size() * input.colorNum), tbestc(conflictNodes.size() * input.colorNum);
				int cnt = 0, tcnt = 0;
				int delta = 100000, tdelta = 100000;

				for(auto i : conflictNodes) {
					for (int k = 0; k < input.colorNum; ++k) {
						if (output[i] == k)
							continue;

						int t = ACTable[i][k] - ACTable[i][output[i]];

						if (TabuTable[i][k] > iter) {
							if (t < tdelta) {
								tdelta = t;

								tcnt = 0;
								tbestu[tcnt] = i;
								tbestc[tcnt] = k;
								++ tcnt;
							}
							else if (t == tdelta) {
								tbestu[tcnt] = i;
								tbestc[tcnt] = k;
								++tcnt;
							}
						}
						else {
							if (t < delta) {
								delta = t;

								cnt = 0;
								bestu[cnt] = i;
								bestc[cnt] = k;
								++cnt;
							}
							else if (t == delta) {
								bestu[cnt] = i;
								bestc[cnt] = k;
								++cnt;
							}
						}
					}
				}

				int u, ci, cj;
				// 禁忌动作是否满足特赦准则
				if (tdelta + f < bestf && tdelta < delta) {
					int choice = rand(tcnt);

					u = tbestu[choice];
					ci = output[u];
					cj = tbestc[choice];
					delta = tdelta;
				}
				else {
					int choice = rand(cnt);

					u = bestu[choice];
					ci = output[u];
					cj = bestc[choice];
				}

				// 执行最佳关键动作
				output[u] = cj;
				f += delta;
				// 找到合法着色方案
				if (f == 0) {
					break;
				}
					 
				if (f < bestf)
					bestf = f;

				TabuTable[u][ci] = iter + f + rand(10); // 更新禁忌表
				// 更新颜色表
				for (auto v : NeighborTable[u]) {
					-- ACTable[v][ci];
					++ ACTable[v][cj];

					// 更新冲突节点
					if (ACTable[v][output[v]] == 0 && pos[v] != -1) {
						int t = pos[v];
						conflictNodes[t] = conflictNodes.back();
						pos[conflictNodes[t]] = t;
						conflictNodes.pop_back();
						pos[v] = -1;
					}
					else if (ACTable[v][output[v]] != 0 && pos[v] == -1) {
						conflictNodes.push_back(v);
						pos[v] = conflictNodes.size() - 1;
					}
				}
			}
		
			cerr << "Iter: " << iter << endl;
		}
	};

	// solver.
	void solveGraphColoring(NodeColors& output, GraphColoring& input, function<long long()> restMilliSec, int seed) {
		Solver().solve(output, input, restMilliSec, seed);
	}

}
