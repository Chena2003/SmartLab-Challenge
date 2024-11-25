#ifndef __WORK_H__
#define __WORK_H__

#include <iostream>
#include <unordered_set>
#include <vector>
#include <random>
#include <functional>
#include <fstream>
#include "doublevector.h"
#include "utils.h"
#include "boost/dynamic_bitset.hpp"

using dynamic_bitset = boost::dynamic_bitset<>;

class Work {
private:
	const int RT0 = 10, MAXRT = 15, MAXACCU = 1000;

	int VNum, CNum; // 节点数量和颜色数量
    std::vector<dynamic_bitset> ColorDomain; // 每个节点的颜色域
	std::vector<doublevector> NeighborTable; // 每个节点的邻居
	std::vector<doublevector> RowNeighborTable, ColNeighborTable; // 节点的行邻居和列邻居
	std::vector<bool> FixNodes; // 固定节点集或删除节点集
	std::vector<doublevector> RowColorDomain; // 行颜色域

	// 解结构
	std::vector<int> Sol, BestSol;
	int conflict;

	// 禁忌搜索
	int ConflictNodeLen = 0;
	std::vector<int> ConflictNodes; // 冲突节点
	std::vector<int> ConflictNodePos; 
	std::vector<std::vector<int>> AdjColorTable; // 邻居颜色表
	std::vector<std::vector<int>> TabuTable; // 禁忌表

	int besthistoryf; // 历史最好的f
	std::vector<int> EqualNontabuDeltaU, EqualNontabuDeltaV; // 最好非禁忌动作节点U和V
	std::vector<int> EqualTabuDeltaU, EqualTabuDeltaV; // 最好禁忌动作节点U和V
	int bestu, bestv, delta; // 最好动作

	// 随机数
	std::mt19937 pseudoRandNumGen;
	void initRand(int seed) { pseudoRandNumGen = std::mt19937(seed); }
	int fastRand(int lb, int ub) { return (pseudoRandNumGen() % (ub - lb)) + lb; }
	int fastRand(int ub) { return pseudoRandNumGen() % ub; }
	int rand(int lb, int ub) { return std::uniform_int_distribution<int>(lb, ub - 1)(pseudoRandNumGen); }
	int rand(int ub) { return std::uniform_int_distribution<int>(0, ub - 1)(pseudoRandNumGen); }


	void InitSol();
	int CheckRule(int Node);
	void FindMove(int iter);
	void MakeMove(int iter);
public:
	Work() = default;
	Work(int vnum, int cnum, const std::vector<Assignment>& fixNode, int seed);
	
	int solve(std::function<long long()> restMilliSec);
	void savesol(std::ostream& os);
	void check();
};




#endif // __WORK_H__
