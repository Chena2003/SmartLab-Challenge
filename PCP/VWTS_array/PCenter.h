////////////////////////////////
/// usage : 1.	SDK for p-center solver.
/// 
/// note  : 1.	
////////////////////////////////

#ifndef CN_HUST_SZX_NPBENCHMARK_P_CENTER_H
#define CN_HUST_SZX_NPBENCHMARK_P_CENTER_H


#include <vector>
#include <functional>
#include <unordered_set>
#include <set>
//#include <boost/dynamic_bitset.hpp>
#include <cstring>

#include "doublevector.h"

namespace szx {

#define INTMAX 1e9 + 10
#define MODNUM 1e8

	using NodeId = int;
	using EdgeId = NodeId;
	using Nodes = std::vector<NodeId>;
	using Flags = std::vector<bool>;
	//using UCenters = std::unordered_set<int>;
	using UCenters = std::set<int>;
	//using dynamic_bitset = boost::dynamic_bitset<>;

	struct SNodes { // 取代unordered_set
		NodeId Num; // 数量
		Flags Nodes; // 点集合
	};

	struct PCenter {
		int U, L, obj;
		NodeId nodeNum;
		NodeId centerNum;
		std::vector<int> coveredNodeNums; // 每个点覆盖点数
		std::vector<Nodes> coverages; // `coverages[n]` are the nodes covered by node `n` if node `n` is a center.
		std::vector<Nodes> nodesWithDrops; // `nodesWithDrops[r]` are the nodes which will drop its farthest covering node in the `r`th radius reduction.

		SNodes fixNodes;
	};

	struct solverNodes {
		struct PCenter &Nodes;
		std::vector<int> weight; // 节点的权重
		std::vector<bool> inCenter; // 节点是否被选中为中心点
		//std::vector<int> delta; // 节点delta
		int* delta;
		std::vector<int> age; // 节点年龄
		//UCenters ucenters; // 当前尚未覆盖节点
		//std::vector<UCenters> seriveNodes; // 节点被服务的中心点集合

		//doublevector ucenters_test;
		doublevector ucenters;
		std::vector<doublevector> seriveNodes; // 节点被服务的中心点集合

		//std::vector<int> ucenters; // 当前尚未覆盖节点
		std::vector<int> seriveNodeNums; // 节点被服务的中心点数量
		std::vector<int> tseriveNodes; // 节点被服务的唯一中心点 

		//solverNodes(struct PCenter& input) :
		//	Nodes(input),
		//	inCenter(input.nodeNum, 0),
		//	seriveNodes(input.nodeNum, doublevector(input.nodeNum)),
		//	ucenters(input.nodeNum),
		//	weight(input.nodeNum, 1),
		//	seriveNodeNums(input.nodeNum, 0) {
		//	delta = new int[input.nodeNum];
		//	memset(delta, 0, input.nodeNum * sizeof(int));
		//};

		solverNodes(struct PCenter& input) : Nodes(input) {};

		//~solverNodes() {
		//	delete[]delta;
		//};
	};


	using Centers = Nodes; // `Centers[k]` is the `k`th picked center.
	using Deltas = Nodes;


	int solvePCenter(Centers& output, PCenter& input, std::function<long long()> restMilliSec, int seed);

}


#endif // CN_HUST_SZX_NPBENCHMARK_P_CENTER_H
