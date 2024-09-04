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
#include <boost/dynamic_bitset.hpp>

namespace szx {

#define INTMAX 1e9 + 10
#define MODNUM 1e8

	using NodeId = int;
	using EdgeId = NodeId;
	using Nodes = std::vector<NodeId>;
	using Flags = std::vector<bool>;
	using UCenters = std::unordered_set<int>;
	using dynamic_bitset = boost::dynamic_bitset<>;

	struct SNodes { // 取代unordered_set
		NodeId Num; // 数量
		dynamic_bitset Nodes; // 点集合
	};

	struct PCenter {
		NodeId nodeNum;
		NodeId centerNum;
		std::vector<int> coveredNodeNums; // 每个点覆盖点数
		std::vector<Nodes> coverages; // `coverages[n]` are the nodes covered by node `n` if node `n` is a center.
		std::vector<Nodes> nodesWithDrops; // `nodesWithDrops[r]` are the nodes which will drop its farthest covering node in the `r`th radius reduction.
		std::vector<dynamic_bitset> serives; // `serivers[i][j]` 表示节点i可被节点j服务

		SNodes fixNodes;
	};

	struct solverNodes {
		struct PCenter Nodes;
		std::vector<int> weight; // 节点的权重
		std::vector<bool> inCenter; // 节点是否被选中为中心点
		std::vector<int> delta; // 节点delta
		std::vector<int> age; // 节点年龄
		std::unordered_set<int> ucenters; // 当前尚未覆盖节点
		std::vector<std::unordered_set<int>> seriveNodes; // 节点被服务的中心点集合
	};


	using Centers = Nodes; // `Centers[k]` is the `k`th picked center.
	using Deltas = Nodes;


	void solvePCenter(Centers& output, PCenter& input, std::function<long long()> restMilliSec, int seed);

}


#endif // CN_HUST_SZX_NPBENCHMARK_P_CENTER_H
