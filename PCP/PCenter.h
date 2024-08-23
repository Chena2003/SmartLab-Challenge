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


namespace szx {

#define INTMAX 1e9 + 10
#define MODNUM 1e8

using NodeId = int;
using EdgeId = NodeId;
using Nodes = std::vector<NodeId>;
using UCenters = unordered_set<int>;


struct PCenter {
	NodeId nodeNum;
	NodeId centerNum;
	std::vector<int> coveredNodeNums; // 每个点覆盖点数
	std::vector<Nodes> coverages; // `coverages[n]` are the nodes covered by node `n` if node `n` is a center.
	std::vector<Nodes> nodesWithDrops; // `nodesWithDrops[r]` are the nodes which will drop its farthest covering node in the `r`th radius reduction.
	std::vector<vector<bool>> serives; // `serivers[n]` 表示节点n可被服务的节点
};

struct solverNodes {
	struct PCenter Nodes;
	std::vector<int> weight; // 节点的权重
	std::vector<bool> inCenter; // 节点是否被选中为中心点
	std::vector<int> delta; // 节点delta
	std::vector<int> age; // 节点年龄
	UCenters ucenters; // 当前尚未覆盖节点
	std::vector<int> seriveNodeNums; // 每个点被服务的中心点数
	std::vector<int> seriveNodeId; // 当节点被服务中心点数为1时，表示为被服务的中心点；否则为-1
};


using Centers = Nodes; // `Centers[k]` is the `k`th picked center.
using Deltas = Nodes;


void solvePCenter(Centers& output, PCenter& input, std::function<long long()> restMilliSec, int seed);

}


#endif // CN_HUST_SZX_NPBENCHMARK_P_CENTER_H
