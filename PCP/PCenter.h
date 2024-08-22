////////////////////////////////
/// usage : 1.	SDK for p-center solver.
/// 
/// note  : 1.	
////////////////////////////////

#ifndef CN_HUST_SZX_NPBENCHMARK_P_CENTER_H
#define CN_HUST_SZX_NPBENCHMARK_P_CENTER_H


#include <vector>
#include <functional>


namespace szx {

using NodeId = int;
using EdgeId = NodeId;
using Nodes = std::vector<NodeId>;

struct PCenter {
	NodeId nodeNum;
	NodeId centerNum;
	std::vector<int> coveredNodeNums; // 每个点覆盖点数
	std::vector<int> weight; // 每个点的权重
	std::vector<Nodes> coverages; // `coverages[n]` are the nodes covered by node `n` if node `n` is a center.
	std::vector<Nodes> nodesWithDrops; // `nodesWithDrops[r]` are the nodes which will drop its farthest covering node in the `r`th radius reduction.
};

using Centers = Nodes; // `Centers[k]` is the `k`th picked center.
using Deltas = Nodes;


void solvePCenter(Centers& output, PCenter& input, std::function<long long()> restMilliSec, int seed);

}


#endif // CN_HUST_SZX_NPBENCHMARK_P_CENTER_H
