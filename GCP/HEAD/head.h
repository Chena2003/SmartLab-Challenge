#ifndef __HEAD_H__
#define __HEAD_H__

#include <vector>
#include <random>
#include <iostream>
#include <algorithm>
#include "staticvector.h"

using NodeColors = int*;

class head {
private:
	int NodeNum, EdgeNum, ColorNum;
	int MAXITER; // ���ɵ�������
	int COPYSIZE;
	int SETSIZE;
	int seed;

	std::vector<std::array<int, 2>> edges;
	//tdarray TabuTable;
	//tdarray ACTable;
	int** TabuTable; // ���ɱ�
	int** ACTable; // �ھ���ɫ��
	int** NeighborTable; // �ڵ��ھӱ�
	int* NeighborNums; // �ڵ��ھӵĸ���

	int* ConflictNodes;
	int* ConflictNodePos;
	int ConflictNodeLen;

	int besthistoryf; // ��ʷ��õ�f
	int* EqualDeltaU; // ��ö����ڵ�
	int* EqualDeltaC; // ��ö�����ɫ
	int bestc, bestu, delta;

	NodeColors p1, p2, e1, e2, best; // ���
	int p1f, p2f, e1f, e2f, bestf; // f

	std::mt19937 pseudoRandNumGen; // �����

	inline int rand(int u) { return pseudoRandNumGen() % u; };
	void InitNodeColor(NodeColors sol);
	int ComputeConflict(NodeColors sol);
	void GPX(NodeColors p1, NodeColors p2, NodeColors son);
	void FindMove(NodeColors sol, int &f, int iter);
	void MakeMove(NodeColors sol, int &f, int iter);
	int TabuCol(NodeColors sol);
	bool check(NodeColors p1, NodeColors p2);
public:
	head() = default;
	head(int input_num_vertex, int input_edge_num, int input_num_color, int tabu_max_iter, std::vector<std::array<int, 2>>& input_edges, int input_seed);
	~head();

	void HybridEvolution();
	void WriteSolution(std::vector<int>& output);
};




#endif
