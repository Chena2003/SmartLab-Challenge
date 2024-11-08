#include "head.h"
#include <iostream>
#include <array>
#include "doublevector.h"
#include <cstring>

using namespace std;

head::head(int input_num_vertex, int input_edge_num, int input_num_color, int tabu_max_iter, vector<array<int, 2>>& input_edges, int input_seed) :
	NodeNum(input_num_vertex),
	EdgeNum(input_edge_num),
	ColorNum(input_num_color),
	MAXITER(tabu_max_iter),
	seed(input_seed),
	edges(input_edges)
{
	COPYSIZE = sizeof(int) * NodeNum;
	SETSIZE = sizeof(int) * ColorNum;
	pseudoRandNumGen = mt19937(input_seed);

	ConflictNodes = new int[NodeNum];
	ConflictNodePos = new int[NodeNum+1];
	NeighborNums = new int[NodeNum];
	EqualDeltaU = new int[4000];
	EqualDeltaC = new int[4000];

	TabuTable = new int* [NodeNum];
	ACTable = new int* [NodeNum];
	NeighborTable = new int* [NodeNum];
	for (int i = 0; i < NodeNum; ++i) {
		TabuTable[i] = new int[ColorNum];
		ACTable[i] = new int[ColorNum];
		NeighborTable[i] = new int[NodeNum];
	}

	p1 = new int[NodeNum];
	p2 = new int[NodeNum];
	e1 = new int[NodeNum];
	e2 = new int[NodeNum];
	best = new int[NodeNum];

	memset(NeighborNums, 0, COPYSIZE);
	for (const auto& edge : input_edges) {
		int i = edge[0], j = edge[1];

		NeighborTable[i][NeighborNums[i]] = j;
		NeighborTable[j][NeighborNums[j]] = i;

		++NeighborNums[i];
		++NeighborNums[j];
	}
}

head::~head() {
	delete p1;
	delete p2;
	delete e1;
	delete e2;
	delete best;
	delete ConflictNodes;
	delete ConflictNodePos;
	delete NeighborNums;
	delete EqualDeltaU;
	delete EqualDeltaC;

	for (int i = 0; i < NodeNum; ++i) {
		delete TabuTable[i];
		delete ACTable[i];
		delete NeighborTable[i];
	}
}

void head::InitNodeColor(NodeColors sol) {
	for (int i = 0; i < NodeNum; ++i) 
		sol[i] = rand(ColorNum);
}

int head::ComputeConflict(NodeColors sol) {
	int f = 0;
	for (int i = 0; i < NodeNum; ++i) {
		int ci = sol[i];
		int adjNum = NeighborNums[i];

		for (int j = 0; j < adjNum; ++j) {
			int adj = NeighborTable[i][j];
			if (ci == sol[adj])
				++f;
		}
	}

	return f / 2;
}

void head::GPX(NodeColors p1, NodeColors p2, NodeColors son) {
	// ÿ����ɫ�Ľڵ㼯��
	NodeColors p[2] = { p1, p2 };
	vector<vector<doublevector>> pnode(2, vector<doublevector>(ColorNum, doublevector(NodeNum)));
	doublevector uncolorNodes(NodeNum); // δȾɫ�ڵ�

	for (int i = 0; i < NodeNum; ++i) {
		// ��ʼ��ÿ����ɫ�Ľڵ�
		pnode[0][p1[i]].insert(i);
		pnode[1][p2[i]].insert(i);

		// ��ʼ��δȾɫ�ڵ�
		uncolorNodes.insert(i);
	}

	// Ͱ�����ҳ��ڵ�����������ɫ
	vector<vector<doublevector>> cnodeNum(2, vector<doublevector>(NodeNum, doublevector(ColorNum)));
	uint32_t maxcolorNum[2] = { 0, 0 };
	for (int i = 0; i < ColorNum; ++i) {
		cnodeNum[0][pnode[0][i].size()].insert(i);
		cnodeNum[1][pnode[1][i].size()].insert(i);

		if (pnode[0][i].size() > maxcolorNum[0])
			maxcolorNum[0] = pnode[0][i].size();
		if (pnode[1][i].size() > maxcolorNum[1])
			maxcolorNum[1] = pnode[1][i].size();
	}

	// �Ӹ�ĸ��ѡȡk����ɫ
	for (int k = 0, id = 0; k < ColorNum; ++k, id ^= 1) {
		doublevector t;
		while (cnodeNum[id][maxcolorNum[id]].size() == 0)
			--maxcolorNum[id];

		int choice = rand(cnodeNum[id][maxcolorNum[id]].size());
		int color = cnodeNum[id][maxcolorNum[id]][choice];

		t = pnode[id][color];

		// ɾ����ɫt�����нڵ�
		for (int i = 0; i < t.size(); ++i) {
			int nodeId = t[i];
			son[nodeId] = k;

			int tid = id ^ 1;
			int tcolor = p[tid][nodeId];
			cnodeNum[tid][pnode[tid][tcolor].size()].erase(tcolor);
			pnode[tid][tcolor].erase(nodeId);
			cnodeNum[tid][pnode[tid][tcolor].size()].insert(tcolor);

			uncolorNodes.erase(nodeId);
		}
		cnodeNum[id][maxcolorNum[id]].erase(color);
		cnodeNum[id][0].insert(color);
	}

	// δ������ɫ�ڵ��������ɫ
	for (int i = 0; i < uncolorNodes.size(); ++i) {
		int nodeId = uncolorNodes[i];
		son[nodeId] = rand(ColorNum);
	}
}

void head::FindMove(NodeColors sol, int &f, int iter) {
	delta = 10000;
	int EqualDeltaLen = 0;
	for(int i = 0; i < ConflictNodeLen; ++i) {
		int ConflictNode = ConflictNodes[i];
		int ci = sol[ConflictNode];

		for (int k = 0; k < ColorNum; ++k) {
			if (ci == k)
				continue;

			int t = ACTable[ConflictNode][k] - ACTable[ConflictNode][ci];

			if (TabuTable[ConflictNode][k] > iter && f + t >= bestf)
				continue;

			if (t < delta) {
				delta = t;

				EqualDeltaLen = 0;
				EqualDeltaU[EqualDeltaLen] = ConflictNode;
				EqualDeltaC[EqualDeltaLen] = k;
				++EqualDeltaLen;
			}
			else if (t == delta) {
				EqualDeltaU[EqualDeltaLen] = ConflictNode;
				EqualDeltaC[EqualDeltaLen] = k;
				++EqualDeltaLen;
			}
		}
	}

	int choice = rand(EqualDeltaLen);
	bestu = EqualDeltaU[choice];
	bestc = EqualDeltaC[choice];
}

void head::MakeMove(NodeColors sol, int &f, int iter) {
	int ci = sol[bestu];
	sol[bestu] = bestc;
	f += delta;

	// ���������ʷ��ֵ
	if (f < besthistoryf)
		besthistoryf = f;

	TabuTable[bestu][ci] = iter + 0.6 * f + rand(10); // ���½��ɱ�

	// ������ɫ��
	int adjNum = NeighborNums[bestu];
	for (int i = 0; i < adjNum; ++i) {
		int adj = NeighborTable[bestu][i];
		
		--ACTable[adj][ci];
		++ACTable[adj][bestc];

		// ���³�ͻ�ڵ�
		if (ACTable[adj][sol[adj]] == 0 && ConflictNodePos[adj] != -1) {
			int t = ConflictNodePos[adj];
			ConflictNodes[t] = ConflictNodes[--ConflictNodeLen];
			ConflictNodePos[ConflictNodes[t]] = t;
			ConflictNodePos[adj] = -1;
		}
		else if (ACTable[adj][sol[adj]] != 0 && ConflictNodePos[adj] == -1) {
			ConflictNodes[ConflictNodeLen] = adj;
			ConflictNodePos[adj] = ConflictNodeLen;
			++ConflictNodeLen;
		}
	}
}

int head::TabuCol(NodeColors sol) {
	ConflictNodeLen = 0;
	memset(ConflictNodePos, 0xff, COPYSIZE + sizeof(int));
	for (int i = 0; i < NodeNum; ++i) {
		memset(TabuTable[i], 0, SETSIZE);
		memset(ACTable[i], 0, SETSIZE);
	}

	int f = 0;
	// ��ʼ����ɫ�����ھӱ��ͳ�ͻ�ڵ�
	for(const auto& edge : edges) {
		int i = edge[0], j = edge[1];
		int ci = sol[i], cj = sol[j];

		++ACTable[i][cj];
		++ACTable[j][ci];

		if (ci == cj) {
			++f;

			if (ConflictNodePos[i] == -1) {
				ConflictNodes[ConflictNodeLen] = i;
				ConflictNodePos[i] = ConflictNodeLen;
				++ConflictNodeLen;
			}
			if (ConflictNodePos[j] == -1) {
				ConflictNodes[ConflictNodeLen] = j;
				ConflictNodePos[j] = ConflictNodeLen;
				++ConflictNodeLen;
			}
		}
	}

	besthistoryf = f;

	for (int iter = 0; iter < MAXITER; ++iter) {
		FindMove(sol, f, iter);

		MakeMove(sol, f, iter);

		if (f == 0) 
			break;
	}

	return f;
}

bool head::check(NodeColors p1, NodeColors p2) {
	vector<vector<int>> bigraph(ColorNum, vector<int>(ColorNum, 0));
	vector<doublevector> pnode(ColorNum, doublevector(NodeNum));
	for (int i = 0; i < NodeNum; ++i)
		pnode[p1[i]].insert(i);

	for (int k = 0; k < ColorNum; ++k) {
		for (int i = 0; i < pnode[k].size(); ++i) {
			int nodeId = pnode[k][i];
			++bigraph[k][p2[nodeId]];
		}
	}

	int sum = 0;
	vector<bool> rvisited(ColorNum, false), cvisited(ColorNum, false);
	for (int k = 0; k < ColorNum; ++k) {
		int max = 0, maxr, maxc;
		int r, c;
		for (r = 0; r < ColorNum; ++r)
			for (c = 0; c < ColorNum; ++c)
				if (!rvisited[r] && !cvisited[c] && bigraph[r][c] >= max) {
					max = bigraph[r][c];
					maxr = r;
					maxc = c;
				}
		sum += max;
		rvisited[maxr] = true;
		cvisited[maxc] = true;
	}

	return (sum + 10) >= NodeNum;
}

void head::HybridEvolution() {
	InitNodeColor(p1);
	InitNodeColor(p2);
	InitNodeColor(e2);
	InitNodeColor(best);

	bestf = ComputeConflict(best);
	e1f = ComputeConflict(e1);
	e2f = ComputeConflict(e2);

	int cycle = 10; // ÿ10������
	NodeColors c1 = new int[NodeNum];
	NodeColors c2 = new int[NodeNum];
	int gen;
	for (gen = 0; bestf > 0; ++gen) {
		GPX(p1, p2, c1);
		GPX(p2, p1, c2);
		
		memcpy(p1, c1, COPYSIZE);
		memcpy(p2, c2, COPYSIZE);

		p1f = TabuCol(p1);
		p2f = TabuCol(p2);

		NodeColors tmp = p1;
		int t = p1f;
		if (p1f > p2f) {
			tmp = p2;
			t = p2f;
		}
		if (t < e1f) {
			e1f = t;
			memcpy(e1, tmp, COPYSIZE);
		}

		if (e1f < bestf) {
			memcpy(best, e1, COPYSIZE);
			bestf = e1f;
		}

		if (gen % cycle == 0) {
			memcpy(p1, e2, COPYSIZE);
			memcpy(e2, e1, COPYSIZE);
			InitNodeColor(e1);
			e1f = ComputeConflict(e1);
		}

		if (check(p1, p2)) {
			cerr << "Gen: " << gen << " p1 == p2" << endl;
			pseudoRandNumGen = mt19937(seed + rand(10000));
			InitNodeColor(p2);
		}
	}

	cerr << "Gen: " << gen << endl;
}

void head::WriteSolution(vector<int>& output) {
	for (int i = 0; i < NodeNum; ++i) {
		output[i] = best[i];
	}
}


