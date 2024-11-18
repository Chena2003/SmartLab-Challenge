#include "doublevector.h"

void doublevector::init(int n) {
	//map.resize(n + 1);
	map = std::vector<int>(n + 1, -1);
}

void doublevector::insert(int d) {
	if (map[d] == -1) {
		data.push_back(d);
		map[d] = data.size() - 1;
	}
}

void doublevector::erase(int d) {
	if (map[d] != -1) {
		int pos = map[d];
		data[pos] = data.back();
		map[data[pos]] = pos;
		data.pop_back();
		map[d] = -1;
	}
}

uint32_t doublevector::size() {
	return data.size();
}