#include "doublevector.h"

void doublevector::init(int n) {
	map.resize(n + 1);
}

void doublevector::insert(int d) {
	data.push_back(d);
	map[d] = data.size() - 1;
}

void doublevector::erase(int d) {
	int pos = map[d];
	data[pos] = data.back();
	map[data[pos]] = pos;
	data.pop_back();
	map[d] = -1;
}

uint32_t doublevector::size() {
	return data.size();
}