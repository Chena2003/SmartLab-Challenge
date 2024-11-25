#include "doublevector.h"
#include <cstring>

// 复制构造函数
doublevector::doublevector(const doublevector& d) : len(d.len), mapsize(d.mapsize) {
	map = new int[mapsize + 1];
	data = new int[mapsize];

	memcpy(map, d.map, INTSIZE * (mapsize + 1));
	memcpy(data, d.data, INTSIZE * mapsize);
}

// 构造函数
doublevector::doublevector(int n) : mapsize(n) {
	map = new int[mapsize + 1];
	data = new int[mapsize];

	memset(map, 0xff, INTSIZE * (mapsize + 1));
	memset(data, 0, INTSIZE * (mapsize));
}

doublevector::~doublevector() {
	delete map;
	delete data;
}

void doublevector::init(int n) {
	mapsize = n;
	map = new int[mapsize + 1];
	data = new int[mapsize];

	memset(map, 0xff, INTSIZE * (mapsize + 1));
	memset(data, 0, INTSIZE * (mapsize));
}

void doublevector::insert(int d) {
	if (map[d] == -1) {
		data[len] = d;
		map[d] = len;
		++len;
	}
}

void doublevector::erase(int d) {
	if (map[d] == -1)
		return;

	int pos = map[d];
	data[pos] = data[len - 1];
	map[data[pos]] = pos;
	-- len;
	map[d] = -1;
}

uint32_t doublevector::size() {
	return len;
}

bool doublevector::find(int d) {
	return map[d] != -1;
}

void doublevector::clear() {
	memset(map, 0xff, INTSIZE * (mapsize + 1));
	memset(data, 0, INTSIZE * mapsize);
	len = 0;
}