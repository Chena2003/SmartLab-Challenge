#ifndef __DOUBLEVECTOR_H__
#define __DOUBLEVECTOR_H__

#include <vector>
#include <cstdint>

class doublevector {
private:
	int mapsize = 0;
	std::vector<int> map;
	std::vector<int> data;

public:
	doublevector() = default;
	doublevector(const doublevector& d) : mapsize(d.mapsize), map(d.map), data(d.data) {};
	doublevector(int n) : mapsize(n), map(n, -1) {};

	void init(int n);
	void insert(int d);
	void erase(int d);
	uint32_t size();
	bool find(int d);
	void clear();

	int operator[](int pos)	const { return data[pos]; };
};

#endif // __DOUBLEVECTOR_H__