#ifndef __DOUBLEVECTOR_H__
#define __DOUBLEVECTOR_H__

#include <iostream>
#include <cstdint>

class doublevector {
private:
	const int INTSIZE = sizeof(int);

	int* map;
	int* data;
	int len = 0;
	int mapsize;

public:
	doublevector() = default;
	doublevector(const doublevector& d);
	doublevector(int n);
	~doublevector();

	void init(int n);
	void insert(int d);
	void erase(int d);
	uint32_t size();
	bool find(int d);
	void clear();

	int operator[](int pos)	const { return data[pos]; };
};

#endif // __DOUBLEVECTOR_H__