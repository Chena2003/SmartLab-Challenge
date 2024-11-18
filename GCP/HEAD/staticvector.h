#ifndef _STATICVECTOR_H__
#define _STATICVECTOR_H__

#include <iostream>
#include <vector>

class staticvector {
	int *data;
	int len = 0;

public:
	staticvector() = default;
	staticvector(int n) { data = new int[n]; };

	void insert(int d) { data[len++] = d; };
	void clear() { len = 0; };
	int size() { return len; };

	int operator[](int pos) { return data[pos]; };

	~staticvector() {
		delete[] data;
	};
};

#endif

