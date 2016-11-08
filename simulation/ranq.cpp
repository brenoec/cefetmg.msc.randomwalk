#ifndef __RANQ2__
#define __RANQ2__

#include <iostream>
#include <limits>
#include <bitset>
#include <math.h>

using namespace std;

class Ranq2 {

	private:
		long long unsigned int _v = 4101842887655102017LL;
		long long unsigned int _w = 1;

	public:
		long long unsigned int next();
		double nextd();
		bool nextb();
		int nexti();

		Ranq2(long long unsigned int seed = 123456789) {
			_v ^= seed;
			_w = next();
			_v = next();
		}
};

long long unsigned int Ranq2::next() {
	_v ^= _v >> 17;
	_v ^= _v << 31;
	_v ^= _v >>  8;

	_w = 4294957665U * (_w & 0xffffffff) + (_w >> 32);

	return _v ^ _w;
}

bool Ranq2::nextb() {
	return (bool) (next() % 2);
}

int Ranq2::nexti() {
	return (int) next();
}

double Ranq2::nextd() {
	return 5.42101086242752217e-20 * next();
}

#endif
