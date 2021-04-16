#ifndef RANDOM_HPP
#define RANDOM_HPP

#include <random>

class Random {
	private:
		std::random_device rd;
		std::mt19937 gen;
		std::uniform_int_distribution<unsigned long long> dis;
	public:
		Random(): gen(rd()), dis(0, (unsigned long long)-1) {
		}
		unsigned long long rand_int(); // this is thread-safe
		unsigned long long rand_int(unsigned long long min, unsigned long long max); // thread-safe
};

#endif
