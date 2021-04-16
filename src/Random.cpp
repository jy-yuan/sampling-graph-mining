#include <thread>
#include <random>

#include "Random.hpp"
#include "Type.hpp"

unsigned long long Random::rand_int() {
	return dis(gen);
}

unsigned long long Random::rand_int(unsigned long long min, unsigned long long max) {
	std::uniform_int_distribution<unsigned long long> distribution(min,max);
	return distribution(gen);
}
