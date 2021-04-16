#include <assert.h>

#include "Regression.hpp"

// UnbiasedLinearRegression

UnbiasedLinearRegression::UnbiasedLinearRegression() {
	is_solved = false;
}

// loss = \sum_{i=1}{n} (k * x_i - y_i)^2 
// d(loss)/d(k) = \sum_{i=1}{n} 2(k * x_i - y_i)  * x_i
// 		= \sum_{i=1}{n} 2 * x_i^2 * k - \sum_{i=1}{n} 2 * x_i y_i
// 		= 0
// =>
// k = \frac{\sum_{i=1}{n}x_i*y_i}{\sum_{i=1}{n}x_i^2}
void UnbiasedLinearRegression::solve(long double *x, long double *y, unsigned int len) {
	long double tmp_0 = 0;
	for (unsigned int i = 0; i < len; ++ i) {
		tmp_0 += x[i] * y[i];
	}
	long double tmp_1 = 0;
	for (unsigned int i = 0; i < len; ++ i) {
		tmp_1 += x[i] * x[i];
	}
	k = tmp_0 / tmp_1;
	is_solved = true;
}

long double UnbiasedLinearRegression::get_k() {
	assert(is_solved);
	return k;
}
