#include <stdio.h>

#include "Regression.hpp"

int main(int argc, char ** argv) {
	const unsigned int len = 7;
	long double true_error[len] = {0.09, 0.26, 0.15, 0.42, 0.57, 0.10, 2.62};
	long double expected_error[len] = {1, 2, 4, 8, 16, 32, 64};
	for (unsigned int i = 0; i < len; ++ i) {
		true_error[i] /= 100.;
		expected_error[i] /= 100.;
	}
	UnbiasedLinearRegression * regression = new UnbiasedLinearRegression();
	regression->solve(true_error, expected_error, len);
	double k = regression->get_k();
	k *= k;

	printf("%.4f\n", k);

	delete regression;
	return 0;
}
