#ifndef REGRESSION_HPP
#define REGRESSION_HPP

class UnbiasedLinearRegression { // y = k * x 
	private:
		bool is_solved;
		long double k;
	public:
		UnbiasedLinearRegression();
		void solve(long double *x, long double *y, unsigned int len); 
		long double get_k();
};

#endif
