#include <iostream>

typedef double(*Haar)(int, int, int);
typedef double(*signal)(double n);

inline double __S(double _n) { return sin(_n); } //looks like discret sinus


int n = 10;
int N = 1024;

double _phi_mkn(int m, int k, int n) //phi_m,k(n)
{
	double _sqrt = 1. / sqrt(pow(2, m));
	double _t = pow(2, -m) * n - k;
	int phi;

	if (_t >= 0.  && _t <= 0.5) phi =  1;
	if (_t >= 0.5 && _t <= 1.)  phi = -1;
	if (_t < 0.	  || _t >  1.)  phi =  0;

	return _sqrt * phi;
}

double _C_mk(int _N, int m, int k, Haar _h, signal _s) {
	double res = 0.;

	for (int n(0); n < _N; ++n) 
	{
		res += __S(n / _N) * _phi_mkn(m, k, n);
	}
	return res;
}

double _approx(Haar h, signal s) {
	double _res = 0.;
	for(int m(0); m < n; m++)
	{
		for (int k(0); k < pow(2, n - m); k++)
		{
			_res += _C_mk(N, m, k, h, s) * _phi_mkn(m, k, n);
		}
	}
	return _res; //finally answer
}

int main(int argc, void* argv[]) 
{

	return 0;
}