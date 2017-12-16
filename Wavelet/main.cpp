#include <iostream>
#include "Wavelet_transform.h"
#include <cmath>

typedef double(*func)	(double t);
typedef double(*_func)	(func, double, int, int);
typedef double(*coef)	(func, _func, double, double, int, int);

int m;//scale
int k;//counter
int N = 8;//amount of discret counts
double step = 1. / N;
double C[3][8];//array of coefs of direct conversion
double sd[3];//discret recovered signal

inline double _input_signal(double _t) // continious signal for transformation
{
	double _res;
	(_t >= 0 && _t <= 1) ? _res = abs(sin(10 * _t) + 5 * _t) : _res = -1.;
	return _res;
}

inline double _phi_0(double _t) //Haar scaling
{
	int _res;
	(_t >= 0 && _t <= 1) ? _res = 1 : _res = 0;
	return _res;
}

double _scaling_func(func _f, double _t, int _m, int _k)
{
	double _tmp = pow(2, _m) * _t - _k;
	return pow(2, 0.5 * _m) * _f(_tmp);
}

void _calc_coefs(func _phi_0, func f, _func _f, double _step, int _m) //C(m,k) and sd[k]
{
	double underint_a, underint_b;
	double _step_int = _step / 128;
	double _trap;
	for (int m(3); m < _m + 1; ++m)
	{
		int k(0);
		for (; k < pow(2, m); ++k)
		{
			double _res = 0.;

			for (int i(0); i < 128; ++i)
			{
				underint_a = f(k * _step + i * _step_int) 
							* _f(_phi_0, k * _step + i * _step_int, m, k);
				underint_b = f(k * _step + (i + 1) * _step_int)
					* _f(_phi_0, k * _step + (i + 1) * _step_int, m, k);
				_trap = 0.5 * (underint_b + underint_a) * _step_int;
				_res += _trap;
			}
			C[m][k] = _res;
			sd[k] = _res * pow(2, 0.5 * m);
		}
	}
}

int main(int argc, void* argv[])
{
		_calc_coefs(_phi_0, _input_signal, _scaling_func, step, 3);

	for(int k(0); k < 8; ++k)
	{
		std::cout << _input_signal(k * step + 0.5 * step) << "    " << sd[k] << std::endl;
	}

	system("pause");
	return 0;
}