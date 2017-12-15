#include <cmath>
namespace Wavelet 
{
	typedef double(*func)	(double t);
	typedef double(*_func)	(func, double, int, int);
	typedef double(*coef)	(func, _func, double, double, int, int);

	int m;//scale
	int k;//counter
	int N;//amount of discret counts
	double step = 1. / N;
	double C[4][16];//array of coefs of direct conversion
	double sd[16];//discret recovered signal

	inline double _input_signal(double _t) // continious signal for transformation
	{
		double _res;
		(_t >= 0 && _t <= 1) ? _res = abs(sin(10 * _t) + 5 * _t) : _res = -1.;
		return _res;
	}

	inline int _phi_0(double _t) //Haar scaling
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

	void _calc_coefs(func _phi_0, func f, _func _f, double _t, double _step, int _m) //C(m,k) and sd[k]
	{
		double underint;
		double _step_int = _step / 128;

		for (int m(0); m < _m; ++m)
		{
			for (int k(0); k < pow(2, m); ++k)
			{
				double _res = 0.;
				for(int i(0); i < 128; ++i)
				{
					underint = f(_t) * _f(_phi_0, _t, m, k);
				}
				_res *= _step_int;
				C[m][k] = _res;
				sd[k] = _res * pow(2, 0.5 * m);
			}
		}
	}

}