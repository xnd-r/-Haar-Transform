#pragma once
#include <cmath>

typedef double(*func)	(double t);
typedef double(*_func)	(func, double, int, long long int);

int _m;//scale
long long int k;//counter
//int N;// = pow(2, m);//amount of discret counts
double _step_;// = 1. / N;
double* C;//array of coefs of direct conversion
double* sd;//discret recovered signal
double* C_seq;
double* sd_seq;


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

double _scaling_func(func _f, double _t, int _m, long long int _k)
{
	double _tmp = pow(2, _m) * _t - _k;
	return pow(2, 0.5 * _m) * _f(_tmp);
}

void _calc_coefs(func _phi_0, func f, _func _f, double _step, int _m)
{
	int ProcNum, ProcRank;
	int number_steps = 128;
	long long int powerM = pow(2., _m);

	MPI_Comm_size(MPI_COMM_WORLD, &ProcNum);
	MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank);

	long long int segment = powerM / ProcNum;
	double _step_int = _step / number_steps;
	double _sdk = 0., _res = 0.;
	double a = (ProcRank*segment)* _step;
	long long int mod = powerM % ProcNum;
	int* displs = new int[ProcNum];
	int* recvcounts = new int[ProcNum];

	if (ProcRank == 0) {
		for (int i(0); i < ProcNum; ++i)
		{
			displs[i] = segment * i;
			recvcounts[i] = 1;
		}
	}

	for (long long int k = ProcRank*segment; k < (1 + ProcRank)*segment; ++k, a += _step)
	{

		for (int i = 0; i <= number_steps; i++)
			_res += (f(a + _step_int * (i + 0.5))
				* _f(_phi_0, a + _step_int * (i + 0.5), _m, k));

		_res *= _step_int;
		MPI_Gatherv(&_res, 1, MPI_DOUBLE, &C[k], recvcounts,
			displs, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		/*	MPI_Gatherv(&_sdk, 1, MPI_DOUBLE, &sd[k], recvcounts,
				displs, MPI_DOUBLE, 0, MPI_COMM_WORLD);*/
				//delete[] tmp;
	}
	if(mod != 0)
	{
		for (long long int k(mod*segment); k < powerM; ++k)
		{
			double a = k * _step;
			double _res = 0.;

			for (long long int i(0); i <= number_steps; i++)
			{
				_res += (f(a + _step_int * (i + 0.5))
					* _f(_phi_0, a + _step_int * (i + 0.5), _m, k));
			}
			_res *= _step_int;

			C[k] = _res;
		}
	}
	for (long long int k(0); k < powerM; ++k)
	{
		sd[k] = C[k] * pow(2, _m * 0.5);
	}
	delete[] displs;
	delete[] recvcounts;


}

void _calc_coefs_seq(func _phi_0, func f, _func _f, double _step, int _m) //C(m,k) and sd[k]
{
	int number_steps = 128;
	double _step_int = _step / number_steps;
	for (long long int k(0); k < pow(2., _m); ++k)
	{
		double a = k * _step;
		double _res = 0.;

		for (long long int i(0); i <= number_steps; i++)
		{
			_res += (f(a + _step_int * (i + 0.5))
				* _f(_phi_0, a + _step_int * (i + 0.5), _m, k));
		}
		_res *= _step_int;

		C_seq[k] = _res;
		sd_seq[k] = _res * pow(2, 0.5 * _m);
	}
}

bool is_sd_equal(double* _sd, double* _sd_seq)
{
	double _eps = 0.001;
	bool _flag = true;
	for(int k(0); k < pow(2, _m); ++k)
	{
		if (abs(_sd[k] - _sd_seq[k]) > _eps) _flag = false;
	}
	return _flag;
}

bool is_Ñ_equal(double* _C, double* _C_seq)
{
	double _eps = 0.001;
	bool _flag = true;
	for (int k(0); k < pow(2, _m); ++k)
	{
		if (abs(_C[k] - _C_seq[k]) > _eps) _flag = false;
	}
	return _flag;
}
