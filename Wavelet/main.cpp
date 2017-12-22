#include "stdafx.h"
#include "mpi.h"
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include "Wavelet_transformh.h"
#include "omp.h"
using namespace Wavelet_transform;

int main(int argc, char* argv[]) {
	
	_m = 16;// atoi(argv[1]);
	k = pow(2, _m);
	N = k;
	_step_ = 1. / N;
	sd = new double[k];
	C = new double[k];

	sd_seq = new double[k];
	C_seq = new double[k];

	double _t1_seq, _t2_seq;
	double _t1_par, _t2_par;
	double t1, t2;
	MPI_Init(&argc, &argv);

	int ProcNum, ProcRank;

	MPI_Comm_size(MPI_COMM_WORLD, &ProcNum);
	MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank);

	if (ProcRank == 0)
	{
		_t1_seq = MPI_Wtime();
		_calc_coefs_seq(_phi_0, _input_signal, _scaling_func, _step_, _m);
		_t2_seq = MPI_Wtime();
		t1 = _t2_seq - _t1_seq;
	}
	_t1_par = MPI_Wtime();
	_calc_coefs(_phi_0, _input_signal, _scaling_func, _step_, _m);
	_t2_par = MPI_Wtime();
	t2 = _t2_par - _t1_par;

	if (ProcRank == 0)//comparation with err == 0.001;
	{
		std::cout << "is_C_equal: " << is_Ñ_equal(C, C_seq) << std::endl;
		std::cout << "is_sd_equal: " << is_sd_equal(sd, sd_seq) << std::endl;
	}

	if (ProcRank == 0)
	{
		std::cout << "seq time" << std::endl;
		std::cout << t1 << std::endl;
		//std::cout.flush();
		std::cout << "par time" << std::endl;
		std::cout << t2 << std::endl;

		//std::cout << "true val" << "   " << "par val" << "		" <<
		//	"seq val" << "		" << "par C[k]" << "	" << "seq_C[k]" << std::endl;
		//
		//for (int i(0); i < k; ++i)
		//{
		//	std::cout << _input_signal(i * _step_ + 0.5 * _step_)
		//		<< "    " << sd[i] << "		" << sd_seq[i] << "		"
		//		<< C[i] << "		" << C_seq[i] << std::endl;
		//	std::cout.flush();
		//}
		//std::cout << std::endl;
	}

	MPI_Finalize();
	delete[] C;
	delete[] sd;
	return 0;

}
