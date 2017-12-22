#pragma once

//	void _calc_coefs_seq(func _phi_0, func f, _func _f, double _step, int _m) //C(m,k) and sd[k]
//	{
//		int number_steps = 128;
//		double _step_int = _step / number_steps;
//		double _ans = 0.;
//		//for (int m(1); m < _m + 1; ++m)
//		//{
//			for (int k(0); k < pow(2., m); ++k)
//			{			
//				double a = k * _step;
////				double b = (k + 1) * _step;
//				double _res = 0.;
//
//				for (int i(0); i <= number_steps; i ++)
//				{
//					_res += (f(a + _step_int * (i + 0.5))
//						  * _f(_phi_0, a + _step_int * (i + 0.5), m, k));
//				}
//				_res *= _step_int;
//			
//				C[m][k] = _res;
//				sd[k] = _res * pow(2, 0.5 * m);
//			}
////		}
//	}
//
