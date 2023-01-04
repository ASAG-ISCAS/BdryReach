/*
 * @Author: your name
 * @Date: 2021-11-11 17:24:27
 * @LastEditTime: 2021-12-28 18:18:30
 * @LastEditors: Please set LastEditors
 * @Description: 打开koroFileHeader查看配置 进行设置:
 * https://github.com/OBKoro1/koro1FileHeader/wiki/%E9%85%8D%E7%BD%AE
 * @FilePath: /Solver/test.cpp
 */

#include <capd/capdlib.h>
#include <ctime>
#include <zonotope/Zonotope.h>
#include <dynamics/LinearSys.h>
#include <dynamics/NonlinearSys.h>
#include <eigen3/Eigen/Core>
#include <functional>
#include <iostream>

#include <plotter/matplotlibcpp.h>
#include <plotter/plotter.h>

namespace plt = matplotlibcpp;

using namespace std;
using namespace reachSolver;
using namespace capd;
using capd::autodiff::Node;

double k  = 0.015;
double k2 = 0.01;
double g  = 9.81;

void _f(Node t, Node in[], int dimIn, Node out[], int dimOut, Node params[], int noParams)
{
  out[0] = params[0]*(in[1]-in[0]);
  out[1] = in[0]*(params[1]-in[2])-in[1];
  out[2] = in[1]*in[2]-params[2]*in[2];
}

int main()
{
	cout << "ha ha ha";
	int dimIn			   = 3;
	int dimOut			   = 3;
	int noParam			   = 3;
	//int MaxDerivativeOrder = 3;
	IMap f(_f, dimIn, dimOut, noParam);

	// NonlinearSys<double> mysys(f, 6, 1, 6);
	// ReachOptions<double> options;
	// Vector_t<double>	 center(6);
	// center << 2., 4., 4., 2., 10., 4.;
	// Matrix_t<double> generators(Eigen::MatrixXd::Identity(6, 6) * 0.2);
	// Zonotope<double> R0_(center, generators);

	// Vector_t<double> uc(1);
	// uc(0, 0) = 0.;
	// Matrix_t<double> ug(1, 1);
	// ug(0, 0) = 0.005;

	// Zonotope<double> U_(uc, ug);

	// options.set_time_step(1.);
	// options.set_taylor_terms(4);
	// options.set_zonotope_order(50);
	// options.set_alg("lin");
	// options.set_tensor_order(2);
	// options.set_tFinal(1);
	// options.set_tStart(0);
	// options.set_R0(R0_);
	// options.set_U(U_);
	// options.set_uTrans_nonlin(Eigen::MatrixXd::Zero(1, 1));
	// options.set_usekrylovError(1);
	// options.set_max_error(DBL_MAX * Eigen::MatrixXd::Ones(6, 1));

	// cout << "ha ha ha";

	// ReachableSet<double> R;

	// mysys.reach(options, R);

	// cout << R.Rtime_point().set()[0][0].rs() << endl;
	// //第10s，第0个ReachableSetElement的Zonotope
	// cout << "ha ha ha";
	
	// return 0;

}
