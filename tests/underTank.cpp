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
#include <underApprox/underApprox.h>

#include <glpk.h>
namespace plt = matplotlibcpp;

using namespace std;
using namespace reachSolver;
using namespace capd;
using capd::autodiff::Node;

double k  = 0.015;
double k2 = 0.01;
double g  = 9.81;

void _f(Node/* t*/, Node in[], int /*dimIn*/, Node out[], int/* dimOut*/, Node params[], int noParams){
	out[0] = 0.1 + k2 * (4 - in[5]) - k * sqrt(2 * g) * sqrt(in[0]);
	out[1] = k * sqrt(2 * g) * (sqrt(in[0]) - sqrt(in[1]));
	out[2] = k * sqrt(2 * g) * (sqrt(in[1]) - sqrt(in[2]));
	out[3] = k * sqrt(2 * g) * (sqrt(in[2]) - sqrt(in[3]));
	out[4] = k * sqrt(2 * g) * (sqrt(in[3]) - sqrt(in[4]));
	out[5] = k * sqrt(2 * g) * (sqrt(in[4]) - sqrt(in[5]));
}

void _fBack(Node/* t*/, Node in[], int /*dimIn*/, Node out[], int/* dimOut*/, Node params[], int noParams){
	out[0] = -(0.1 + k2 * (4 - in[5]) - k * sqrt(2 * g) * sqrt(in[0]));
	out[1] = -k * sqrt(2 * g) * (sqrt(in[0]) - sqrt(in[1]));
	out[2] = -k * sqrt(2 * g) * (sqrt(in[1]) - sqrt(in[2]));
	out[3] = -k * sqrt(2 * g) * (sqrt(in[2]) - sqrt(in[3]));
	out[4] = -k * sqrt(2 * g) * (sqrt(in[3]) - sqrt(in[4]));
	out[5] = -k * sqrt(2 * g) * (sqrt(in[4]) - sqrt(in[5]));
}


//be care for
int dimIn = 7;

int dimOut = 6;
int noParam = 0;
int MaxDerivativeOrder = 3;

IMap f(_f, dimIn,dimOut,noParam,MaxDerivativeOrder);
IMap fBack(_fBack, dimIn,dimOut,noParam,MaxDerivativeOrder);
int main()
{
    NonlinearSys<double> mysys(f, 6, 0, 6);
	NonlinearSys<double> mysysBack(fBack, 6, 0, 6);

    ReachOptions<double> options;

    //create R0
	Vector_t<double>	 center(6);
	center << 2., 4., 4., 2., 10., 4.;
	Matrix_t<double> generators(Eigen::MatrixXd::Identity(6, 6) * 0.2);
	Zonotope<double> R0_(center, generators);

    options.set_taylor_terms(4);
    options.set_zonotope_order(50);
    options.set_intermediate_order(50);
    options.set_error_order(20);
    options.set_alg("lin");
    options.set_tensor_order(3);

	options.set_tFinal(6);
    options.set_tStart(0);

    options.set_R0(R0_);

    options.set_usekrylovError(1);
    options.set_max_error(DBL_MAX*Eigen::MatrixXd::Ones(6,1));
	
	
	plt::figure_size(1200, 780);

	// vector<Zonotope<double>> underR = UnderApprox::underReach(mysys, options, R0_, 0.5, 6, 0.001);
	
	clock_t start, end;
	start = clock();

	vector<Zonotope<double>> underR = UnderApprox::underReachClp(mysys, mysysBack, options, R0_, 20, 4, 100, 1, 1, 3);

	end = clock();
	cout << "time cost: " << (double) (end - start) / CLOCKS_PER_SEC << endl;

	// vector<Zonotope<double>> underR = UnderApprox::underReachPlus(mysys, options, R0_, 2, 1, 1);

	for(int i = 1; i < underR.size(); i++){
		Plotter::plotZonotope(underR[i], 1, 2, "g");
	}
	// Plotter::plotZonotope(newZ, 1, 2, "g");
	Plotter::plotZonotope(R0_, 1, 2, "k");
	plt::show();
}
