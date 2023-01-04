/*
 * @Author: your name
 * @Date: 2021-11-11 17:24:27
 * @LastEditTime: 2021-12-30 22:41:19
 * @LastEditors: Please set LastEditors
 * @Description: 打开koroFileHeader查看配置 进行设置: https://github.com/OBKoro1/koro1FileHeader/wiki/%E9%85%8D%E7%BD%AE
 * @FilePath: /Solver/test.cpp
 */

#include <overApprox/overApprox.h>
#include <dynamics/NonlinearSys.h>
#include <dynamics/LinearSys.h>

#include <zonotope/Zonotope.h>
#include <zonotope/polyZonotope.h>
#include <zonotope/globalfunctions.h>

#include <iostream>
#include <eigen3/Eigen/Core>
#include <ctime>

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

void _f(Node/* t*/, Node in[], int /*dimIn*/, Node out[], int/* dimOut*/, Node params[], int noParams){
	out[0] = in[6] + 0.1 + k2 * (4 - in[5]) - k * sqrt(2 * g) * sqrt(in[0]);
	out[1] = k * sqrt(2 * g) * (sqrt(in[0]) - sqrt(in[1]));
	out[2] = k * sqrt(2 * g) * (sqrt(in[1]) - sqrt(in[2]));
	out[3] = k * sqrt(2 * g) * (sqrt(in[2]) - sqrt(in[3]));
	out[4] = k * sqrt(2 * g) * (sqrt(in[3]) - sqrt(in[4]));
	out[5] = k * sqrt(2 * g) * (sqrt(in[4]) - sqrt(in[5]));
}

//be care for
int dimIn			   = 7;
int dimOut			   = 6;
int noParam			   = 0;
int MaxDerivativeOrder = 3;

IMap f(_f, dimIn,dimOut,noParam,MaxDerivativeOrder);

int main(){
    double step = 1;
    double tFinal = 400;
    int N = tFinal/step;
    int dim1 = 1;
    int dim2 = 2;


    NonlinearSys<double> mysys(f, 6, 1, 6);
    ReachOptions<double> options;

    //create R0
	Vector_t<double>	 center(6);
	center << 2., 4., 4., 2., 10., 4.;
	Matrix_t<double> generators(Eigen::MatrixXd::Identity(6, 6) * 0.2);
	Zonotope<double> R0_(center, generators);

	Vector_t<double> uc(1);
	uc(0, 0) = 0.;
	Matrix_t<double> ug(1, 1);
	ug(0, 0) = 0.005;

    Zonotope<double> U_(uc, ug);

    options.set_time_step(step);
    options.set_taylor_terms(4);
    options.set_zonotope_order(50);
    options.set_intermediate_order(50);
    options.set_error_order(20);
    options.set_alg("lin");
    options.set_tensor_order(3);

    options.set_tFinal(tFinal);
    options.set_tStart(0);

    options.set_R0(R0_);
    options.set_U(U_);

    options.set_uTrans_nonlin(Eigen::MatrixXd::Zero(1, 1));
    options.set_usekrylovError(1);
    options.set_max_error(DBL_MAX * Eigen::MatrixXd::Ones(6, 1));
    ReachableSet<double> R;

    vector<Zonotope<double>> R0_border = Zono_border(R0_);
    
    vector<ReachableSet<double>> BdReachset = OverApprox::BdReach(mysys, options, R0_);

    mysys.reach(options, R);

    plt::figure_size(1200, 780);
    Plotter::plotReach(R, 1, 2, "r");

    for(int i = 0; i < BdReachset.size(); i++){

        Plotter::plotReach(BdReachset[i], 1, 2, "b");
    }
    plt::show();
                
}

