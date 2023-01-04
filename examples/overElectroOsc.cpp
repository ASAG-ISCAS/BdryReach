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


void _f(Node/* t*/, Node in[], int /*dimIn*/, Node out[], int/* dimOut*/, Node params[], int noParams){
    out[0] = in[1];
    out[1] = (0.2 - 0.7*sin(in[0]) - 0.05 * in[1]);
}

//be care for
int dimIn			   = 3;
int dimOut			   = 2;
int noParam			   = 0;
int MaxDerivativeOrder = 3;

IMap f(_f, dimIn,dimOut,noParam,MaxDerivativeOrder);

int main(){
    double step = 0.01;
    double tFinal = 6;
    int N = tFinal/step;
    int dim1 = 1;
    int dim2 = 2;


    NonlinearSys<double> mysys(f, 2, 0, 2);

    ReachOptions<double> options;

    //create R0
	Vector_t<double>	 center(2);
    Matrix_t<double> generators(2,2);
	center << 0,3;
    generators<< 0.1,0,
        0,0.1;

	Zonotope<double> R0_(center, generators);

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

    options.set_uTrans_nonlin(Eigen::MatrixXd::Zero(1, 1));
    options.set_usekrylovError(1);
    options.set_max_error(DBL_MAX * Eigen::MatrixXd::Ones(2, 1));

    ReachableSet<double> R;
    
    vector<ReachableSet<double>> BdReachset = OverApprox::BdReach(mysys, options, R0_);
    
    mysys.reach(options, R);

    plt::figure_size(1200, 780);
    Plotter::plotReach(R, 1, 2, "r");

    for(int i = 0; i < BdReachset.size(); i++){

        Plotter::plotReach(BdReachset[i], 1, 2, "b");
    }
    plt::show();
}

