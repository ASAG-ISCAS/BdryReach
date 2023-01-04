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


void _f(Node/* t*/, Node in[], int /*dimIn*/, Node out[], int/* dimOut*/, Node params[], int noParams){
    out[0] = 0.4*in[0] - 5*in[2]*in[3];
    out[1] = -0.4*in[0] + in[1];
	out[2] = -in[1] + 5*in[2]*in[3];
	out[3] = -5*in[4]*in[5] + 5*in[2]*in[3];
	out[4] = 5*in[4]*in[5] - 5*in[2]*in[3];
	out[5] = -0.5*in[6] + 5*in[4]*in[5];
	out[6] = 0.5*in[6] - 5*in[4]*in[5];

}

//be care for
int dimIn = 8;

int dimOut = 7;
int noParam = 0;
int MaxDerivativeOrder = 3;

IMap f(_f, dimIn,dimOut,noParam,MaxDerivativeOrder);

int main()
{
    NonlinearSys<double> mysys(f, 7, 0, 7);
    ReachOptions<double> options;

    //create R0
    Vector_t<double> center(7);
    center << -0.007, -0.007, -0.007, -0.007, -0.007, -0.007, -0.007;
    Matrix_t<double> generators(Eigen::MatrixXd::Identity(7,7) * 0.008);
    Zonotope<double> R0_(center,generators);

	// cout << "good" << endl;

    options.set_time_step(0.2);
    options.set_taylor_terms(4);
    options.set_zonotope_order(50);
    options.set_intermediate_order(50);
    options.set_error_order(20);
    options.set_alg("lin");
    options.set_tensor_order(3);

	options.set_tFinal(0.2);
    options.set_tStart(0);

    options.set_R0(R0_);

    options.set_usekrylovError(1);
    options.set_max_error(DBL_MAX*Eigen::MatrixXd::Ones(7,1));



	// ReachableSet<double> R;

	// clock_t start, end;
	// start = clock();
	// mysys.reach(options, R);
	// end = clock();
	// cout << "time cost: " << (double) (end - start) / CLOCKS_PER_SEC << endl;

	// vector<Zonotope<double>> R0_border = Zono_border(R0_);
    // vector<ReachableSet<double>> R_border(R0_border.size());
	// // cout << R0_border << endl;

    // for(int i = 0; i < R0_border.size(); i++){
    //     options.set_R0(R0_border[i]);
    //     mysys.reach(options, R_border[i]);
    // }

	// //cout << R.Rtime_point().set()[10][0].rs() << endl;
	// // plt::figure_size(1200, 780);
	// // for(auto i = 0; i < 1; i++){
	// // 	//cout << R.Rtime_point().set()[i][0].rs() << endl;
	// // 	Plotter::plotZonotope(R.Rtime_point().set()[i][0].rs(), 1, 2, "b");
	// // 	for(int j = 0; j < R0_border.size(); j++){
	// // 		Plotter::plotZonotope(R_border[j].Rtime_point().set()[i][0].rs(), 1, 2, "r");
	// // 	}
	// // }
	// // plt::show();
	// // cout << R.Rtime_point().set()[0][0].rs() << endl;

	// // Plotter::plotReach(R, 1, 2, "b");
	// // Plotter::plotOvertime(R, 1 , "b");
	// //第10s，第0个ReachableSetElement的Zonotope
	

	options.set_R0(R0_);
	// Zonotope<double> newZ = UnderApprox::getUnder(mysys, options, R0_);
	plt::figure_size(1200, 780);

	// vector<Zonotope<double>> underR = UnderApprox::underReach(mysys, options, R0_, 0.5, 6, 0.001);
	
	vector<Zonotope<double>> underR = UnderApprox::underReachClp(mysys, mysys, options, R0_, 0.2, 1, 0.016, 0.01,0.1);

	// vector<Zonotope<double>> underR = UnderApprox::underReachPlus(mysys, options, R0_, 2, 1, 1);

	cout << underR.size() << endl;
	for(int i = 1; i < underR.size(); i++){
		Plotter::plotZonotope(underR[i], 1, 2, "g");
	}
	// Plotter::plotZonotope(newZ, 1, 2, "g");
	plt::show();


	// Vector_t<double> c1(2),c2(2);
	// Matrix_t<double> G1(2,2),G2(2,2);

	// c1 << 0,0;
	// G1 << 1,0,
	// 	0,1;
	// c2 << 3,0;
	// G2 << 0.7, -0.7,
	// 	0.7,0.7;
	
	// Zonotope<double> Z1(c1,G1),Z2(c2,G2);

	// vector<double> al(2,-1),au(2,1);

	// UnderApprox::getScale(Z1,Z2,al,au);

	// cout << al << endl;
	// cout << au << endl;

	// Vector_t<double> Cadd(al.size());
	// Vector_t<double> Gsacle(al.size());

	// for(int i = 0; i < al.size(); i++){

	// 	Cadd(i) = 0.5 * (au[i] + al[i]);

	// 	Gsacle(i) = 0.5 * (au[i] - al[i]);
	// }

	// Zonotope<double> res(Z2.center() + Z2.generators() * Cadd, Z2.generators() * Gsacle.asDiagonal());

	// Plotter::plotZonotope(Z1,1,2,"r");
	// Plotter::plotZonotope(Z2,1,2,"b");

	// Plotter::plotZonotope(res,1,2,"y");
	// plt::show();
}
