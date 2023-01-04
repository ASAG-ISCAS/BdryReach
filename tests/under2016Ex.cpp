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


// void _f(Node/* t*/, Node in[], int /*dimIn*/, Node out[], int/* dimOut*/, Node params[], int noParams){
//     out[0] = -in[1];
//     out[1] = -(0.2 - 0.7*sin(in[0]) - 0.05 * in[1]);
// }

// void _fBack(Node/* t*/, Node in[], int /*dimIn*/, Node out[], int/* dimOut*/, Node params[], int noParams){
//     out[0] = in[1];
//     out[1] = (0.2 - 0.7*sin(in[0]) - 0.05 * in[1]);
// }

void _f(Node/* t*/, Node in[], int /*dimIn*/, Node out[], int/* dimOut*/, Node params[], int noParams){
    out[0] = in[1];
    out[1] = (0.2 - 0.7*sin(in[0]) - 0.05 * in[1]);
}

void _fBack(Node/* t*/, Node in[], int /*dimIn*/, Node out[], int/* dimOut*/, Node params[], int noParams){
    out[0] = -in[1];
    out[1] = -(0.2 - 0.7*sin(in[0]) - 0.05 * in[1]);
}

//be care for
int dimIn = 3;

int dimOut = 2;
int noParam = 0;
int MaxDerivativeOrder = 3;

IMap f(_f, dimIn,dimOut,noParam,MaxDerivativeOrder);
IMap fBack(_fBack, dimIn,dimOut,noParam,MaxDerivativeOrder);
int main()
{
    NonlinearSys<double> mysys(f, 2, 0, 2);
	NonlinearSys<double> mysysBack(fBack, 2, 0, 2);

    ReachOptions<double> options;

    // create R0
    Vector_t<double> center(2);
    center << 0, 3;
	// center << -7.5, 2.8;
    Matrix_t<double> generators(2,2);
    generators<< 0.1,0,
                 0,0.1;
    // generators<< 0.07,0.07,0.1,
    //              0.07,-0.07,0;

//     Matrix_t<double> generators(2,100);
// 	center << -7.9721, 2.4942;
// generators <<
// 0.012647, 3.4143e-05, -0.11387, -2.4092e-06, -0.00010345, -2.1572e-06, -0.00010242, -1.9109e-06, -0.00010131, -0.00010012, -1.6705e-06, -0.0001004, -9.9248e-05, -9.8021e-05, -9.6727e-05, -9.5368e-05, -9.395e-05, -9.2475e-05, -9.0949e-05, -8.9375e-05, -8.7756e-05, -8.6097e-05, -8.4401e-05, -8.2672e-05, -8.0914e-05, -7.9128e-05, -7.732e-05, -7.549e-05, -7.3645e-05, -7.1785e-05, -6.9912e-05, -6.8033e-05, -6.615e-05, -6.4263e-05, -6.2379e-05, -6.0493e-05, -5.8604e-05, -5.6726e-05, -5.4856e-05, -0.065113, -4.8403e-05, -3.0648e-05, -5.1098e-05, -3.0277e-05, -5.3629e-05, -2.9892e-05, -5.5992e-05, -2.9492e-05, -8.7261e-05, -8.885e-05, -9.0249e-05, -9.1453e-05, -9.2446e-05, -9.3227e-05, -9.3798e-05, -9.4162e-05, -9.4316e-05, -9.4265e-05, -9.4012e-05, -9.3561e-05, -9.2914e-05, -9.2075e-05, -9.1048e-05, 1.8123e-06, -8.9838e-05, 1.8505e-06, -8.8454e-05, 2.4889e-06, -8.6897e-05, -8.5174e-05, 3.623e-06, -8.3288e-05, -8.1248e-05, -7.9057e-05, -7.6726e-05, -7.4258e-05, 9.0357e-06, -7.1663e-05, -6.8944e-05, -6.6106e-05, -6.3165e-05, -6.0124e-05, -5.699e-05, -5.3769e-05, -5.0469e-05, -4.7097e-05, -4.3661e-05, -4.0166e-05, -3.6622e-05, -3.3035e-05, -2.9414e-05, -2.5764e-05, -2.2093e-05, -1.8407e-05, -1.4741e-05, -1.1143e-05, -7.482e-06, -3.7658e-06, 5.7119e-05, 0,
// 0.0036926, 6.037e-06, 0.0097159, 3.5558e-07, 1.5268e-05, 3.2184e-07, 1.528e-05, 2.8819e-07, 1.5279e-05, 1.5264e-05, 2.5467e-07, 1.5474e-05, 1.5464e-05, 1.544e-05, 1.5404e-05, 1.5355e-05, 1.5294e-05, 1.522e-05, 1.5136e-05, 1.5039e-05, 1.4932e-05, 1.4814e-05, 1.4685e-05, 1.4546e-05, 1.4398e-05, 1.424e-05, 1.4072e-05, 1.3896e-05, 1.3711e-05, 1.3518e-05, 1.3317e-05, 1.3108e-05, 1.2892e-05, 1.267e-05, 1.2441e-05, 1.2206e-05, 1.1963e-05, 1.1716e-05, 1.1463e-05, 0.014146, 9.8955e-05, 6.2656e-05, 0.00010669, 6.3219e-05, 0.00011441, 6.3772e-05, 0.0001221, 6.4315e-05, 0.0001946, 0.00020272, 0.00021077, 0.00021873, 0.00022656, 0.00023424, 0.00024177, 0.00024914, 0.00025634, 0.00026335, 0.00027018, 0.00027681, 0.00028324, 0.00028946, 0.00029547, 1.0001e-07, 0.00030126, 9.3988e-08, 0.00030684, 1.159e-07, 0.00031219, 0.00031731, 1.3996e-07, 0.00032219, 0.00032684, 0.00033124, 0.0003354, 0.00033932, 1.9733e-07, 0.000343, 0.00034642, 0.00034958, 0.00035251, 0.0003552, 0.00035765, 0.00035985, 0.0003618, 0.00036352, 0.00036498, 0.0003662, 0.00036718, 0.00036793, 0.00036845, 0.00036873, 0.00036879, 0.00036862, 0.00036891, 0.0003717, 0.00037429, 0.00037668, 0, 0.00037886;
// ;
// 	Vector_t<double> center(2);
//     Matrix_t<double> generators(2,44);

// center <<9.9214, 2.748;
// generators <<
// 0.033115, 4.0883e-06, 1.3122e-06, 1.5048e-06, 3.8311e-06, 1.2475e-06, 1.6461e-06, 4.5574e-06, 3.4981e-06, 9.3004e-07, 3.7638e-06, 2.2547e-06, 4.2061e-06, 3.5626e-06, 8.8482e-07, 7.0209e-06, 9.2125e-07, 9.571e-07, 7.1689e-07, 9.0508e-07, 1.0406e-06, 9.3603e-07, 1.0633e-06, 9.6611e-07, 1.0853e-06, 9.9527e-07, 1.1065e-06, 1.0235e-06, 1.127e-06, 1.0506e-06, 1.1467e-06, 1.0767e-06, 1.1656e-06, 1.1017e-06, 1.1836e-06, 9.4251e-06, 1.1256e-06, 1.2007e-06, 6.7501e-06, 6.5725e-06, 9.4506e-06, 9.7324e-06, 1.6901e-05, 0.25327,
// 0.0040409, 6.527e-07, 2.2658e-07, 2.695e-07, 6.995e-07, 2.3133e-07, 3.1557e-07, 9.0196e-07, 7.1371e-07, 1.9536e-07, 8.2124e-07, 5.0677e-07, 9.7221e-07, 8.4557e-07, 2.1125e-07, 1.7086e-06, 2.2446e-07, 2.377e-07, 1.8093e-07, 2.2887e-07, 2.6686e-07, 2.4075e-07, 2.7681e-07, 2.5248e-07, 2.8658e-07, 2.6403e-07, 2.9616e-07, 2.7535e-07, 3.0551e-07, 2.8641e-07, 3.146e-07, 2.9716e-07, 3.2342e-07, 3.0758e-07, 3.3193e-07, 2.6435e-06, 3.1762e-07, 3.4011e-07, 2.022e-06, 2.0165e-06, 2.9626e-06, 3.1374e-06, 5.5753e-06, 0.10566;


    Zonotope<double> R0_(center,generators);

    options.set_time_step(0.2);
    options.set_taylor_terms(4);
    options.set_zonotope_order(50);
    options.set_intermediate_order(50);
    options.set_error_order(20);
    options.set_alg("lin");
    options.set_tensor_order(3);

	options.set_tFinal(3);
    options.set_tStart(0);

    options.set_R0(R0_);

    options.set_usekrylovError(1);
    options.set_max_error(DBL_MAX*Eigen::MatrixXd::Ones(2,1));



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
	
	clock_t start, end;
	start = clock();

	// vector<Zonotope<double>> underR = UnderApprox::underReachClp(mysys, mysysBack, options, R0_, 0.5, 8, 0.01, 0.01, 0.01, 50);

	vector<Zonotope<double>> underR = UnderApprox::underReachClp(mysys, mysysBack, options, R0_, 0.5, 6, 0.01, 0.01, 0.1, 50);
	end = clock();
	cout << "time cost: " << (double) (end - start) / CLOCKS_PER_SEC << endl;
	
	// vector<Zonotope<double>> underR = UnderApprox::underReachPlus(mysys, options, R0_, 2, 1, 1);

	for(int i = 1; i < underR.size() ; i++){
		Plotter::plotZonotope(underR[i], 1, 2, "g");
	}
	// Plotter::plotZonotope(newZ, 1, 2, "g");
	Plotter::plotZonotope(R0_, 1, 2, "k");
	plt::show();

	cout << underR[3] << endl; 

	// Vector_t<double> c1(2), c2(2);

	// // c1 << 0, 1;
	// // c2 << 0,1.01;

	// 	c1 << 0.1, 1;
	// c2 << 0,1;
	// // Matrix_t<double> G1(2,5),G2(2,6);
	// Matrix_t<double> G1(2,2),G2(2,2);

	// // G1 << 1,0,0,1,1,
	// // 0,-1,0,-1,-3;

	// // G1<<1,0,
	// // 0,1;

	// // G2 <<2,0,
	// // 1,1;
	// // G2 << 1, 0, 1, 1, 1, 2,
	// // 0, 1, 1, -1, 3, -2;

	// // G2 << 1, 0, 1, 1, 1,
	// // 0, 1, 1, -1, 3;

	// G1 << 1,0,
	// 0,0.4;
	// G2 << 2,0,
	// 1,1;

	// Zonotope<double> Z1(c1,G1), Z2(c2,G2);

	// // cout << Z2 << endl;

	// cout << UnderApprox::zonotopeIsIn(Z1, Z2) << endl;


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
