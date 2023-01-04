/*
 * @Author: your name
 * @Date: 2021-11-11 17:24:27
 * @LastEditTime: 2021-12-28 18:18:30
 * @LastEditors: Please set LastEditors
 * @Description: 打开koroFileHeader查看配置 进行设置:
 * https://github.com/OBKoro1/koro1FileHeader/wiki/%E9%85%8D%E7%BD%AE
 * @FilePath: /Solver/test.cpp
 */

#include <plotter/matplotlibcpp.h>
#include <plotter/plotter.h>
#include <underApprox/underApprox.h>


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

    options.set_tStart(0);

    options.set_R0(R0_);

    options.set_usekrylovError(1);
    options.set_max_error(DBL_MAX*Eigen::MatrixXd::Ones(6,1));
	
	
	plt::figure_size(1200, 780);
	
	clock_t start, end;
	start = clock();

	// vector<Zonotope<double>> underR = UnderApprox::underReachClp(mysys, mysysBack, options, R0_, 20, 4, 100, 1, 1, 3, 3);

	//160s 171s
	// vector<Zonotope<double>> underR = UnderApprox::underReachClp(mysys, mysysBack, options, R0_, 20, 8, 100, 1, 1, 3, 50);

// 	[1.862966380918353 -0.001473058126840164 -0.01254175447648878 -0.007609894906006537 -0.01019663255244353 -0.004787715770853222 -0.003199724763632948 ;
// 1.807867417364385 0.009882707389619161 -0.01243258816786255 -0.006145740607046082 -0.01700157243742399 -0.003132829097064877 0.000220783697146118 ;
// 2.038046729812934 0.02614825654396373 -0.01081257114496257 -0.004452993215694236 -0.02178721211247605 0.0007799937239454394 0.007588616443787138 ;
// 2.462516704623012 0.03456682973276262 -0.007334216940736357 0.0003068709889923792 -0.01967780945567174 0.008070633717943271 0.01491544666381639 ;
// 3.278776888111691 0.0306406060789433 0.006485087160192328 0.01153344132611478 -0.01320946606412981 0.01553231269889925 0.01731566265429951 ;
// 4.366131911502078 0.01828778548090716 0.02548068416674968 0.01961044432389498 0.008249996460672985 0.01570792562457117 0.01270617048009444 ;]

	//175s
	// vector<Zonotope<double>> underR = UnderApprox::underReachClp(mysys, mysysBack, options, R0_, 25, 7, 100, 1, 1, 3, 50);

	//120s 221s
	// vector<Zonotope<double>> underR = UnderApprox::underReachClp(mysys, mysysBack, options, R0_, 40, 3, 100, 0.5, 0.25, 10, 50);
	
	//还不如 上面的 175s 
	// vector<Zonotope<double>> underR = UnderApprox::underReachClp(mysys, mysysBack, options, R0_, 25, 7, 100, 0.25, 0.5, 5, 50);

	// vector<Zonotope<double>> underR = UnderApprox::underReachClp(mysys, mysysBack, options, R0_, 25, 7, 100, 0.25, 0.5, 5, 50);

	//80s final time 1
	// vector<Zonotope<double>> underR = UnderApprox::underReachClp(mysys, mysysBack, options, R0_, 40, 2, 100, 0.5, 1, 10, 10);

	// vector<Zonotope<double>> underR = UnderApprox::underReachClp(mysys, mysysBack, options, R0_, 40, 3, 100, 0.5, 0.25, 10, 50);

	vector<Zonotope<double>> underR = UnderApprox::underReachClp(mysys, mysysBack, options, R0_, 20, 6, 100, 1, 1, 5, 20);
	end = clock();
	cout << "time cost: " << (double) (end - start) / CLOCKS_PER_SEC << endl;

	for(int i = 1; i < underR.size(); i++){
		Plotter::plotZonotope(underR[i], 1, 2, "g");
	}
	// Plotter::plotZonotope(newZ, 1, 2, "g");
	Plotter::plotZonotope(R0_, 1, 2, "k");
	plt::show();
}
