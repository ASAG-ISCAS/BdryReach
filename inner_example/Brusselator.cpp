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

double mu = 1.;

void _f(Node/* t*/, Node in[], int /*dimIn*/, Node out[], int/* dimOut*/, Node params[], int noParams){
    out[0] = 1+in[0]*in[0]*in[1]-1.5*in[0]-in[0];
    out[1] = 1.5*in[0]-in[0]*in[0]*in[1];
}

void _fBack(Node/* t*/, Node in[], int /*dimIn*/, Node out[], int/* dimOut*/, Node params[], int noParams){
    out[0] = -(1+in[0]*in[0]*in[1]-1.5*in[0]-in[0]);
    out[1] = -(1.5*in[0]-in[0]*in[0]*in[1]);
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
    center << 0.75, 0.7;
    Matrix_t<double> generators(2,2);
    generators<< 0.25,0,
                 0,0.3;

    Zonotope<double> R0_(center,generators);

    options.set_time_step(0.001);
    options.set_taylor_terms(4);
    options.set_zonotope_order(50);
    options.set_intermediate_order(50);
    options.set_error_order(20);
    options.set_alg("lin");
    options.set_tensor_order(3);

    options.set_tStart(0);

    options.set_R0(R0_);

    options.set_usekrylovError(1);
    options.set_max_error(DBL_MAX*Eigen::MatrixXd::Ones(2,1));
	
	plt::figure_size(1200, 780);
	
	clock_t start, end;
	start = clock();

	vector<Zonotope<double>> underR = UnderApprox::underReachClp(mysys, mysysBack, options, R0_, 0.55, 6, 0.02, 0.01, 0.11, 50);

	end = clock();
	cout << "time cost: " << (double) (end - start) / CLOCKS_PER_SEC << endl;
	
	for(int i = 1; i < underR.size() ; i++){
		Plotter::plotZonotope(underR[i], 1, 2, "g");
	}
	Plotter::plotZonotope(R0_, 1, 2, "k");
	plt::show();

}
