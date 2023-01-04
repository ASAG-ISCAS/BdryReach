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

double a = 0.2;
double b = 0.2;
double c = 5.7;
void _f(Node/* t*/, Node in[], int /*dimIn*/, Node out[], int/* dimOut*/, Node params[], int noParams){
    out[0] = -in[1] - in[2];
    out[1] = in[0] + a * in[1];
    out[2] = b + in[2]*(in[0] - c);

}

void _fBack(Node/* t*/, Node in[], int /*dimIn*/, Node out[], int/* dimOut*/, Node params[], int noParams){
    out[0] = -(-in[1] - in[2]);
    out[1] = -(in[0] + a * in[1]);
    out[2] = -(b + in[2]*(in[0] - c));


}

//be care for
int dimIn = 4;

int dimOut = 3;
int noParam = 0;
int MaxDerivativeOrder = 3;

IMap f(_f, dimIn,dimOut,noParam,MaxDerivativeOrder);
IMap fBack(_fBack, dimIn,dimOut,noParam,MaxDerivativeOrder);
int main()
{
    NonlinearSys<double> mysys(f, 3, 0, 3);
    NonlinearSys<double> mysysBack(fBack, 3, 0, 3);

    ReachOptions<double> options;

    //create R0
    Vector_t<double> center(3);
    center << 0.05, -8.35, 0.05;
    Matrix_t<double> generators(Eigen::MatrixXd::Identity(3,3) * 0.15);
    Zonotope<double> R0_(center,generators);

    options.set_taylor_terms(4);
    options.set_zonotope_order(50);
    options.set_intermediate_order(50);
    options.set_error_order(20);
    options.set_alg("lin");
    options.set_tensor_order(3);

    options.set_tStart(0);

    options.set_R0(R0_);

    options.set_usekrylovError(1);
    options.set_max_error(DBL_MAX*Eigen::MatrixXd::Ones(3,1));

	plt::figure_size(1200, 780);

    clock_t start, end;
	start = clock();
	
   
	// vector<Zonotope<double>> underR = UnderApprox::underReachClp(mysys, mysysBack, options, R0_, 1.5, 1, 0.05, 0.01,0.05, 50, 20);
    // 目前最优参数
     vector<Zonotope<double>> underR = UnderApprox::underReachClp(mysys, mysysBack, options, R0_, 1.5, 1, 0.04, 0.01,0.015, 50, 20);
    end = clock();
	cout << "time cost: " << (double) (end - start) / CLOCKS_PER_SEC << endl;

    for(int i = 1; i < underR.size() ; i++){
		cout << "BZF{" << i << "} = zonotope(" << underR[i] << ");" << endl;
	}

	for(int i = 1; i < underR.size(); i++){
		Plotter::plotZonotope(underR[i], 1, 2, "g");
	}
	plt::show();
}
