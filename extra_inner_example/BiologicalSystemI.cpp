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


void _f(Node/* t*/, Node in[], int /*dimIn*/, Node out[], int/* dimOut*/, Node params[], int noParams){
    out[0] = -(0.4*in[0] - 5*in[2]*in[3]);
    out[1] = 0.4*in[0] - in[1];
	out[2] = in[1] - 5*in[2]*in[3];
	out[3] = 5*in[4]*in[5] - 5*in[2]*in[3];
	out[4] = -(5*in[4]*in[5] - 5*in[2]*in[3]);
	out[5] = 0.5*in[6] - 5*in[4]*in[5];
	out[6] = -(0.5*in[6] - 5*in[4]*in[5]);

}

void _fBack(Node/* t*/, Node in[], int /*dimIn*/, Node out[], int/* dimOut*/, Node params[], int noParams){
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
IMap fBack(_fBack, dimIn,dimOut,noParam,MaxDerivativeOrder);
int main()
{
    NonlinearSys<double> mysys(f, 7, 0, 7);
    NonlinearSys<double> mysysBack(fBack, 7, 0, 7);

    ReachOptions<double> options;

    //create R0
    Vector_t<double> center(7);
    center << 0.10875, 0.10875, 0.10875, 0.10875, 0.10875, 0.10875, 0.10875;
    Matrix_t<double> generators(Eigen::MatrixXd::Identity(7,7) * 0.00875);
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
    options.set_max_error(DBL_MAX*Eigen::MatrixXd::Ones(7,1));

	plt::figure_size(1200, 780);

    clock_t start, end;
	start = clock();
	
    //0.2 final time 1
	// vector<Zonotope<double>> underR = UnderApprox::underReachClp(mysys, mysysBack, options, R0_, 0.2, 1, 1, 0.01,0.1, 2, 2);

    //1.3 
    // vector<Zonotope<double>> underR = UnderApprox::underReachClp(mysys, mysysBack, options, R0_, 1.3, 1, 0.005, 0.01,0.05, 5, 5);

    //1.3 final time 2 this good 466s 0.41
    //  vector<Zonotope<double>> underR = UnderApprox::underReachClp(mysys, mysysBack, options, R0_, 1.3, 1, 0.005, 0.01,0.1, 5, 5);

    //  vector<Zonotope<double>> underR = UnderApprox::underReachClp(mysys, mysysBack, options, R0_, 1.3, 1, 0.005, 0.01,0.26, 5, 5);

    vector<Zonotope<double>> underR = UnderApprox::underReachClp(mysys, mysysBack, options, R0_, 0.65, 2, 1, 0.005,0.0025, 2, 10);

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
