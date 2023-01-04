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
    out[0] = 3*in[2] - in[0]*in[5] ; 
    out[1] = in[3] - in[1]*in[5] ; 
    out[2] = in[0]*in[5] - 3*in[2]; 
    out[3] = in[1]*in[5] - in[3]; 
    out[4] = 3*in[2] + 5*in[0] - in[4] ; 
    out[5] = 5*in[4] + 3*in[2] + in[3] - in[5] * (in[0] + in[1] + 2*in[7] + 1) ; 
    out[6] = 5*in[3] + in[1] - 0.5*in[6] ; 
    out[7] = 5*in[6] - 2*in[5]*in[7] + in[8] - 0.2*in[7] ; 
    out[8] = 2*in[5]*in[7] - in[8];

}

void _fBack(Node/* t*/, Node in[], int /*dimIn*/, Node out[], int/* dimOut*/, Node params[], int noParams){
    out[0] = -(3*in[2] - in[0]*in[5] )   ;
    out[1] = -(in[3] - in[1]*in[5] )   ;
    out[2] = -(in[0]*in[5] - 3*in[2])   ;
    out[3] = -(in[1]*in[5] - in[3])   ;
    out[4] = -(3*in[2] + 5*in[0] - in[4] )   ;
    out[5] = -(5*in[4] + 3*in[2] + in[3] - in[5] * (in[0] + in[1] + 2*in[7] + 1) )   ;
    out[6] = -(5*in[3] + in[1] - 0.5*in[6] )   ;
    out[7] = -(5*in[6] - 2*in[5]*in[7] + in[8] - 0.2*in[7] )   ;
    out[8] = -(2*in[5]*in[7] - in[8])   ;

}

//be care for
int dimIn = 10;

int dimOut = 9;
int noParam = 0;
int MaxDerivativeOrder = 3;

IMap f(_f, dimIn,dimOut,noParam,MaxDerivativeOrder);
IMap fBack(_fBack, dimIn,dimOut,noParam,MaxDerivativeOrder);
int main()
{
    NonlinearSys<double> mysys(f, 9, 0, 9);
    NonlinearSys<double> mysysBack(fBack, 9, 0, 9);

    ReachOptions<double> options;

    //create R0
    Vector_t<double> center(9);
    center << 1, 1, 1, 1, 1, 1, 1, 1, 1;
    Matrix_t<double> generators(Eigen::MatrixXd::Identity(9,9) * 0.01);
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
    options.set_max_error(DBL_MAX*Eigen::MatrixXd::Ones(9,1));

	plt::figure_size(1200, 780);

    clock_t start, end;
	start = clock();
	
	// vector<Zonotope<double>> underR = UnderApprox::underReachClp(mysys, mysysBack, options, R0_, 0.2, 1, 1, 0.002,0.005, 3, 2);
    
    //目前最优 time final 1
    // vector<Zonotope<double>> underR = UnderApprox::underReachClp(mysys, mysysBack, options, R0_, 0.2, 1, 1, 0.002,0.005, 2, 10);

    //能算到0.375s 1000s
    // vector<Zonotope<double>> underR = UnderApprox::underReachClp(mysys, mysysBack, options, R0_, 0.025, 15, 1, 0.005,0.001, 2, 15);
    
    //能算到 0.35s 127.416s
    // vector<Zonotope<double>> underR = UnderApprox::underReachClp(mysys, mysysBack, options, R0_, 0.025, 15, 1, 0.005,0.005, 2, 10);

    //能算到 0.375s 270s time final 2
    vector<Zonotope<double>> underR = UnderApprox::underReachClp(mysys, mysysBack, options, R0_, 0.025, 15, 1, 0.005,0.0025, 2, 10);

    //这都算不到 0.3s
    //vector<Zonotope<double>> underR = UnderApprox::underReachClp(mysys, mysysBack, options, R0_, 0.025, 16, 1, 0.005,0.00025, 2, 20);
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
