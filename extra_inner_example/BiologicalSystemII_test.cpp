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
    center << 
1.205529535570056  ,  
0.8301124684740934 ,
0.7944893105167011 ,
1.169886238658259  ,
3.2984050612386    ,
2.135093993661071  ,
2.975549615668121  ,
2.791996906132318  ,
2.810657126051584 ;
 Matrix_t<double> generators(9,9);

 generators <<
0.005194235501322227, 0.000476280413435716, 0.0002452359923434782, 0.0003241446740046792 ,0.0001056989906948412, -9.836549149978793e-05 ,0.003368312491843222, 4.274328339627506e-05, -0.001068163932012729 ,
-0.0004823864333136317, 0.0004106152082608192,0.0002318407878424212, 0.005518583501895711, 0.002053472646357613, -0.0001153184313934239, -0.0008023855331421122, 3.695143419883567e-05, -0.00101376103527014 ,
0.003163654926256969, -0.0004762804134357025, -0.0002452359923434667, -0.0003241446740046653, -0.0001056989906948371, 9.836549149979752e-05, 0.004876742586267243, -4.274328339627393e-05, 0.001068163932012787,
0.0004823864333136437, -0.0004106152082608058, -0.0002318407878424096, 0.003181589169794427, 0.00660774145169918, 0.0001153184313934317 ,0.0008023855331421428, -3.695143419883451e-05, 0.001013761035270185 ,
0.0120245773955466, 9.983725832841613e-05, 7.670959817200691e-05, 8.700087164701691e-05, -1.142459218644073e-05 ,-5.793111321619528e-05, 0.00932511314321929, 9.085210240121367e-06, 0.006206585976369049 ,
0.005470401964000133, -0.004112954166680764, -0.001417790743916263, -0.002398815794916746, -0.002114116395813115, 0.0002874138692582876, 0.006335514880530054, -0.0003627012476345231, 0.006219815680668204,
0.0001796840682894773, 0.007654404519967284, -0.0001523813061399372, 0.005316243324529729, 0.01316632796209107 ,0.0001353588806509463, 0.0004648475997802013, -1.656150532526665e-05, 0.0006772036408841903,
-0.002226950430360349, 0.01066816170066687, 0.002460699799081726 ,0.004181613034083179, 0.009065582081214673 ,-0.0002403445889932815, -0.003122673898375732, 0.0009202228714314135, -0.003620622246973415 ,
0.002327926863381532, 0.00494608641882128, 0.001860711294127165 ,9.161997573077946e-05, 0.003878086374607554, 0.0003651055969272838 ,0.003451281165947576, 0.003529387506495147, 0.004120223321866977 ;
    // Vector_t<double> center(9);
    // center << 1, 1, 1, 1, 1, 1, 1, 1, 1;
    // Matrix_t<double> generators(Eigen::MatrixXd::Identity(9,9) * 0.01);
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
    
    //目前最优
    // vector<Zonotope<double>> underR = UnderApprox::underReachClp(mysys, mysysBack, options, R0_, 0.2, 1, 1, 0.002,0.005, 2, 10);

    //能算到0.375s 1000s
    // vector<Zonotope<double>> underR = UnderApprox::underReachClp(mysys, mysysBack, options, R0_, 0.025, 15, 1, 0.005,0.001, 2, 15);
    
    //能算到 0.35s 127.416s
    // vector<Zonotope<double>> underR = UnderApprox::underReachClp(mysys, mysysBack, options, R0_, 0.025, 15, 1, 0.005,0.005, 2, 10);

    //能算到 0.375s 270s
    // vector<Zonotope<double>> underR = UnderApprox::underReachClp(mysys, mysysBack, options, R0_, 0.025, 15, 1, 0.005,0.0025, 2, 10);

    vector<Zonotope<double>> underR = UnderApprox::underReachClp(mysys, mysysBack, options, R0_, 0.00625, 4, 1, 0.00005,0.00005, 2, 20);
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
