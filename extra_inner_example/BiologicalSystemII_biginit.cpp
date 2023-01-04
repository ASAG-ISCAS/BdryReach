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
    Vector_t<double> center(Eigen::MatrixXd::Ones(9,1));
    // center << 1, 1, 1, 1, 1, 1, 1, 1, 1;
    Matrix_t<double> generators(Eigen::MatrixXd::Identity(9,9) * 0.05);
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
	

    //能算到 0.3s 137s 算不到了
    // vector<Zonotope<double>> underR = UnderApprox::underReachClp(mysys, mysysBack, options, R0_, 0.15, 2, 100, 0.015,0.005, 2, 10);

    //能算到 0.36s 165s 算不到了
    // vector<Zonotope<double>> underR = UnderApprox::underReachClp(mysys, mysysBack, options, R0_, 0.12, 3, 100, 0.024,0.004, 2, 10);

    //0.28s 146s 都不能算出来
    // vector<Zonotope<double>> underR = UnderApprox::underReachClp(mysys, mysysBack, options, R0_, 0.14, 2, 100, 0.028,0.0035, 2, 10);

    // vector<Zonotope<double>> underR = UnderApprox::underReachClp(mysys, mysysBack, options, R0_, 0.14, 2, 100, 0.014,0.0014, 10, 20);

    //0.28s 140.55s  能算  0.53 this good
    vector<Zonotope<double>> underR = UnderApprox::underReachClp(mysys, mysysBack, options, R0_, 0.07, 4, 100, 0.007,0.0035, 2, 10);
    //0.28s 149s  
    // vector<Zonotope<double>> underR = UnderApprox::underReachClp(mysys, mysysBack, options, R0_, 0.28, 1, 100, 0.028,0.0035, 2, 10);

    //0.26s 142s 
    // vector<Zonotope<double>> underR = UnderApprox::underReachClp(mysys, mysysBack, options, R0_, 0.13, 2, 100, 0.026,0.00325, 2, 10);

    //0.26s 142s this good 0.62593
    // vector<Zonotope<double>> underR = UnderApprox::underReachClp(mysys, mysysBack, options, R0_, 0.065, 4, 100, 0.0065,0.00325, 2, 10);

    //0.36s 359s
    // vector<Zonotope<double>> underR = UnderApprox::underReachClp(mysys, mysysBack, options, R0_, 0.12, 3, 100, 0.012,0.004, 2, 10);

    //0.26s 141s
    //  vector<Zonotope<double>> underR = UnderApprox::underReachClp(mysys, mysysBack, options, R0_, 0.26, 1, 100, 0.013,0.00325, 5, 10);
    
    //0.3s 1503s
    // vector<Zonotope<double>> underR = UnderApprox::underReachClp(mysys, mysysBack, options, R0_, 0.025, 12, 0.05, 0.005,0.0025, 2, 10);

    //0.3s 1160s
    // vector<Zonotope<double>> underR = UnderApprox::underReachClp(mysys, mysysBack, options, R0_, 0.05, 6, 0.05, 0.005,0.0025, 2, 10);

    //0.3s 2267s
    // vector<Zonotope<double>> underR = UnderApprox::underReachClp(mysys, mysysBack, options, R0_, 0.1, 3, 0.05, 0.01,0.0025, 2, 50);
//     BZF{3} = zonotope([1.215662011080988 0.02808441915126846 0.0008122054952680081 -9.66959947297827e-05 0.001441843713639277 0.001136712119190665 -0.004381961652756991 -4.928472770895583e-05 8.153674485345538e-05 0.01578155508441586 ;
// 0.8764765098795303 -0.001315066136476022 0.0007630659628156076 -0.0001087977337378559 0.001254797559948245 0.02928199784347753 -0.004139535853540989 0.008266457734023722 7.108621064261527e-05 -0.002835855687237778 ;
// 0.7842943343820273 0.01220333562826489 -0.0008122054952679938 9.669599472978587e-05 -0.00144184371363926 -0.001136712119190645 0.004381961652757076 4.928472770896352e-05 -8.153674485345441e-05 0.02352351507051774 ;
// 1.123709989965242 0.001315066136476037 -0.0007630659628155929 0.0001087977337378587 -0.001254797559948229 0.01277613831047486 0.004139535853541056 0.03267387052102734 -7.108621064261436e-05 0.00283585568723782 ;
// 2.899359181666473 0.0489303865757908 0.0001904079477807256 -4.126236592323558e-05 0.0002315276037170325 0.0002427536521591029 0.03389180769665316 -8.954569027653606e-05 1.322981268893048e-05 0.03624252705489166 ;
// 1.989448420918206 0.01830278488445548 -0.00557546746071205 0.000321030592169082 -0.01489717285070686 -0.009233793375747488 0.02968662288098421 -0.00401953430151649 -0.0008308318837001356 0.02475639087539381 ;
// 2.585780604996934 0.0003465216411707128 -0.0003735533989625072 9.28009926489331e-05 0.03919250561895638 0.01915240486155582 0.002073895452024242 0.05171977223242705 -2.414714181755633e-05 0.001281351450985415 ;
// 2.408404085808584 -0.005406775147459044 0.01060345718914483 -0.0002407618888824733 0.0449697098143051 0.01317539493214853 -0.0135673626652681 0.02935919826914984 0.002459805376115289 -0.009791696761487922 ;
// 2.19112857655824 0.005556270234958719 0.007043368624890843 0.0003084234351623942 0.01744773554460484 -0.0006610998856932027 0.01479396812152273 0.01097036467810125 0.01106601311760344 0.01052681163046144 ;]);

    // vector<Zonotope<double>> underR = UnderApprox::underReachClp(mysys, mysysBack, options, R0_, 0.1, 4, 0.05, 0.01,0.0025, 2, 50);
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
