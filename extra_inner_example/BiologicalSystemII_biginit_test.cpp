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
1.216573899056914 ,
0.8916869452344006,
0.7834408718105307,
1.108372595441522 ,
2.759176445318996 ,
1.931046899135812 ,
2.455720951906271 ,
2.281665613961302 ,
2.014070311541179;

    // center << 1, 1, 1, 1, 1, 1, 1, 1, 1;
    Matrix_t<double> generators(9,9);
generators << 
0.03074950117082809, 0.01631409864342757, -0.004044365605425073, -0.001062026971074542 ,-0.0001658860865607907, 0.001202222372379981 ,0.001026292747709442, 0.001002088070749567, 0.0001107159346542772,
 -0.001083548559115314, -0.002638010927325588, -0.003816579705671562, -0.001178658603102678 ,0.008067586963264916, 0.001050738614287564, 0.03105494078696434, 0.0009403350786188515, 9.6911940592711e-05,
 0.01178170688876852, 0.02505467064407382, 0.004044365605425207, 0.001062026971074602 ,0.0001658860865608057, -0.001202222372379956 ,-0.001026292747709414, -0.001002088070749533, -0.0001107159346542749,
0.001083548559115333, 0.002638010927325664, 0.003816579705671672, 0.00117865860310273, 0.0344620872266883 ,-0.001050738614287541 ,0.01213949874888123,-0.0009403350786188177, -9.691194059270909e-05,
0.04820424460143667, 0.03513533658459009, 0.03526496334670005, -0.00040206893439028 ,-8.974245407542865e-05, 0.0001746191337394943, 0.0002000632099355064, 0.0002109179093739291, 1.623155963667378e-05,
0.01654696140922944, 0.02416328630865737, 0.02928394082510272, 0.003733253999085734 ,-0.002486640209935811, -0.01338161354787081, -0.008743224803353135, -0.007397660681504518, -0.001216854590854982,
0.0002463350838369035, 0.001082636434404617, 0.001715701691693238, 0.0008931073023190398, 0.04989191999066934, 0.0406138300933599, 0.01755828914122462 ,-0.0004124178059794841, -2.967158798024467e-05,
-0.004284820119818835, -0.008750447797031137, -0.01213092599162908, -0.002640529673335155 ,0.02630061142426849, 0.04335912282193612 ,0.01144181577891717, 0.01479701292582204, 0.003871418568819538,
0.00437982398097074, 0.009324883572366145, 0.01306586319127726, 0.003237400472740221, 0.009259202902013283 ,0.01551914354487106,-0.0008536846578429158, 0.009213113139147935,0.01853357745547676 ;
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
	

    

    //0.275s 的初始区域 算0.3s 的

    vector<Zonotope<double>> underR = UnderApprox::underReachClp(mysys, mysysBack, options, R0_, 0.025, 1, 0.05, 0.005,0.0025, 2, 10);

    end = clock();
	cout << "time cost: " << (double) (end - start) / CLOCKS_PER_SEC << endl;

    for(int i = 1; i < underR.size() ; i++){
		cout << "BZF{" << i << "} = zonotope(" << underR[i] << ");" << endl;
	}

	for(int i = 0; i < underR.size(); i++){
		Plotter::plotZonotope(underR[i], 1, 2, "g");
	}
	plt::show();
}
