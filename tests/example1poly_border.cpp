/*
 * @Author: your name
 * @Date: 2021-11-11 17:24:27
 * @LastEditTime: 2021-12-30 22:41:19
 * @LastEditors: Please set LastEditors
 * @Description: 打开koroFileHeader查看配置 进行设置: https://github.com/OBKoro1/koro1FileHeader/wiki/%E9%85%8D%E7%BD%AE
 * @FilePath: /Solver/test.cpp
 */

#include <dynamics/NonlinearSys.h>
#include <dynamics/LinearSys.h>
//#include "include/dynamics/ContDynamics.h"
#include <zonotope/Zonotope.h>
#include <zonotope/polyZonotope.h>
#include <zonotope/globalfunctions.h>

// #include <plotter/matplotlibcpp.h>
// #include <plotter/plotter.h>

/*
#include "../src/dynamics/NonlinearSys.tpp"
#include "../src/dynamics/LinearSys.tpp"
#include "../src/dynamics/ContDynamics.tpp"
#include "../src/zonotope/Zonotope.tpp"

#include "../src/dynamics/ReachableSet.tpp"
#include "../src/dynamics/ReachOptions.tpp"
#include "../src/dynamics/polyReachOptions.tpp"
#include "../src/dynamics/TimeRelated.tpp"
#include "../src/dynamics/ReachableSetElement.tpp"
#include "../src/dynamics/polyReachableSetElement.tpp"
#include "../src/zonotope/polyZonotope.tpp"
*/

#include <iostream>
/*#include <NLF/NLF.h>
#include <Functions.h>
#include <Constants.h>
#include <Utilities.h>*/
#include <eigen3/Eigen/Core>
#include <ctime>
//#include <eigen3/Eigen/Dense>
//#include <eigen3/Eigen/StdVector>

#include <plotter/matplotlibcpp.h>
#include <plotter/plotter.h>

namespace plt = matplotlibcpp;

using namespace std;
using namespace reachSolver;
using namespace capd;
using capd::autodiff::Node;
//using namespace autodiff

double k  = 0.015;
double k2 = 0.01;
double g  = 9.81;

void _f(Node/* t*/, Node in[], int /*dimIn*/, Node out[], int/* dimOut*/, Node params[], int noParams){
	out[0] = in[6] + 0.1 + k2 * (4 - in[5]) - k * sqrt(2 * g) * sqrt(in[0]);
	out[1] = k * sqrt(2 * g) * (sqrt(in[0]) - sqrt(in[1]));
	out[2] = k * sqrt(2 * g) * (sqrt(in[1]) - sqrt(in[2]));
	out[3] = k * sqrt(2 * g) * (sqrt(in[2]) - sqrt(in[3]));
	out[4] = k * sqrt(2 * g) * (sqrt(in[3]) - sqrt(in[4]));
	out[5] = k * sqrt(2 * g) * (sqrt(in[4]) - sqrt(in[5]));
}
/*static INTERVAL_AUTODIFF vanderPolEq(CONST INTERVAL_AUTODIFF &x, int n){
    // differential equation
    switch(n){
        case 1:
            return x(2);
        case 2:
            return mu * (1-x(1)*x(1))*x(2)-x(1);
    }
    return 0.0;
}*/

//be care for
int dimIn			   = 7;
int dimOut			   = 6;
int noParam			   = 0;
int MaxDerivativeOrder = 3;

IMap f(_f, dimIn,dimOut,noParam,MaxDerivativeOrder);

int main(){
    double step = 1;
    double tFinal = 400;
    int N = tFinal/step;
    int dim1 = 1;
    int dim2 = 2;


    NonlinearSys<double> mysys(f, 6, 1, 6);
    NonlinearSys<double> mysys2(f, 6, 1, 6);
    polyReachOptions<double> options;

    ReachOptions<double> options2;

    //create R0
	Vector_t<double>	 center(6);
	center << 2., 4., 4., 2., 10., 4.;
	Matrix_t<double> generators(Eigen::MatrixXd::Identity(6, 6) * 0.2);
	Zonotope<double> R0_(center, generators);

	Vector_t<double> uc(1);
	uc(0, 0) = 0.;
	Matrix_t<double> ug(1, 1);
	ug(0, 0) = 0.005;

    Zonotope<double> U_(uc, ug);

    options.set_time_step(step);
    options.set_taylor_terms(4);
    options.set_zonotope_order(50);
    options.set_intermediate_order(50);
    options.set_error_order(20);
    options.set_alg("poly");
    options.set_tensor_order(3);

    options2.set_time_step(step);
    options2.set_taylor_terms(4);
    options2.set_zonotope_order(50);
    options2.set_intermediate_order(50);
    options2.set_error_order(20);
    options2.set_alg("lin");
    options2.set_tensor_order(3);

    options2.set_maxDepGenOrder(50);
    options2.set_maxPolyZonoRatio(0.01);
    options2.set_restructureTechnique("reducePca");

    options2.set_tFinal(tFinal);
    options2.set_tStart(0);

    options2.set_R0(R0_);
    options2.set_U(U_);

    options2.set_uTrans_nonlin(Eigen::MatrixXd::Zero(1, 1));
    options2.set_usekrylovError(1);
    options2.set_max_error(DBL_MAX * Eigen::MatrixXd::Ones(6, 1));
    //ploy
    options.set_maxDepGenOrder(50);
    options.set_maxPolyZonoRatio(0.01);
    options.set_restructureTechnique("reducePca");

    options.set_tFinal(tFinal);
    options.set_tStart(0);

    options.set_R0(polyZonotope<double>(R0_));
    options.set_U(U_);

    options.set_uTrans_nonlin(Eigen::MatrixXd::Zero(1, 1));
    options.set_usekrylovError(1);
    options.set_max_error(DBL_MAX * Eigen::MatrixXd::Ones(6, 1));
    polyReachableSet<double> R;

    vector<Zonotope<double>> R0_border = Zono_border(R0_);
    
    vector<polyReachableSet<double>> R_border(R0_border.size());
    // cout << R0_border.size() << endl;
    // cout << R0_border << endl;
    // mysys2.reach(options2, R2);

    // cout << "Goof" << endl;
    // mysys.reach(options, R);

    clock_t start, end;
    start = clock();
    mysys.reach(options, R);
    end = clock();

    for(int i = 0; i < R0_border.size(); i++){
        options.set_R0(polyZonotope<double>(R0_border[i]));
        mysys.reach(options, R_border[i]);
    }

    // int select = 0;

    // options.set_R0(polyZonotope<double>(R0_border[select]));
    // mysys.reach(options, R_border[select]);

    // cout << isFullDim(polyZonotope<double>(R0_));
    // cout << R0_border[select] << endl;
    /*cout<<"time cost: "<<(double)(end-start)/CLOCKS_PER_SEC<<endl;
    int a[5] = {10,100,200,300};

    //Matrix_t<double> abc;
    //cout<<abc.size()<<endl;
    cout<<R.Rtime_point().set().size()<<endl;
    for(auto i =0;i<4;i++){
        cout<<a[i]<<endl;
        cout<<R.Rtime_point().set()[a[i]][0].rs().center()<<endl;
        cout<<"Genertors:"<<"double "<<R.Rtime_point().set()[a[i]][0].rs().generators().cols()<<endl;
    }*/
    for(int k = 1; k < 3; k++){
        plt::figure_size(1200, 780); 

        Plotter::plotReach(R,2, 1+2*k, 2+2*k, "r");

        for(int j = 0; j < R_border.size(); j++){
           Plotter::plotReach(R_border[j],2, 1+2*k, 2+2*k, "b");
        }

        plt::show();
    }

    // Plotter::plotReach(R,2, 1, 2, "r");
    // Plotter::plotReach(R_border[0],2, 1, 2, "b");
    // Plotter::plotReach(R_border[1],2, 1, 2, "b");
    // Plotter::plotOvertime(R_border[0], 1, "b");
    // cout << R_border[0].Rtime_point().set()[399][0].rs();

            // for(int i = 0; i < R0_border.size(); i++){
            //     cout << R0_border[i] << endl;
            //     cout << R_border[i].Rtime_point().set()[399][0].rs();
            // }
    //     for(int k = 0; k < 3; k++){
    //     plt::figure_size(1200, 780);        
    //     Plotter::plotZonotope(R0_,dim1 + 2*k,dim2 + 2*k, "r");
        // for(int i = 0; i < R_border.size(); i++){
        //     // plt::figure_size(1200, 780);
        //     // Plotter::plotZonotope(R0_,dim1 + 2*k,dim2 + 2*k, "r");
        //     // Plotter::plotZonotope(R0_border[i],dim1 + 2*k,dim2 + 2*k, "b");
        //     // plt::show();
        // }
    //     plt::show();
    // }
    

    // polyZonotope<double> polyZ = R.Rtime_point().set()[1347][0].rs();
    // Zonotope<double> Z = R2.Rtime_point().set()[1347][0].rs();
    // Plotter::plotpolyZonotope(polyZ,10);
    // Plotter::plotZonotope(Z);
    // plt::show();
    
    // cout << R_border[select].Rtime_point().set()[N-1][0].rs();
    for(int j = 0; j < N; j++){
        // cout << R_border[select].Rtime_point().set()[j][0].rs();
            }

    // cout<<"time cost: "<<(double)(end-start)/CLOCKS_PER_SEC<<endl;
    /*INTERVAL_VECTOR x(34);
    for(auto i=0; i<34; i++)
    {
        x.theElements[i].inf = (double)(i*1.1);
        x.theElements[i].sup = (double)(i*2.1);
    }
    //x.theElements[6] = 2.0;
    cout<<x<<endl;
    cout<<VectorFunction(FunctionTest,x)<<endl;*/
       
                
}

