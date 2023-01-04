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
    out[0] = in[1];
    out[1] = (0.2 - 0.7*sin(in[0]) - 0.05 * in[1]);
}

//be care for
int dimIn			   = 2;
int dimOut			   = 2;
int noParam			   = 0;
int MaxDerivativeOrder = 3;

IMap f(_f, dimIn,dimOut,noParam,MaxDerivativeOrder);

int main(){
    double step = 0.01;
    double tFinal = 6;
    int N = tFinal/step;
    int dim1 = 1;
    int dim2 = 2;


    NonlinearSys<double> mysys(f, 2, 0, 2);

    ReachOptions<double> options;

    //create R0
	Vector_t<double>	 center(2);
    Matrix_t<double> generators(2,2);
	center << 0,3;
    generators<< 0.1,0,
        0,0.1;

	Zonotope<double> R0_(center, generators);

    options.set_time_step(step);
    options.set_taylor_terms(4);
    options.set_zonotope_order(50);
    options.set_intermediate_order(50);
    options.set_error_order(20);
    options.set_alg("lin");
    options.set_tensor_order(3);

    options.set_tFinal(tFinal);
    options.set_tStart(0);

    options.set_R0(R0_);

    options.set_uTrans_nonlin(Eigen::MatrixXd::Zero(1, 1));
    options.set_usekrylovError(1);
    options.set_max_error(DBL_MAX * Eigen::MatrixXd::Ones(2, 1));
    ReachableSet<double> R;

    vector<Zonotope<double>> R0_border = Zono_border(R0_);
    
    vector<ReachableSet<double>> R_border(R0_border.size());
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
        options.set_R0(R0_border[i]);
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
    for(int k = 0; k < 1; k++){
        plt::figure_size(1200, 780); 

        Plotter::plotReach(R, 1+2*k, 2+2*k, "r");

        for(int j = 0; j < R_border.size(); j++){
           Plotter::plotReach(R_border[j], 1+2*k, 2+2*k, "b");
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

