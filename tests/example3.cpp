/*
 * @Author: your name
 * @Date: 2021-11-11 17:24:27
 * @LastEditTime: 2021-12-30 17:48:35
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

// #include <"../src/>dynamics/NonlinearSys.tpp"
// #include "../src/dynamics/LinearSys.tpp"
// #include "../src/dynamics/ContDynamics.tpp"
// #include "../src/zonotope/Zonotope.tpp"
// #include "../src/zonotope/polyZonotope.tpp"

// #include "../src/dynamics/ReachableSet.tpp"
// #include "../src/dynamics/ReachOptions.tpp"
// #include "../src/dynamics/TimeRelated.tpp"
// #include "../src/dynamics/ReachableSetElement.tpp"

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


using namespace std;
using namespace reachSolver;
using namespace capd;
using capd::autodiff::Node;
//using namespace autodiff
double mu = 1.;

void _f(Node/* t*/, Node in[], int /*dimIn*/, Node out[], int/* dimOut*/, Node params[], int noParams){
    out[0] = in[1];
    out[1] = mu * (1-in[0]*in[0])*in[1] - in[0];
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
int dimIn = 3;

int dimOut = 2;
int noParam = 0;
int MaxDerivativeOrder = 3;

IMap f(_f, dimIn,dimOut,noParam,MaxDerivativeOrder);

int main(){
    NonlinearSys<double> mysys(f, 2, 0, 2);
    ReachOptions<double> options;

    //create R0
    Vector_t<double> center(2);
    center << 1.4, 2.3;
    // center << 1.25, 2.3;
    Matrix_t<double> generators(2,2);
    generators<< 0.3,0,
                 0,0.05;
        // generators<< 0.15,0,
        //          0,0.05;
    Zonotope<double> R0_(center,generators);

    options.set_time_step(0.02);
    options.set_taylor_terms(4);
    options.set_zonotope_order(10);
    options.set_intermediate_order(10);
    options.set_error_order(5);
    options.set_alg("lin");
    options.set_tensor_order(3);
    options.set_tFinal(3.2);
    options.set_tStart(0);
    options.set_R0(R0_);
    options.set_verbose(true);
    options.set_reductionInterval(100);

    options.set_usekrylovError(1);
    options.set_max_error(0.05*Eigen::MatrixXd::Ones(2,1));

    //options.set_uTrans_nonlin(0);
    //options.set_usekrylovError(1);
    //options.set_max_error(DBL_MAX*Eigen::MatrixXd::Ones(6,1));


    // ValidateOptions
    
    //ReachableSet<double> params;
    //params.set_t
    //cout<<"1"<<endl;

    // ReachableSet<double> R;
    vector<ReachableSet<double>> R;

    //mysys.reach(options, R);

    clock_t start, end;
    start = clock();
    //mysys.reach(options, R);
    mysys.reachSplit(options, R);
    end = clock();
    cout<<"time cost: "<<(double)(end-start)/CLOCKS_PER_SEC<<endl;
    int a[5] = {10,100,200,300};

    Plotter::plotVecReach(R, 1, 2, "b");


    //Matrix_t<double> abc;
    //cout<<abc.size()<<endl;
    //cout<<R.Rtime_point().set().size()<<endl;
    /*for(auto i =0;i<4;i++){
        cout<<a[i]<<endl;
        cout<<R.Rtime_point().set()[a[i]][0].rs().center()<<endl;
        cout<<"Genertors:"<<"Number "<<R.Rtime_point().set()[a[i]][0].rs().generators().cols()<<endl;
    }*/

    //cout<<R.Rtime_point().set()[0][0].rs()<<endl;

    
    /*INTERVAL_VECTOR x(34);
    for(auto i=0; i<34; i++)
    {
        x.theElements[i].inf = (double)(i*1.1);
        x.theElements[i].sup = (double)(i*2.1);
    }
    //x.theElements[6] = 2.0;
    cout<<x<<endl;
    cout<<VectorFunction(FunctionTest,x)<<endl;*/
    

   /* clock_t start, end;
    start = clock();
    mysys.reach(options, R);
    end = clock();
    cout<<"time cost: "<<(double)(end-start)/CLOCKS_PER_SEC<<endl;
    int a[5] = {10,100,200,300};
    //cout<<R.Rtime_point().set().size()<<endl;
    for(auto i =0;i<4;i++){
        cout<<a[i]<<endl;
        cout<<R.Rtime_point().set()[a[i]][0].rs()<<endl;
    }*/



}

