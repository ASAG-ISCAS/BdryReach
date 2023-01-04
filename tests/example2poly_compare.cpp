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
double mu = 1.;

void _f(Node/* t*/, Node in[], int /*dimIn*/, Node out[], int/* dimOut*/, Node params[], int noParams){
    out[0] = in[1];
    out[1] = mu * (1-in[0]*in[0])*in[1] - in[0] /*+ in[2]*/;
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
int dimIn = 2;

int dimOut = 2;
int noParam = 0;
int MaxDerivativeOrder = 3;

IMap f(_f, dimIn,dimOut,noParam,MaxDerivativeOrder);

int main(){
    NonlinearSys<double> mysys(f, 2, 0, 2);
    NonlinearSys<double> mysys2(f, 2, 0, 2);
    polyReachOptions<double> options;

    ReachOptions<double> options2;

    //create R0
    Vector_t<double> center(2);
    center << 1.4, 2.4;
    Matrix_t<double> generators(2,2);
    generators<< 0.17,0,
                 0,0.06;
    Zonotope<double> R0_(center,generators);

    options.set_time_step(0.005);
    options.set_taylor_terms(4);
    options.set_zonotope_order(50);
    options.set_intermediate_order(50);
    options.set_error_order(20);
    options.set_alg("poly");
    options.set_tensor_order(3);

    options2.set_time_step(0.005);
    options2.set_taylor_terms(4);
    options2.set_zonotope_order(50);
    options2.set_intermediate_order(50);
    options2.set_error_order(20);
    options2.set_alg("lin");
    options2.set_tensor_order(3);

    options2.set_maxDepGenOrder(50);
    options2.set_maxPolyZonoRatio(0.01);
    options2.set_restructureTechnique("reducePca");

    options2.set_tFinal(6.74);
    options2.set_tStart(0);

    options2.set_R0(R0_);

    options2.set_usekrylovError(1);
    options2.set_max_error(DBL_MAX*Eigen::MatrixXd::Ones(2,1));

    //ploy
    options.set_maxDepGenOrder(50);
    options.set_maxPolyZonoRatio(0.01);
    options.set_restructureTechnique("reducePca");

    options.set_tFinal(6.74);
    options.set_tStart(0);

    options.set_R0(polyZonotope<double>(R0_));

    options.set_usekrylovError(1);
    options.set_max_error(DBL_MAX*Eigen::MatrixXd::Ones(2,1));
    polyReachableSet<double> R;

    ReachableSet<double> R2;

    mysys2.reach(options2, R2);

    // cout << "Goof" << endl;
    mysys.reach(options, R);

    /*clock_t start, end;
    start = clock();
    mysys.reach(options, R);
    end = clock();
    cout<<"time cost: "<<(double)(end-start)/CLOCKS_PER_SEC<<endl;
    int a[5] = {10,100,200,300};

    //Matrix_t<double> abc;
    //cout<<abc.size()<<endl;
    cout<<R.Rtime_point().set().size()<<endl;
    for(auto i =0;i<4;i++){
        cout<<a[i]<<endl;
        cout<<R.Rtime_point().set()[a[i]][0].rs().center()<<endl;
        cout<<"Genertors:"<<"double "<<R.Rtime_point().set()[a[i]][0].rs().generators().cols()<<endl;
    }*/
    plt::figure_size(1200, 780);
    // for(int i = 0; i < 1348; i++){
    //     Zonotope<double> Z = R2.Rtime_point().set()[i][0].rs();
    //     Plotter::plotZonotope(Z,1,2);
    //     // cout << polyZ << endl;
    // }
    // for(int i = 0; i < 1348; i++){
    //     polyZonotope<double> polyZ = R.Rtime_point().set()[i][0].rs();
    //     Plotter::plotpolyZonotope(polyZ,2,1,2);
    //     // cout << polyZ << endl;
    // }
    Plotter::plotReach(R2, 1, 2,"r");
    Plotter::plotReach(R, 2, 1, 2,"b");
    plt::show();
    // polyZonotope<double> polyZ = R.Rtime_point().set()[1347][0].rs();
    // Zonotope<double> Z = R2.Rtime_point().set()[1347][0].rs();
    // Plotter::plotpolyZonotope(polyZ,10);
    // Plotter::plotZonotope(Z);
    // plt::show();
    // cout << R.Rtime_point().set()[1346][0].rs();

    
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

    // Matrix_t<double> A(2,2);
    // A<< 0.17,0,
    //     0,0.06;

    // Matrix_t<double> O;
    // O.setZero(2,2);

    // Matrix_t<bool> res= (A.array() == O.array());

    // // res = (A.array() == O.array());
    // cout << findZero(A);

    // vector<int> a={1,2,3,4,5};
    // vector<int> b={2,4};

    // cout << setdiff(a,b);
    
    // Matrix_t<double> A(7,3);
    // A << 1, 2, 3,
    // 2,3,4,
    // 2,1,5,
    // 3,4,1,
    // 5,4,2,
    // 1,2,3,
    // 3,4,2;
    // cout << A % 2;
    //   vector<double> res = {1,2,3,4,5,6};
    //   vector<double> b = {4,6};
    // cout << setdiff(res,b);
}

