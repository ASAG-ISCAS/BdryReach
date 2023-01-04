/*
 * @Author: your name
 * @Date: 2021-11-11 17:24:27
 * @LastEditTime: 2021-12-30 22:41:19
 * @LastEditors: Please set LastEditors
 * @Description: 打开koroFileHeader查看配置 进行设置: https://github.com/OBKoro1/koro1FileHeader/wiki/%E9%85%8D%E7%BD%AE
 * @FilePath: /Solver/test.cpp
 */

#include <overApprox/overApprox.h>

#include <plotter/matplotlibcpp.h>
#include <plotter/plotter.h>

namespace plt = matplotlibcpp;

using namespace std;
using namespace reachSolver;
using namespace capd;
using capd::autodiff::Node;

double mu = 1.;

void _f(Node/* t*/, Node in[], int /*dimIn*/, Node out[], int/* dimOut*/, Node params[], int noParams){
    out[0] = in[1];
    out[1] = mu * (1-in[0]*in[0])*in[1] - in[0] /*+ in[2]*/;
}

int dimIn = 3; //微分方程的输入维度，由于没有控制 u, 因此输入维度和输出维度相同，都为 2

int dimOut = 2; //微分方程的输出维度
int noParam = 0; //微分方程的参数设置，由于该微分方程没有参数，因此为 0
int MaxDerivativeOrder = 3; //对微分方程进行泰勒展开进行到的最大阶数

IMap f(_f, dimIn,dimOut,noParam,MaxDerivativeOrder); //构成 IMap, 以便于进行区间计算

int main(){
    NonlinearSys<double> mysys(f, 2, 0, 2);
    ReachOptions<double> options;

    //create R0
    Vector_t<double> center(2);
    Vector_t<double> c(2);
    c << 1.4, 2.4;
    Matrix_t<double> generators(2,1);
    Matrix_t<double> G(2,2);
    G<< 0.17,0,
                 0,0.06;
    Zonotope<double> R0_(c,G);

    options.set_R0(R0_);

    options.set_time_step(0.005);
    options.set_taylor_terms(4);
    options.set_zonotope_order(50);
    options.set_intermediate_order(50);
    options.set_error_order(20);
    options.set_alg("lin");
    options.set_tensor_order(3);

    options.set_tFinal(6.74);
    options.set_tStart(0);

    options.set_usekrylovError(1);
    options.set_max_error(DBL_MAX*Eigen::MatrixXd::Ones(2,1));

    ReachableSet<double> R;
    
    vector<ReachableSet<double>> BdReachset = OverApprox::BdReach(mysys, options, R0_);
    
    mysys.reach(options, R);

    plt::figure_size(1200, 780);
    Plotter::plotReach(R, 1, 2, "r");

    for(int i = 0; i < BdReachset.size(); i++){

        Plotter::plotReach(BdReachset[i], 1, 2, "b");
    }
    plt::show();
}

