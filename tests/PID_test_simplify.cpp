/*
 * @Author: your name
 * @Date: 2021-11-11 17:24:27
 * @LastEditTime: 2021-12-28 18:18:30
 * @LastEditors: Please set LastEditors
 * @Description: 打开koroFileHeader查看配置 进行设置:
 * https://github.com/OBKoro1/koro1FileHeader/wiki/%E9%85%8D%E7%BD%AE
 * @FilePath: /Solver/test.cpp
 */

#include <capd/capdlib.h>
#include <ctime>
#include <zonotope/Zonotope.h>
#include <dynamics/LinearSys.h>
#include <dynamics/NonlinearSys.h>
#include <eigen3/Eigen/Core>
#include <functional>
#include <iostream>

#include <plotter/matplotlibcpp.h>
#include <plotter/plotter.h>

namespace plt = matplotlibcpp;

using namespace std;
using namespace reachSolver;
using namespace capd;
using capd::autodiff::Node;

double k  = 0.015;
double k2 = 0.01;
double g  = 9.81;
double u1 = 2;
double u2 = 2;
void _f(Node /* t*/, Node in[], int /*dimIn*/, Node out[], int /* dimOut*/, Node params[], int noParams)
{
out[0]=in[1];
out[1]=-3*(4*sin(in[2])*in[1]*in[3] + 3*cos(in[2])*sin(in[2])*in[1]*in[1] + 2*sin(in[2])*in[1]*in[1] + 2*sin(in[2])*in[3]*in[3] - 6*cos(in[2])*u2 + 4*u1 - 4*u2)/(9*cos(in[2])*cos(in[2]) - 16);
out[2]=in[3];
out[3]=3*(6*cos(in[2])*sin(in[2])*in[1]*in[3] + 4*sin(in[2])*in[1]*in[3] + 6*cos(in[2])*sin(in[2])*in[1]*in[1] + 3*cos(in[2])*sin(in[2])*in[3]*in[3] + 10*sin(in[2])*in[1]*in[1] + 2*sin(in[2])*in[3]*in[3] + 6*cos(in[2])*u1 - 12*cos(in[2])*u2 + 4*u1 - 20*u2)/(9*cos(in[2])*cos(in[2]) - 16);
}

/*static INTERVAL_AUTODIFF tank6Eq(CONST INTERVAL_AUTODIFF &x, int n){
    // differential equation
    switch(n){
        case 1:
            return x(7)+0.1+k2*(4-x(6))-k*Sqrt(2*g)*Sqrt(x(1));
        case 2:
            return k*Sqrt(2*g)*(Sqrt(x(1))-Sqrt(in[1]));
        case 3:
            return k*Sqrt(2*g)*(Sqrt(in[1])-Sqrt(in[2]));
        case 4:
            return k*Sqrt(2*g)*(Sqrt(in[2])-Sqrt(in[3]));
        case 5:
            return k*Sqrt(2*g)*(Sqrt(in[3])-Sqrt(x(5)));
        case 6:
            return k*Sqrt(2*g)*(Sqrt(x(5))-Sqrt(x(6)));
    }
    return 0.0;
}

FUNCTION FunctionTest (7, 6, tank6Eq);*/
int dimIn			   = 4;
int dimOut			   = 4;
int noParam			   = 0;
int MaxDerivativeOrder = 3;

IMap f(_f, dimIn, dimOut, noParam, MaxDerivativeOrder);

int main()
{
	NonlinearSys<double> mysys(f, 4, 0, 4);
	ReachOptions<double> options;
	Vector_t<double>	 center(4);
	center << 0., 0., 0., 0.;
	Matrix_t<double> generators(Eigen::MatrixXd::Identity(4, 4) * 0);
	Zonotope<double> R0_(center, generators);

	Vector_t<double> uc(1);
	uc(0, 0) = 0.;
	Matrix_t<double> ug(1, 1);
	ug(0, 0) = 0.005;

	Zonotope<double> U_(uc, ug);

	options.set_time_step(0.01);
	options.set_taylor_terms(4);
	options.set_zonotope_order(50);
	options.set_alg("lin");
	options.set_tensor_order(3);
	options.set_tFinal(0.7);
	options.set_tStart(0);
	options.set_R0(R0_);
	options.set_uTrans_nonlin(Eigen::MatrixXd::Zero(1, 1));
	options.set_usekrylovError(1);
	options.set_max_error(DBL_MAX * Eigen::MatrixXd::Ones(4, 1));

	// ValidateOptions

	// ReachableSet<double> params;
	// params.set_t
	// cout<<"1"<<endl;

	ReachableSet<double> R;

	/* VECTOR x(7);
       for(auto i=0; i<6; i++)
       {
           x.theElements[i] = center(i,0);
       }
       x.theElements[6] = 2.0;*/
	// cout<<x<<endl;
	// cout<<VectorFunction(FunctionTest,x)<<endl;
	//
	clock_t start, end;
	start = clock();
	mysys.reach(options, R);
	end = clock();
	cout << "time cost: " << (double) (end - start) / CLOCKS_PER_SEC << endl;
	int a[5] = {10, 100, 200, 300};

	Matrix_t<double> abc;
	cout << abc.size() << endl;
	// cout<<R.Rtime_point().set().size()<<endl;
	for (auto i = 0; i < 4; i++)
	{
		// cout<<a[i]<<endl;
		// cout<<R.Rtime_point().set()[a[i]][0].rs()<<endl;
	}
	//cout << R.Rtime_point().set()[10][0].rs() << endl;
	plt::figure_size(1200, 780);
	for(auto i = 0; i < 70; i++){
		//cout << R.Rtime_point().set()[i][0].rs() << endl;
		Plotter::plotZonotope(R.Rtime_point().set()[i][0].rs());
	}
	plt::show();
	//cout << R.Rtime_point().set()[10][0].rs() << endl;
	//第10s，第0个ReachableSetElement的Zonotope
	return 0;
	//--------------------------------------------------------------------------------------------------------------

	/*
  IntervalMatrix m1;
  IntervalMatrix m2;
  IntervalMatrix MA1;
  IntervalMatrix MA2;
  Matrix_t<double> m1sup(6,1);
  m1.sup = m1sup;
  Matrix_t<double> m1inf(6,1);
  m1.inf = m1inf;
  m1.inf<<1.,2.,3.,4.,5.,6.;
  m1.sup<<1.,2.,3.,4.,5.,6.;
  Matrix_t<double> m2sup(1,1);
  m2.sup = m2sup;
  Matrix_t<double> m2inf(1,1);
  m2.inf = m2inf;
  m2.sup<<8.;
  m2.inf<<7.;

  Matrix_t<double> MA1sup(6,6);
  MA1.sup = MA1sup;
  Matrix_t<double> MA1inf(6,6);
  MA1.inf = MA1inf;

  Matrix_t<double> MA2sup(6,1);
  MA2.sup = MA2sup;
  Matrix_t<double> MA2inf(6,1);
  MA2.inf = MA2inf;
  //m1.inf(6,2);
  cout<<m1.sup(0,0)<<endl;
  cout<<m1.inf(0,0)<<endl;

  //m1.inf(0,0) = 100.;
  //cout<<m1.inf(0,0)<<endl;
  cout<<"eadadad"<<endl;

  double nmd = k*Sqrt(2*g)/Sqrt(1.)/Sqrt(1.)/Sqrt(1.)/4.;
  cout<<nmd<<endl;*/
}
