#include <dynamics/NonlinearSys.h>
#include <dynamics/LinearSys.h>
#include <zonotope/Zonotope.h>
#include <zonotope/polyZonotope.h>
#include <zonotope/globalfunctions.h>

#include <iostream>
#include <eigen3/Eigen/Core>
#include <ctime>

#include <plotter/matplotlibcpp.h>
#include <plotter/plotter.h>

using namespace std;
using namespace reachSolver;
using namespace capd;
using capd::autodiff::Node;
//using namespace autodiff

double g = 9.81;
double m = 1;
double l = 1;
double u;
double q = 3.14/2;
double kp = 100;
double ki = 0;
double kd = 10;

void _f(Node /* t*/, Node in[], int /*dimIn*/, Node out[], int /* dimOut*/, Node params[], int noParams)
{
	out[0] = in[1];
	out[1] = 3*(u - 0.5*m*g*cos(in[0]))/(m*l*l);
}



//be care for
int dimIn = 2;
int dimOut = 2;
int noParam = 0;
int MaxDerivativeOrder = 3;


int main(){
	double tStep = 0.0005;
	double uStep = 0.01;
	int N = uStep/tStep; 
	vector<Zonotope<double>> R0(51);
    ReachOptions<double> options;//set ReachOptions ***

    //create R0
    Vector_t<double> center(2);
    center << 0, 0;
    Matrix_t<double> generators(2,2);
    generators<< 0,0.1,
                 0.1,0;
    Zonotope<double> R0_(center,generators);
	R0[0] = R0_;

    options.set_time_step(tStep);
    options.set_taylor_terms(4);
    options.set_zonotope_order(50);
    options.set_intermediate_order(50);
    options.set_error_order(20);
    options.set_alg("lin"); // set alg ***
    options.set_tensor_order(3);

    //ploy
    options.set_maxDepGenOrder(50);
    options.set_maxPolyZonoRatio(0.01);
    options.set_restructureTechnique("reducePca");

    options.set_tFinal(uStep);
    options.set_tStart(0);

    //options.set_R0(R0_);//set zonotope R0 ***

    options.set_usekrylovError(1);
    options.set_max_error(DBL_MAX*Eigen::MatrixXd::Ones(2,1));
    plt::figure_size(1200, 780);
    double Isum = 0;
	for(auto i = 0; i < 50; i++){
        Isum += R0[i].center()(0);
        u = kp * (q - R0[i].center()(0)) + ki * (q*(i+1) - Isum) * uStep + kd * (- R0[i].center()(1));
        IMap f(_f, dimIn,dimOut,noParam,MaxDerivativeOrder);
        NonlinearSys<double> mysys(f, 2, 0, 2);
		options.set_R0(R0[i]);//set zonotope R0 ***
		ReachableSet<double> R;//set zonotope R ***

		mysys.reach(options, R);

		for(auto j = 0; j < N; j++){
			Plotter::plotZonotope(R.Rtime_point().set()[j][0].rs());
		}
		R0[i+1] = R.Rtime_point().set()[N-1][0].rs();
	}
	plt::show();
	//cout << R.Rtime_point().set()[799][0].rs() << endl;
	//第10s，第0个ReachableSetElement的Zonotope
}

