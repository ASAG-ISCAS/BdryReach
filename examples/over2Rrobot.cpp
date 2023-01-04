/*
 * @Author: your name
 * @Date: 2021-11-11 17:24:27
 * @LastEditTime: 2021-12-28 18:18:30
 * @LastEditors: Please set LastEditors
 * @Description: 打开koroFileHeader查看配置 进行设置:
 * https://github.com/OBKoro1/koro1FileHeader/wiki/%E9%85%8D%E7%BD%AE
 * @FilePath: /Solver/test.cpp
 */

#include <overApprox/overApprox.h>
#include <plotter/matplotlibcpp.h>
#include <plotter/plotter.h>

#include <ctime>
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

int dimIn			   = 5;
int dimOut			   = 4;
int noParam			   = 0;
int MaxDerivativeOrder = 3;

IMap f(_f, dimIn, dimOut, noParam, MaxDerivativeOrder);

int main()
{
	NonlinearSys<double> mysys(f, 4, 0, 4);
	ReachOptions<double> options;
	Vector_t<double>	 center(4);
	// center << 0.5, 0.02, 0.6, 0.;
	center << 0, 0, 0, 0.;
	// Matrix_t<double> generators(Eigen::MatrixXd::Identity(4, 4) * 0.02);
	Vector_t<double> gdrag(4);
	gdrag << 0.1,0.1,0.1,0.1;

	Matrix_t<double> generators = gdrag.asDiagonal();

	Zonotope<double> R0_(center, generators);

	options.set_time_step(0.01);
	options.set_taylor_terms(4);
	options.set_zonotope_order(50);
	options.set_alg("lin");
	options.set_tensor_order(2);
	options.set_tFinal(0.5);
	options.set_tStart(0);
	options.set_R0(R0_);
	options.set_uTrans_nonlin(Eigen::MatrixXd::Zero(1, 1));
	options.set_usekrylovError(1);
	options.set_max_error(DBL_MAX * Eigen::MatrixXd::Ones(4, 1));

	ReachableSet<double> R;

	// mysys.reach(options, R);

	vector<Zonotope<double>> RStemp;

	RStemp.push_back(R0_);

	vector<Zonotope<double>> R0Split = UnderApprox::boundSplit(RStemp, 0.02);

	cout << "All splits: " << R0Split.size() << endl;

	vector<ReachableSet<double>> RSplit(R0Split.size());

	clock_t start, end;
	start = clock();

	for(int i = 0; i < R0Split.size(); i++){

		options.set_R0(R0Split[i]);
		mysys.reach(options, RSplit[i]);
	}

	end = clock();

	double allTime = (double) (end - start) / CLOCKS_PER_SEC;

	cout << "All time cost: " << (double) (end - start) / CLOCKS_PER_SEC << endl;

	start = clock();

	vector<ReachableSet<double>> BdReachset = OverApprox::BdReachSplit(mysys, options, R0_, 0.03);

	end = clock();

	double boundTime = (double) (end - start) / CLOCKS_PER_SEC;

	cout << "Bound time cost: " << (double) (end - start) / CLOCKS_PER_SEC << endl;

    plt::figure_size(1200, 780);
    // Plotter::plotReach(R, 1, 3, "r");

	for(int i = 0; i < BdReachset.size(); i++){

		Plotter::plotReach(BdReachset[i], 1, 3, "b");
	}

	for(int i = 0; i < RSplit.size(); i++){

        Plotter::plotReach(RSplit[i], 1, 3, "r");
    }

	// Vector_t<double>	 c(2);
	// c << -0.6,2.25;
	// Matrix_t<double> G(2,2);
	// G << 0.05, 0,
	// 	0, 0.25;

	Vector_t<double>	 c(2);
	c << -0.6,1.5;
	Matrix_t<double> G(2,2);
	G << 0.05, 0,
		0, 0.25;


	Zonotope<double> safe(c,G);

	Plotter::plotZonotope(safe, 1, 2, "tab:orange", "2");

    plt::show();

	cout << "All splits: " << R0Split.size() << endl;
	cout << "Bound splits: " << BdReachset.size() << endl;
	cout << "All time cost: " << allTime << endl;
	cout << "Bound time cost: " << boundTime << endl;

}
