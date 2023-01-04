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
#include <underApprox/underApprox.h>

#include <glpk.h>
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

//be care for
int dimIn = 3;

int dimOut = 2;
int noParam = 0;
int MaxDerivativeOrder = 3;

IMap f(_f, dimIn,dimOut,noParam,MaxDerivativeOrder);

int main()
{
    NonlinearSys<double> mysys(f, 2, 0, 2);
    ReachOptions<double> options;

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
    options.set_alg("lin");
    options.set_tensor_order(3);

	options.set_tFinal(0.005);
    options.set_tStart(0);

    options.set_R0(R0_);

    options.set_usekrylovError(1);
    options.set_max_error(DBL_MAX*Eigen::MatrixXd::Ones(2,1));



	// ReachableSet<double> R;

	// clock_t start, end;
	// start = clock();
	// mysys.reach(options, R);
	// end = clock();
	// cout << "time cost: " << (double) (end - start) / CLOCKS_PER_SEC << endl;

	// vector<Zonotope<double>> R0_border = Zono_border(R0_);
    // vector<ReachableSet<double>> R_border(R0_border.size());
	// // cout << R0_border << endl;

    // for(int i = 0; i < R0_border.size(); i++){
    //     options.set_R0(R0_border[i]);
    //     mysys.reach(options, R_border[i]);
    // }

	// //cout << R.Rtime_point().set()[10][0].rs() << endl;
	// // plt::figure_size(1200, 780);
	// // for(auto i = 0; i < 1; i++){
	// // 	//cout << R.Rtime_point().set()[i][0].rs() << endl;
	// // 	Plotter::plotZonotope(R.Rtime_point().set()[i][0].rs(), 1, 2, "b");
	// // 	for(int j = 0; j < R0_border.size(); j++){
	// // 		Plotter::plotZonotope(R_border[j].Rtime_point().set()[i][0].rs(), 1, 2, "r");
	// // 	}
	// // }
	// // plt::show();
	// // cout << R.Rtime_point().set()[0][0].rs() << endl;

	// // Plotter::plotReach(R, 1, 2, "b");
	// // Plotter::plotOvertime(R, 1 , "b");
	// //第10s，第0个ReachableSetElement的Zonotope
	

	// Zonotope<double> Z_in = R.Rtime_point().set()[0][0].rs();
	// Zonotope<double> Z_bound = R_border[0].Rtime_point().set()[0][0].rs();

	// cout << Z_in << endl;
	// cout << Z_bound << endl;

	// Matrix_t<double> G1 = Z_in.generators();
	// Matrix_t<double> G2 = Z_bound.generators();
	

	// int d = G1.rows();
	// int m = G1.cols();
	// int n = G2.cols();

	// cout << d << endl;
	// cout << m << endl;
	// cout << n << endl;


	// Matrix_t<double> m1(d, m + n + 1);
	// Matrix_t<double> m2(m, m + n + 1);
	// Matrix_t<double> m3(m, m + n + 1);

	// m1 << G1, -G2, Eigen::MatrixXd::Zero(d, 1); 

	// m2 << Eigen::MatrixXd::Identity(m, m), Eigen::MatrixXd::Zero(m, n), Eigen::MatrixXd::Ones(m,1);

	// m3 << Eigen::MatrixXd::Identity(m, m), Eigen::MatrixXd::Zero(m, n), -1*Eigen::MatrixXd::Ones(m,1);

	// // cout << m1 << endl;
	// // cout << m2 << endl;
	// // cout << m3 << endl;
	// // cout << m4 << endl;

	// Matrix_t<double> M(d + 2*m, m + n + 1);

	// M << m1, m2, m3;

	// // cout << M << endl;

	// //glpk
	//  glp_prob *lp;
	// int ia[1+ 10000], ja[1+ 10000];
	// double ar[1+ 10000], z, x1, x2, x3;
	// lp = glp_create_prob();

	// glp_set_obj_dir(lp, GLP_MIN);
	// glp_add_rows(lp, d+2*m);


	// for(int i = 1; i <= d; i++){
	// 	glp_set_row_bnds(lp, i, GLP_FX,Z_in.center()(i-1) - Z_bound.center()(i-1) ,0.0);
	// }

	// for(int i = d + 1; i <= d + m; i++){
	// 	glp_set_row_bnds(lp, i, GLP_LO, 0.0, 0.0);
	// }

	// for(int i = d + m + 1; i <= d + 2*m; i++){
	// 	glp_set_row_bnds(lp, i, GLP_UP, 0.0, 0.0);
	// }

	// glp_add_cols(lp, m + n + 1);

	// for(int i = 1; i <= m; i++){
	// 	glp_set_col_bnds(lp, i, GLP_FR, 0.0, 0.0);
	// }

	// for(int i = m + 1; i <= m + n; i++){
	// 	glp_set_col_bnds(lp, i, GLP_DB, -1.0, 1.0);
	// }

	// for(int i = m + n + 1; i <= m + n + 1; i++){
	// 	glp_set_col_bnds(lp, i, GLP_DB, 0.0, 1.0);
	// }

	// glp_set_obj_coef(lp, m + n + 1, 1.0);

	// matrix2glpkForm(M,ia,ja,ar);

	// glp_load_matrix(lp, (m + n + 1)*(d+2*m), ia, ja, ar);

	// glp_simplex(lp, NULL);
	// z = glp_get_obj_val(lp);
 	// cout << z << endl;
	// glp_delete_prob(lp);

	// options.set_R0(R0_);
	// Zonotope<double> newZ = UnderApprox::getUnder(mysys, options, R0_);

	// plt::figure_size(1200, 780);
	// for(auto i = 0; i < 1; i++){
	// 	//cout << R.Rtime_point().set()[i][0].rs() << endl;
	// 	Plotter::plotZonotope(R.Rtime_point().set()[i][0].rs(), 1, 2, "b");
	// 	for(int j = 0; j < R0_border.size(); j++){
	// 		Plotter::plotZonotope(R_border[j].Rtime_point().set()[i][0].rs(), 1, 2, "r");
	// 	}
	// }
	// Plotter::plotZonotope(newZ, 1, 2, "g");
	// plt::show();

	Vector_t<double> c(3);
    c << 0, 0, 0;
    Matrix_t<double> g(3,3);
    g<< 1,0,0,
	1,1,0,
	0,0,1;

	Zonotope<double> zono(c,g);

	vector<Zonotope<double>> zbound = Zono_border(zono);

	// vector<Zonotope<double>> zbound1;

	// zbound1.push_back(zbound[0]);

	// cout << zbound1 << endl;

	vector<Zonotope<double>> zboundsplit = UnderApprox::boundSplit(zbound, 0.05);

	plt::figure_size(1200, 780);
	for(int i = 0; i < zboundsplit.size(); i++){
		 Plotter::plotZonotope(zboundsplit[i], 2, 3, "g");
	}
	plt::show();
}
