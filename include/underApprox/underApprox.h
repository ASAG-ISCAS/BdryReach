/*
 * @Author: RenDejin
 * @Date: 2021-11-11 17:24:27
 * @LastEditTime: 2021-12-28 18:18:30
 * @LastEditors: Please set LastEditors
 * @Description: 打开koroFileHeader查看配置 进行设置:
 * https://github.com/OBKoro1/koro1FileHeader/wiki/%E9%85%8D%E7%BD%AE
 * @FilePath: /Solver/test.cpp
 */

#include <capd/capdlib.h>
#include <zonotope/Zonotope.h>
#include <dynamics/LinearSys.h>
#include <dynamics/NonlinearSys.h>
#include <eigen3/Eigen/Core>
#include <functional>
#include <iostream>
#include <cmath>

#include <glpk.h>

using namespace std;
using namespace capd;
using capd::autodiff::Node;

int ia[1+ 10000000], ja[1+ 10000000];
double ar[1+ 10000000];

namespace reachSolver
{
	class UnderApprox
	{
	public:
		/**
		 * @brief get the scale that translate the over-approximation into the under-approximation.
		 */
		template <typename Number>
		static double getScale(const Zonotope<Number> & Zbound, const Zonotope<Number> & Zover, double scale){

			Matrix_t<Number> G1 = Zover.generators();
			Matrix_t<Number> G2 = Zbound.generators();
			

			int d = G1.rows();
			int m = G1.cols();
			int n = G2.cols();

			// cout << d << endl;
			// cout << m << endl;
			// cout << n << endl;


			Matrix_t<Number> m1(d, m + n + 1);
			Matrix_t<Number> m2(m, m + n + 1);
			Matrix_t<Number> m3(m, m + n + 1);

			m1 << G1, -G2, Eigen::MatrixXd::Zero(d, 1); 

			m2 << Eigen::MatrixXd::Identity(m, m), Eigen::MatrixXd::Zero(m, n), Eigen::MatrixXd::Ones(m,1);

			m3 << Eigen::MatrixXd::Identity(m, m), Eigen::MatrixXd::Zero(m, n), -1*Eigen::MatrixXd::Ones(m,1);

			Matrix_t<Number> M(d + 2*m, m + n + 1);

			M << m1, m2, m3;

			// cout << M << endl;

			//glpk
			glp_prob *lp;
			int ia[1+ 100000], ja[1+ 100000];
			double ar[1+ 100000], z;

			lp = glp_create_prob();

			glp_set_obj_dir(lp, GLP_MIN);
			
			//约束
			glp_add_rows(lp, d+2*m);

			for(int i = 1; i <= d; i++){
				glp_set_row_bnds(lp, i, GLP_FX, Zbound.center()(i-1) - Zover.center()(i-1) ,0.0);
			}

			for(int i = d + 1; i <= d + m; i++){
				glp_set_row_bnds(lp, i, GLP_LO, 0.0, 0.0);
			}

			for(int i = d + m + 1; i <= d + 2*m; i++){
				glp_set_row_bnds(lp, i, GLP_UP, 0.0, 0.0);
			}

			//变量约束
			glp_add_cols(lp, m + n + 1);

			for(int i = 1; i <= m; i++){
				glp_set_col_bnds(lp, i, GLP_FR, 0.0, 0.0);
			}

			for(int i = m + 1; i <= m + n; i++){
				glp_set_col_bnds(lp, i, GLP_DB, -1.0, 1.0);
			}

			for(int i = m + n + 1; i <= m + n + 1; i++){
				glp_set_col_bnds(lp, i, GLP_DB, 0.0, scale);
			}

			//目标函数
			glp_set_obj_coef(lp, m + n + 1, 1.0);

			matrix2glpkForm(M,ia,ja,ar);

			glp_load_matrix(lp, (m + n + 1)*(d+2*m), ia, ja, ar);

			glp_simplex(lp, NULL);
			z = glp_get_obj_val(lp);

			glp_delete_prob(lp);

			return z;
        }

		/**
		 * @brief get the interval which contains the intersection between Zbound and Zover.
		 */
		template <typename Number>
		static void getScale(const Zonotope<Number> & Zbound, const Zonotope<Number> & Zover, vector<double> &al, vector<double> &au){

			Vector_t<Number> c1 = Zover.center();
			Vector_t<Number> c2 = Zbound.center();

			Matrix_t<Number> G1 = Zover.generators();
			Matrix_t<Number> G2 = Zbound.generators();

			int d = G1.rows();
			int m = G1.cols();
			int n = G2.cols();

			for(int k = 1; k <= m; k++){

				Matrix_t<Number> M(d, m + n);

				M << G1, -G2;

				// cout << M << endl;

				//glpk

				//min
				glp_prob *lp;

				glp_smcp parm;

				glp_init_smcp(&parm);

				parm.presolve = GLP_ON;

				//迭代次数限制

				parm.it_lim = 1e4;

				//重要参数！！！尽可能设置低一些，防止数值误差导致 Zover 过度收缩！！！！

				//只用调这个
				parm.tol_bnd = 1e-32;
				//Tolerance used to check if the basic solution is primal feasible. (Do not change this parameter
				//without detailed understanding its purpose.)

				// parm.tol_dj = 1e-11;
				// Tolerance used to check if the basic solution is dual feasible. (Do not change this parameter
 				// without detailed understanding its purpose.)

				// parm.tol_piv = 1e-32;
				// Tolerance used to choose eligble pivotal elements of the simplex table. (Do not change this
				// parameter without detailed understanding its purpose.)

				parm.msg_lev = GLP_MSG_OFF;

				int ia[1+ 100000], ja[1+ 100000];
				double ar[1+ 100000], min, max;

				lp = glp_create_prob();

				glp_set_obj_dir(lp, GLP_MIN);
				
				//约束
				glp_add_rows(lp, d);

				for(int i = 1; i <= d; i++){
					glp_set_row_bnds(lp, i, GLP_FX, Zbound.center()(i-1) - Zover.center()(i-1) ,0.0);
				}

				//变量约束
				glp_add_cols(lp, m + n);

				for(int i = 1; i <= m; i++){

					if(al[i-1] == au[i-1]){
						// glp_set_col_bnds(lp, i, GLP_FX, al[i-1], 0.0);
						glp_set_col_bnds(lp, i, GLP_FX, 0.0, 0.0);
					}else{
						glp_set_col_bnds(lp, i, GLP_DB, al[i-1], au[i-1]);
						// glp_set_col_bnds(lp, i, GLP_DB, -1, 1);
					}
				}

				for(int i = m + 1; i <= m + n; i++){
					glp_set_col_bnds(lp, i, GLP_DB, -1.0, 1.0);
				}

				//目标函数
				glp_set_obj_coef(lp, k, 1.0);

				matrix2glpkForm(M,ia,ja,ar);

				glp_load_matrix(lp, (m + n)*d, ia, ja, ar);

				int glpReturn = glp_simplex(lp, &parm);

				min = glp_get_obj_val(lp);

				// if(glpReturn == GLP_ENOPFS){
				// 	min = 1;
				// }

				//max
				glp_set_obj_dir(lp, GLP_MAX);

				glpReturn = glp_simplex(lp, &parm);

				max = glp_get_obj_val(lp);
				
				// cout << "原来" << endl;
				// cout << max << endl;
				// cout << min << endl<< endl;;

				//缩小 max min, 防止数值误差

				double roundloc = 1e-8;
				if(max != 1){
					if(max > 0){
						max = double(long((max + 0.5*roundloc) / roundloc)) * roundloc;
					}else{
						max = - double(long((abs(max) - 0.5*roundloc) / roundloc)) * roundloc;
					}
				}
				
				if(min != -1){
					if(min > 0){
						min = double(long((min - 0.5*roundloc) / roundloc)) * roundloc;
					}else{
						min = - double(long((abs(min) + 0.5*roundloc) / roundloc)) * roundloc;
					}
				}
				
				if(max > 1){
					max = 1;
				}

				if(min < -1){
					min = -1;
				}

				// cout << "现在" << endl;
				// cout << max << endl;
				// cout << min << endl;

				if(glpReturn != GLP_ENOPFS){

					intervalDiv(al, au, k - 1, min, max);
				}else{

					break;
				}

				// cout << al << endl;
				// cout << au << endl;

				glp_delete_prob(lp);
			}
		}

		/**
		 * @brief get the interval which contains the intersection between Zbound and Zover.
		 * 计算出全部相交的区间，不再在程序中去除相交区间
		 * 如果相交返回 true, 不相交返回 false.
		 */
		template <typename Number>
		static int getScale2(const Zonotope<Number> & Zbound, const Zonotope<Number> & Zover, vector<double> &al, vector<double> &au){

			Vector_t<Number> c1 = Zover.center();
			Vector_t<Number> c2 = Zbound.center();

			Matrix_t<Number> G1 = Zover.generators();
			Matrix_t<Number> G2 = Zbound.generators();

			int d = G1.rows();
			int m = G1.cols();
			int n = G2.cols();

			for(int k = 1; k <= m; k++){

				Matrix_t<Number> M(d, m + n);

				M << G1, -G2;

				// cout << M << endl;

				//glpk

				//min
				glp_prob *lp;

				glp_smcp parm;

				glp_init_smcp(&parm);

				parm.presolve = GLP_ON;

				//重要参数！！！尽可能设置低一些，防止数值误差导致 Zover 过度收缩！！！！

				//只用调这个
				parm.tol_bnd = 1e-12;
				//Tolerance used to check if the basic solution is primal feasible. (Do not change this parameter
				//without detailed understanding its purpose.)

				// parm.tol_dj = 1e-11;
				// Tolerance used to check if the basic solution is dual feasible. (Do not change this parameter
 				// without detailed understanding its purpose.)

				// parm.tol_piv = 1e-32;
				// Tolerance used to choose eligble pivotal elements of the simplex table. (Do not change this
				// parameter without detailed understanding its purpose.)

				parm.msg_lev = GLP_MSG_OFF;

				int ia[1+ 100000], ja[1+ 100000];
				double ar[1+ 100000], min, max;

				lp = glp_create_prob();

				glp_set_obj_dir(lp, GLP_MIN);
				
				//约束
				glp_add_rows(lp, d);

				for(int i = 1; i <= d; i++){
					glp_set_row_bnds(lp, i, GLP_FX, Zbound.center()(i-1) - Zover.center()(i-1) ,0.0);
				}

				//变量约束
				glp_add_cols(lp, m + n);

				for(int i = 1; i <= m; i++){

					if(al[i-1] == au[i-1]){
						glp_set_col_bnds(lp, i, GLP_FX, al[i-1], 0.0);
					}else{
						glp_set_col_bnds(lp, i, GLP_DB, al[i-1], au[i-1]);
						// glp_set_col_bnds(lp, i, GLP_DB, -1, 1);
					}
				}

				for(int i = m + 1; i <= m + n; i++){
					glp_set_col_bnds(lp, i, GLP_DB, -1.0, 1.0);
				}

				//目标函数
				glp_set_obj_coef(lp, k, 1.0);

				matrix2glpkForm(M,ia,ja,ar);

				glp_load_matrix(lp, (m + n)*d, ia, ja, ar);

				int glpReturn = glp_simplex(lp, &parm);

				min = glp_get_obj_val(lp);

				//如果无解，说明不相交
				if(glpReturn == GLP_ENOPFS){
					return 0;
				}

				//max
				glp_set_obj_dir(lp, GLP_MAX);

				glpReturn = glp_simplex(lp, &parm);

				max = glp_get_obj_val(lp);
				
				al[k-1] = min;
				au[k-1] = max;

				// cout << al << endl;
				// cout << au << endl;


				glp_delete_prob(lp);
			}

			return 1;
		}
		
		/**
		 * @brief get the interval that Zover does not intersect with all Zbounds.
		 */
		static void intervalDiv(vector<double>& al, vector<double>& au, int i, double min, double max){

			double alnew, aunew;

			if((min - al[i]) >= (au[i] - max)){

				alnew = al[i];
				aunew = min;
			}else{
				
				alnew = max;
				aunew = au[i];
			}

			al[i] = alnew;
			au[i] = aunew;
		}

		/**
		 * @brief get the interval that Zover does not intersect with all Zbounds.
		 * 使用 Cora under 里面 interval divide 方法。
		 */
		static void intervalDiv(vector<double>& al, vector<double>& au, const vector<double>& alold, const vector<double>& auold, const vector<double>& norm2){

			int idx = 0;
			double lengthmax = 0;
			int directmax;
			double newal, newau;
			for(int i = 0; i < al.size(); i++){

				double intervalmax;
				int direct;

				if((auold[i] - au[i]) > (al[i] - alold[i])){
					
					// intervalmax = (auold[i] - au[i]) * norm2[i];
					intervalmax = (auold[i] - au[i]);
					direct = 1;
					
				}else{

					// intervalmax = (al[i] - alold[i]) * norm2[i];
					intervalmax = (al[i] - alold[i]);
					direct = -1;

				}

				if(intervalmax > lengthmax){

					lengthmax = intervalmax;
					idx = i;

					if(direct == 1){
						
						newal = au[idx];
						newau = auold[idx];

					}else{

						newal = alold[idx];
						newau = al[idx];

					}
				}
			}

			cout << "in scale" << endl;
			cout << alold << endl;
			cout << auold << endl;
			cout << al << endl;
			cout << au << endl;
			cout << au[idx] << " " << auold[idx] << endl;
			cout << alold[idx] << " " << al[idx] << endl;
			cout << idx << endl;

			al = alold;
			au = auold;

			al[idx] = newal;
			au[idx] = newau;

		}

		/**
		 * @brief get the under-approximation of the reachable set.
		 */
		template <typename Number>
		static Zonotope<Number> getUnder(NonlinearSys<Number> mysys, NonlinearSys<Number> mysysBack, ReachOptions<Number> options, Zonotope<Number> R0, double radius, double over_step, double bound_step){

				ReachableSet<Number> R;

				options.set_R0(R0);

				// options.set_tFinal(3);

				options.set_time_step(over_step);

				mysys.reach(options, R);

				vector<Zonotope<Number>> R0_border = Zono_border(R0);

				vector<Zonotope<Number>> R0_borderSplit = boundSplit(R0_border, radius);

				cout << R0_borderSplit.size() << endl;

				// cout << R0_borderSplit[0] << endl;

				vector<ReachableSet<Number>> R_border(R0_borderSplit.size());

				options.set_time_step(bound_step);

				for(int i = 0; i < R0_borderSplit.size(); i++){
					options.set_R0(R0_borderSplit[i]);
					mysys.reach(options, R_border[i]);
				}

				// cout << R_border[0].Rtime_point().set()[R_border[0].Rtime_point().set().size() - 1][0].rs() << endl;

				double scale = 1.0;

				Zonotope<Number> Zover = R.Rtime_point().set()[R.Rtime_point().set().size() - 1][0].rs();

				Plotter::plotZonotope(Zover, 1, 2, "b");

				for(int i = 0; i < R_border.size(); i++){

					Zonotope<Number> Zbound = R_border[i].Rtime_point().set()[R_border[i].Rtime_point().set().size() - 1][0].rs();

					Plotter::plotZonotope(Zbound, 1, 2, "r");

					scale = getScale(Zbound, Zover, scale);

					cout << scale << endl;
				}

			Zonotope<Number> res(Zover.center(), scale * Zover.generators());

			Zonotope<Number> ZunderC(Zover.center(),Eigen::MatrixXd::Zero(Zover.center().rows(), 1));

			//如果验证成功，返回结果，否则程序停止
			if(verificationUnder(mysysBack, options, ZunderC, R0, over_step)){

				return res;

			}else{
				
				cerr<<"can not do the under process.\n";
        		exit(1);
				
			}

		}
		
		/**
		 * @brief get the under-approximation of the reachable set.(Use Contrator LP)
		 */
		template <typename Number>
		static Zonotope<Number> getUnderClp(NonlinearSys<Number> mysys, NonlinearSys<Number> mysysBack, ReachOptions<Number> options, Zonotope<Number> R0, double radius, double over_step, double bound_step, int Zover_order, int Zbound_order){

			ReachableSet<Number> R;

			options.set_R0(R0);

			options.set_time_step(over_step);

			mysys.reach(options, R);

			vector<Zonotope<Number>> R0_border = Zono_border(R0);

			vector<Zonotope<Number>> R0_borderSplit = boundSplit(R0_border, radius);

			cout << R0_borderSplit.size() << endl;

			vector<ReachableSet<Number>> R_border(R0_borderSplit.size());

			options.set_time_step(bound_step);
			options.set_zonotope_order(Zbound_order);

			for(int i = 0; i < R0_borderSplit.size(); i++){
				options.set_R0(R0_borderSplit[i]);
				mysys.reach(options, R_border[i]);
			}

			Zonotope<Number> Zover = R.Rtime_point().set()[R.Rtime_point().set().size() - 1][0].rs();

			// Zover.Reduce(Zover_order);

			int dim = R0.center().rows();

			// //用各个 Zbound 的 generator 和 center 的平均值构建“凸包”。

			// Zonotope<Number> Zover;

			// int dim = R0.center().rows();

			// // cout << Zover << endl;

			// int Zovercols = 0;

			// Matrix_t<Number> addBoundG(dim, Zovercols);
			// Vector_t<Number> Zboundcentermean;
			// Zboundcentermean.setZero(dim);
			
			// // addBoundG = Zover.generators();

			// for(int i = 0; i < R_border.size(); i++){

			// 	Zonotope<Number> Zbound = R_border[i].Rtime_point().set()[R_border[i].Rtime_point().set().size() - 1][0].rs();

			// 	// Plotter::plotZonotope(Zbound, 1, 2, "b");

			// 	Zovercols += Zbound.generators().cols();
		
			// 	Zboundcentermean = Zboundcentermean + Zbound.center();

			// 	Matrix_t<Number> addtemp = addBoundG;

			// 	addBoundG.conservativeResize(dim, Zovercols);

			// 	addBoundG << addtemp, Zbound.generators();
			// }

			// cout << Zboundcentermean << endl;

			// Zover.set_generators(addBoundG*0.5);
			// Zover.set_center(Zboundcentermean/R_border.size());

			// Plotter::plotZonotope(Zover, 1, 2, "r");

			// Zover.Reduce(200);

			// cout << Zover << endl;

			// // Plotter::plotZonotope(Zover, 1, 2, "r");

			// // plt::show();

			vector<double> al(Zover.generators().cols(), -1);
			vector<double> au(Zover.generators().cols(), 1);

			Plotter::plotZonotope(Zover, 1, 2, "b");

			for(int i = 0; i < R_border.size() /*R_border.size()*/; i++){
			// for(int i = R_border.size() - 1; i >= R_border.size() - 1 /*R_border.size()*/; i--){
				
				// cout << "compute " << i << " th boundary" << endl; 
				
				Zonotope<Number> Zbound = R_border[i].Rtime_point().set()[R_border[i].Rtime_point().set().size() - 1][0].rs();

				// cout << Zbound.generators().cols() << endl;
				// Plotter::plotZonotope(Zbound, 1, 2, "m");

				//计算 Zbound 前 dim-1 大的 generator(dim 为维数)
				Matrix_t<Number> ZboundMaxg = sortGenDes(Zbound.generators());
			
				Matrix_t<Number> maxg = ZboundMaxg.leftCols(dim - 1);

				//计算 maxg 所代表面的法向量
				Vector_t<Number> normalVec(dim);
				Matrix_t<Number> tempB(dim, dim - 1);

				for(int j = 0; j < dim; j++){
					tempB = maxg;
					RemoveRow(tempB, j);
					normalVec(j) = pow(-1, j) * tempB.determinant();
				}

				//按 Zover generator 与法向量夹角从小到大排序(cos 值从大到小)
				vector<int> idx(Zover.generators().cols());

				Matrix_t<Number> newG = sortGenCos(Zover.generators(), normalVec, idx);

				// Matrix_t<Number> newGDes = sortGenDes(newG);

				Zover.set_generators(newG);

				// cout << maxg << endl;
				// cout << newG << endl;

				vector<double> altemp(al.size()), autemp(au.size());

				for(int k = 0; k < al.size(); k++){

					altemp[k] = al[idx[k]];
					autemp[k] = au[idx[k]];
				}

				al = altemp;
				au = autemp;

				// cout << al.size() << endl;

				// cout << "begin scale" << endl;
				getScale(Zbound, Zover, al, au);
				// cout << "end scale" << endl;
				// cout << Zover << endl;
				// cout << Zbound << endl;

				// //glpk 解出来有可能超过[-1,1]

				// for(int k = 0; k < al.size(); k++){
				// 	if(al[k] < -1){
				// 		al[k] = -1;
				// 	}

				// 	if(au[k] > 1){
				// 		au[k] = 1;
				// 	}
				// }

				// cout << al << endl;
				// cout << au << endl;

				Vector_t<double> Cadd(al.size());
				Vector_t<double> Gsacle(al.size());

				for(int k = 0; k < al.size(); k++){
					
					if(abs(au[k] - al[k]) > 1e-12){
						Cadd(k) = 0.5 * (au[k] + al[k]);
						Gsacle(k) = 0.5 * (au[k] - al[k]);
					}else{
						Cadd(k) = 0;
						Gsacle(k) = 0;
					}
					// Cadd(k) = 0.5 * (au[k] + al[k]);

					// Gsacle(k) = 0.5 * (au[k] - al[k]);
				}

				Matrix_t<Number> genMatrix = Zover.generators() * Gsacle.asDiagonal();

				deleteZeros(genMatrix);

				Zonotope<Number> res(Zover.center() + Zover.generators() * Cadd, genMatrix);

				for(int k = 0; k < al.size(); k++){
					if(abs(al[k] - au[k]) < 1e-18){
						al.erase(al.begin() + k);
						au.erase(au.begin() + k);
						k--;
					}
				}

				for(int k = 0; k < al.size(); k++){
					al[k] = -1;
					au[k] = 1;
				}

				Zover = res;
			}
			
			cout << Zover << endl;

			Zonotope<Number> ZunderC(Zover.center(),Eigen::MatrixXd::Zero(Zover.center().rows(), 1));

			// return Zover;

			//如果验证成功，返回结果，否则调用整体缩减的 getUnder 方法
			if(verificationUnder(mysysBack, options, ZunderC, R0, over_step)){

				cout << "Verification success!!!" << endl;
				return Zover;

			}else{
				
				cout << "getUnderClp doesn't work!!!" << endl;
				Zover = getUnder(mysys, mysysBack, options, R0, radius, over_step, bound_step);

				return Zover;

			}

			return Zover;
			
			// Zonotope<Number> Zbound = R_border[4].Rtime_point().set()[R_border[4].Rtime_point().set().size() - 1][0].rs();

			// Plotter::plotZonotope(Zbound, 1, 2, "y");

		}

		/**
		 * @brief get the under-approximation of the reachable set.(Use Contrator LP)
		 * 不再求一个 generator 的范围，而是求所有的范围再除去
		 */
		template <typename Number>
		static Zonotope<Number> getUnderClp2(NonlinearSys<Number> mysys, NonlinearSys<Number> mysysBack, ReachOptions<Number> options, Zonotope<Number> R0, double radius, double over_step, double bound_step){

			ReachableSet<Number> R;

			options.set_R0(R0);

			options.set_time_step(over_step);

			mysys.reach(options, R);

			vector<Zonotope<Number>> R0_border = Zono_border(R0);

			vector<Zonotope<Number>> R0_borderSplit = boundSplit(R0_border, radius);

			cout << R0_borderSplit.size() << endl;

			vector<ReachableSet<Number>> R_border(R0_borderSplit.size());

			options.set_time_step(bound_step);

			for(int i = 0; i < R0_borderSplit.size(); i++){
				options.set_R0(R0_borderSplit[i]);
				mysys.reach(options, R_border[i]);
			}

			Zonotope<Number> Zover = R.Rtime_point().set()[R.Rtime_point().set().size() - 1][0].rs();

			int dim = R0.center().rows();

			Matrix_t<Number> ZoverG = Zover.generators();

			vector<double> al(Zover.generators().cols(), -1);
			vector<double> au(Zover.generators().cols(), 1);

			Plotter::plotZonotope(Zover, 1, 2, "b");

			vector<double> norm2(ZoverG.cols(), 0);

			for(int i = 0; i < norm2.size(); i++){
				for(int j = 0; j < ZoverG.rows(); j++){
					norm2[i] += ZoverG(j,i) * ZoverG(j,i);
				}
				norm2[i] = sqrt(norm2[i]);
			}

			for(int i = 0; i < R_border.size() /*R_border.size()  48*/; i++){

				Zonotope<Number> Zbound = R_border[i].Rtime_point().set()[R_border[i].Rtime_point().set().size() - 1][0].rs();				
				
				cout << Zover << endl;
				cout << Zbound << endl;

				Matrix_t<Number> ZoverG = Zover.generators();

				Plotter::plotZonotope(Zbound, 1, 2, "r");

				vector<double> alold = al;
				vector<double> auold = au;

				int ifInsect = getScale2(Zbound, Zover, al, au);

				cout << "getScale2" << endl;
				cout << al << endl;
				cout << au << endl;
				
				cout << "ifInsect " << ifInsect << endl;
				if(ifInsect){
					//检查是否全部相交
					int insectAll = 1;

					for(int k = 0; k < al.size(); k++){

						if(abs(al[k] - alold[k]) > 1e-12){
							insectAll = 0;
							break;
						}

						if(abs(au[k] - auold[k]) > 1e-12){
							insectAll = 0;
							break;
						}
					}

					if(insectAll == 1){
						
						getScale(Zbound, Zover, al, au);
					}else{

						intervalDiv(al, au, alold, auold, norm2);
					}
				}

				cout << "intervalDiv" << endl;
				cout << al << endl;
				cout << au << endl;
			}

			Vector_t<double> Cadd(al.size());
			Vector_t<double> Gsacle(al.size());

			for(int k = 0; k < al.size(); k++){
				
				if(abs(au[k] - al[k]) > 1e-12){
					Cadd(k) = 0.5 * (au[k] + al[k]);
				}
				
				Gsacle(k) = 0.5 * (au[k] - al[k]);
			}

			Matrix_t<Number> genMatrix = Zover.generators() * Gsacle.asDiagonal();

			deleteZeros(genMatrix);

			Zonotope<Number> res(Zover.center() + Zover.generators() * Cadd, genMatrix);

			// cout << res << endl;

			return res;

			// for(int k = 0; k < al.size(); k++){
			// 	if(abs(al[k] - au[k]) < 1e-18){
			// 		al.erase(al.begin() + k);
			// 		au.erase(au.begin() + k);
			// 		k--;
			// 	}
			// }

			// for(int k = 0; k < al.size(); k++){
			// 	al[k] = -1;
			// 	au[k] = 1;
			// }

			// Zover = res;
			
			// cout << Zover << endl;

			// Zonotope<Number> ZunderC(Zover.center(),Eigen::MatrixXd::Zero(Zover.center().rows(), 1));

			// //如果验证成功，返回结果，否则调用整体缩减的 getUnder 方法
			// if(verificationUnder(mysysBack, options, ZunderC, R0, over_step)){

			// 	return Zover;

			// }else{
				
			// 	Zover = getUnder(mysys, mysysBack, options, R0, radius, over_step, bound_step);

			// 	return Zover;

			// }
		}

		/**
		 * @brief get the series of under-approximation of the reachable set.
		 */
		template <typename Number>
		static vector<Zonotope<Number>> underReach(NonlinearSys<Number> mysys, NonlinearSys<Number> mysysBack, ReachOptions<Number> options, Zonotope<Number> R0, double overRtime, int steps, double radius, double over_step, double bound_step){
			
			options.set_tFinal(overRtime);

			vector<Zonotope<Number>> res;

			res.push_back(R0);

			for(int i = 0; i < steps; i++){
				res.push_back(getUnder(mysys, mysysBack, options, res[i], radius, over_step, bound_step));
			}

			return res;
		}

		/**
		 * @brief get the series of under-approximation of the reachable set.(Use Clp)
		 */
		template <typename Number>
		static vector<Zonotope<Number>> underReachClp(NonlinearSys<Number> mysys, NonlinearSys<Number> mysysBack, ReachOptions<Number> options, Zonotope<Number> R0, double overRtime, int steps, double radius, double over_step, double bound_step, int Zover_order, int Zbound_order){
			
			options.set_zonotope_order(Zover_order);

			options.set_tFinal(overRtime);

			vector<Zonotope<Number>> res;

			res.push_back(R0);

			for(int i = 0; i < steps; i++){
				res.push_back(getUnderClp(mysys, mysysBack, options, res[i], radius, over_step, bound_step, Zover_order, Zbound_order));

				cout << "迭代次数 " << i+1 << endl; 
			}

			return res;
		}

		
		/**
		 * @brief Split the bound such that radius of zonotope smaller than desired radius.
		 */
		template <typename Number>
		static vector<Zonotope<Number>> boundSplit(vector<Zonotope<Number>> bound, double radius){

			double dSquare = radius * radius;

			vector<Zonotope<Number>> res;

			for(int i = 0; i < bound.size(); i++){
				
				// int flag = 1;
				vector<Vector_t<Number>> allCenter;
				vector<Vector_t<Number>> allCentertemp;

				Matrix_t<Number> G = bound[i].generators();
				Matrix_t<Number> newG = bound[i].generators();
				Vector_t<Number> center = bound[i].center();

				allCentertemp.push_back(center);

				for(int j = 0; j < G.cols(); j++){
					
					Vector_t<Number> g = G.col(j);

					double normSquare = vector2NormSquare(g);

					if(normSquare > dSquare){

						// flag = 0;

						int n = ceil(sqrt(normSquare/dSquare));

						Vector_t<Number> newg = g/n;
						newG.col(j) = newg;

						for(int m = 0; m < allCentertemp.size(); m++){
							
							Vector_t<Number> newcenterAdd = - g + newg;

							for(int k = 0; k < n; k++){
								
								allCenter.push_back(allCentertemp[m] + newcenterAdd);

								newcenterAdd += 2*newg;
							}
						}

						allCentertemp = allCenter;
						allCenter.clear();

						// cout << allCentertemp << endl;
					}
				}

				// if(flag == 1){
				// 	allCenter.push_back(center);
				// }

				// cout << allCentertemp << endl;

				for(int j = 0; j < allCentertemp.size(); j++){

					res.push_back(Zonotope<Number>(allCentertemp[j], newG));
				}
			}

			return res;
		}

		/**
		 * @brief determines if Z1 is contained in Z2.(sufficient condition)
		 */
		template <typename Number>
		static bool zonotopeIsIn(const Zonotope<Number>& Z1, const Zonotope<Number>& Z2){
			
			int d = Z1.center().rows();

			Vector_t<Number> c1 = Z1.center();
			Vector_t<Number> c2 = Z2.center();

			Matrix_t<Number> G1 = Z1.generators();
			Matrix_t<Number> G2 = Z2.generators();

			int m = G1.cols();
			int n = G2.cols();

			Matrix_t<Number> M1;
			Matrix_t<Number> M2;
			Matrix_t<Number> M3;
			Matrix_t<Number> M(m*d+n+d, 2*n*(m+1));

			M1.setZero(m*d, 2*n*(m+1));

			for(int i = 0; i < d; i++){
				for(int j = 0; j < m; j++){
					for(int l = 0; l < n; l++){

						M1(i*m+j,j*n+l) = G2(i,l);
					}

					for(int l = 0; l < n; l++){

						M1(i*m+j,n*(m+1)+j*n+l) = -G2(i,l);
					}					
				}
			}

			// cout << M1 << endl;

			M2 = Eigen::MatrixXd::Identity(n, n);

			Matrix_t<Number> M2temp;
			
			for(int i = 0; i < 2*(m + 1) - 1; i++){

				M2temp = M2;

				M2.conservativeResize(n, (i+2)*n);

				M2 << M2temp, Eigen::MatrixXd::Identity(n, n);
			}

			// cout << M2 << endl;
			
			M3.setZero(d, 2*n*(m+1));

			for(int i = 0; i < d; i++){

				for(int j = 0; j < n; j++){

					M3(i, m*n+j) = G2(i,j);
					M3(i, (m+1)*n + m*n + j) = -G2(i,j);
				}
			}
			M << M1, M2, M3;

			// cout << M << endl;

			//glpk

			glp_prob *lp;

			glp_smcp parm;

			glp_init_smcp(&parm);

			parm.presolve = GLP_ON;

			parm.msg_lev = GLP_MSG_OFF;

			//只用调这个
			parm.tol_bnd = 1e-14;
			// //Tolerance used to check if the basic solution is primal feasible. (Do not change this parameter
			// //without detailed understanding its purpose.)

			// parm.tol_dj = 1e-32;
			// // Tolerance used to check if the basic solution is dual feasible. (Do not change this parameter
			// // without detailed understanding its purpose.)

			// parm.tol_piv = 1e-32;
			// // Tolerance used to choose eligble pivotal elements of the simplex table. (Do not change this
			// // parameter without detailed understanding its purpose.)

			// int ia[1+ 100000], ja[1+ 100000];
			// double ar[1+ 100000];

			lp = glp_create_prob();

			glp_set_obj_dir(lp, GLP_MIN);
			
			//约束
			glp_add_rows(lp, m*d+n+d);

			for(int i = 1; i <= d; i++){

				for(int j = 1; j<= m; j++){

					glp_set_row_bnds(lp, (i-1)*m+j, GLP_FX, G1(i-1, j-1) ,0.0);		
				}
			}

			for(int i = m*d + 1; i <= m*d + n; i++){

				glp_set_row_bnds(lp, i, GLP_UP, 0, 1);
			}

			for(int i = m*d + n + 1; i <= m*d + n + d; i++){

				glp_set_row_bnds(lp, i, GLP_FX, c2(i-(m*d + n + 1)) - c1(i-(m*d + n + 1)), 0);
			}

			//变量约束
			glp_add_cols(lp, 2*n*(m+1));

			for(int i = 1; i <= 2*n*(m+1); i++){

				glp_set_col_bnds(lp, i, GLP_LO, 0, 0);

			}

			//目标函数(无目标函数)
			// glp_set_obj_coef(lp, 1, 1.0);

			matrix2glpkForm(M,ia,ja,ar);

			glp_load_matrix(lp, (m*d+n+d)*2*n*(m+1), ia, ja, ar);

			int glpReturn = glp_simplex(lp, &parm);

			Vector_t<Number> ans(2*n*(m+1));

			for(int i = 1; i <= 2*n*(m+1); i++){

				ans[i-1] = glp_get_col_prim(lp, i);
				// cout << glp_get_col_prim(lp, i) << endl;
			}

			// Vector_t<Number> realans(n*(m+1));
			
			Matrix_t<Number> Gamma(n,m);
			Vector_t<Number> Beda(n);
			for(int i = 0; i < n*m; i++){

				// cout << i/m << " " << i%m << endl;
				Gamma(i%n, i/n) = ans[i] - ans[i + n*(m+1)];
			}

			for(int i = 0; i < n; i++){

				Beda(i) = ans[i+n*m] - ans[i + n*(m+1)+n*m];
			}

			// cout << Gamma << endl;
			// cout << Beda << endl;

			// for(int i = 0; i < n*(m+1); i++){

			// 	cout << realans[i] << endl;
			// 	if((i+1)%6 == 0){
			// 		cout << endl;
			// 	}
			// }

			if(glpReturn != GLP_ENOPFS){

				return true;
			}else{

				return false;
			}
		}

		/**
		 * @brief verify whether the "under-approximation" is really the under.
		 */
		template <typename Number>
		static bool verificationUnder(NonlinearSys<Number> mysysBack, ReachOptions<Number> options, Zonotope<Number> Zunder, Zonotope<Number> R0, double over_step){

			ReachableSet<Number> Rback;

			options.set_R0(Zunder);

			// cout << options.tFinal() << endl;
			// cout << options.time_step() << endl;

			mysysBack.reach(options, Rback);

			Zonotope<Number> ZunderR0 = Rback.Rtime_point().set()[Rback.Rtime_point().set().size() - 1][0].rs();

			// cout << ZunderR0 << endl;

			ZunderR0.Reduce(1);

			// cout << Zunder << endl;
			// cout << ZunderR0 << endl;
			// cout << R0 << endl;

			// Plotter::plotZonotope(Zunder, 1, 2, "k","4");
			// Plotter::plotZonotope(ZunderR0, 1, 2, "k","4");
			// Plotter::plotZonotope(R0, 1, 2, "c","4");

			return zonotopeIsIn(ZunderR0, R0);
		}
	};
} // namespace reachSolver
