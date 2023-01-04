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
#include <plotter/matplotlibcpp.h>
#include <plotter/plotter.h>

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

using namespace std;
using namespace reachSolver;
using namespace capd;
using namespace Eigen;
using capd::autodiff::Node;
//using namespace autodiff

#include <glpk.h>
int main()
{ glp_prob *lp;
int ia[1+1000], ja[1+1000];
double ar[1+1000], z, x1, x2, x3;
lp = glp_create_prob();

glp_set_obj_dir(lp, GLP_MAX);
glp_add_rows(lp, 3);

glp_set_row_bnds(lp, 1, GLP_UP, 0.0, 100.0);

glp_set_row_bnds(lp, 2, GLP_UP, 0.0, 600.0);

 glp_set_row_bnds(lp, 3, GLP_UP, 0.0, 300.0);
 glp_add_cols(lp, 3);

 glp_set_col_bnds(lp, 1, GLP_LO, 0.0, 0.0);
 glp_set_obj_coef(lp, 1, 10.0);

 glp_set_col_bnds(lp, 2, GLP_LO, 0.0, 0.0);
 glp_set_obj_coef(lp, 2, 6.0);

 glp_set_col_bnds(lp, 3, GLP_LO, 0.0, 0.0);
 glp_set_obj_coef(lp, 3, 4.0);
//  ia[1] = 1, ja[1] = 1, ar[1] = 1.0; /* a[1,1] = 1 */
//  ia[2] = 1, ja[2] = 2, ar[2] = 1.0; /* a[1,2] = 1 */
//  ia[3] = 1, ja[3] = 3, ar[3] = 1.0; /* a[1,3] = 1 */
//  ia[4] = 2, ja[4] = 1, ar[4] = 10.0; /* a[2,1] = 10 */
//  ia[5] = 3, ja[5] = 1, ar[5] = 2.0; /* a[3,1] = 2 */
//  ia[6] = 2, ja[6] = 2, ar[6] = 4.0; /* a[2,2] = 4 */
//  ia[7] = 3, ja[7] = 2, ar[7] = 2.0; /* a[3,2] = 2 */
//  ia[8] = 2, ja[8] = 3, ar[8] = 5.0; /* a[2,3] = 5 */
//  ia[9] = 3, ja[9] = 3, ar[9] = 6.0; /* a[3,3] = 6 */

Matrix_t<double> m(3,3);

m << 1,1,1,
    10,4,5,
    2,2,6;

matrix2glpkForm(m,ia,ja,ar);
 glp_load_matrix(lp, 9, ia, ja, ar);
 glp_simplex(lp, NULL);
 z = glp_get_obj_val(lp);
 x1 = glp_get_col_prim(lp, 1);
 x2 = glp_get_col_prim(lp, 2);
 x3 = glp_get_col_prim(lp, 3);
s36: printf("\nz = %g; x1 = %g; x2 = %g; x3 = %g\n",
z, x1, x2, x3);
 glp_delete_prob(lp);
return 0;
}

