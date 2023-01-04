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


int main()
{
    
    Vector_t<double> center(3);
    center << 0.658644099322084,
0.325810579325642,
0.442178212373429;
    Matrix_t<double> generators(3,5);
    generators<< -0.0468997642128547,	-0.0378868401340827,	0.204529313418159,	0,	0,
0.109344637273894,	0.00304141877654800, 0,	0.0840386893938892,	0,
-0.0923083994650241,	-0.0840387697281924,	0,	0,	0.106055621969682;

   Zonotope<double> Z(center, generators);

   vector<Zonotope<double>> R0_border = Zono_border(Z);

	cout << R0_border.size() << endl;
   cout << R0_border << endl;
}
