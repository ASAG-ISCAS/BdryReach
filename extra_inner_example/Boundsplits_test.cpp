/*
 * @Author: your name
 * @Date: 2021-11-11 17:24:27
 * @LastEditTime: 2021-12-28 18:18:30
 * @LastEditors: Please set LastEditors
 * @Description: 打开koroFileHeader查看配置 进行设置:
 * https://github.com/OBKoro1/koro1FileHeader/wiki/%E9%85%8D%E7%BD%AE
 * @FilePath: /Solver/test.cpp
 */

#include <plotter/matplotlibcpp.h>
#include <plotter/plotter.h>
#include <underApprox/underApprox.h>
#include <zonotope/Zonotope.h>

namespace plt = matplotlibcpp;

using namespace std;
using namespace reachSolver;
using namespace capd;
using capd::autodiff::Node;



int main()
{
    Vector_t<double> c(3);

    Matrix_t<double> G(3,4);

    c << 0,0,0;
    G << 1, 0,1,0,
        0, 1,1, 0,
        0,0,0,0.2;

    Zonotope<double> Z(c,G);

    vector<Zonotope<double>> R0_border = Zono_border(Z);

	vector<Zonotope<double>> R0_borderSplit = UnderApprox::boundSplit(R0_border, 0.5);

    cout << "Z0 = zonotope(" << Z << ");" << endl; 
    for(int i = 0; i < R0_borderSplit.size(); i++){
        cout << "Z{" << i+1 <<"} = zonotope(" << R0_borderSplit[i] << ");" << endl; 
    }
	plt::figure_size(1200, 780);

    Plotter::plotZonotope(Z, 3, 2, "r");
    for(int i = 0; i < R0_borderSplit.size(); i++){
        Plotter::plotZonotope(R0_borderSplit[i], 3, 2, "g");
    }
	plt::show();
}
