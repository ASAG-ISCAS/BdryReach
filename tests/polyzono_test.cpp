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

int main(){
    //create R0
//     Vector_t<double> center(2);
//     center << 4, 4;
//     Matrix_t<double> generators(2,3);
//     generators<< 1,2,0,
//                  1,-2,0;
//     Matrix_t<double> Grest(2,3);
//     Grest<< 1,0,2,
//             0,0,-1;
//     Matrix_t<int> expMat(3,3);
//     expMat<< 1,0,1,
//             0,0,1,
//             1,2,1;      
    // Vector_t<double> c(3);
    // c << 4, 4, 2;
    // Matrix_t<double> G(3,4);
    // G<< 2,1,2,3,
    //              0,2,2,4,
    //              2,3,1,4;
    Vector_t<double> c(2);
    c << 2,3;
    // Matrix_t<double> G(2,6);
    // G<< 1,0,1,1,2,3,
    // 0,1,1,-1,-3,2;
        Matrix_t<double> G(2,4);
    G<< 1,0,1,1,
    0,1,1,-1;

    // Vector_t<double> c(3);
    // c << 0,0,0;
    // Matrix_t<double> G(3,4);
    // G<< 0,1,0,1,
    // 0,0,1,1,
    // 1,0,0,0;
    // G<< 1,0,0,1,
    // 0,0,1,1,
    // 0,1,0,0;

    // G<< 1,0,1,0,
    // 0,1,1,0,
    // 0,0,0,1;

    // Matrix_t<double> G(3,5);
    // G<< 1,0,-1,0,4,
    // 0,1,-1,0,3,
    // 0,0,0,1,-3;
    // Eigen::FullPivLU<Matrix_t<double>> lu(G);

    // cout << "Here is, up to permutations, its LU decomposition matrix:"
    //  << endl << lu.matrixLU() << endl;

    // // cout << "Here is the L part:" << endl;
    // // Matrix5x5 l = Matrix5x5::Identity();
    // // l.block<5,3>(0,0).triangularView<StrictlyLower>() = lu.matrixLU();
    // // cout << l << endl;
    // cout << "Here is the U part:" << endl;
    // Matrix_t<double> u = lu.matrixLU().triangularView<Eigen::Upper>();
    // cout << u << endl;
    // // cout << "Let us now reconstruct the original matrix m:" << endl;
    // // cout << lu.permutationP().inverse() * l * u * lu.permutationQ().inverse() << endl;

    // Matrix_t<double> Q = lu.permutationQ();
    // cout << Q << endl;
    // cout << G * Q << endl;


    Matrix_t<double> Gr(3,1);
    Gr<< 1,
            0,
            9;
    Matrix_t<int> exp(4,4);
    exp<< 1,0,3,1,
             0,1,1,0,
             1,0,0,0,
             2,5,0,0;
    // polyZonotope<double> polyZono(c,G,exp,Gr);
    // Zonotope<double> zono = zonotope(polyZono);
    Zonotope<double> zono(c,G);
    // zonotope(polyZono);
        Matrix_t<double> V(2,9);
    V<< -2, 2,	4,	8,	10,	6,	4,	0,	-2,
0,	0,	0,	4,	8,	8,	8,	4,	0;
    Matrix_t<double> Vtemp = V.adjoint();
        Matrix_t<double> polyPrev(6,2);
//     polyPrev<< -2,	0,
// 0,	4,
// 4,	8,
// 10,	8,
// 8,	4,
// 4,	0;
    polyPrev<< 0,	0,
1,	1,
2,	2,
3,	3,
4,	4,
1,	0;
    vector<Zonotope<double>> aa;
    // cout << zono.generators() << endl;
    cout << zono << endl;
    aa = Zono_border(zono);

    cout << Zono_border_matrix(zono)  << endl;
    cout << aa << endl;
    // cout << Zono_border_matrix(zono) << endl;
    // cout << Zono_tiling(zono) << endl;

    vector<Zonotope<double>> tiling = Zono_tiling(zono);

    // cout << tiling << endl;

    plt::figure_size(1200, 780); 

    for(int i = 0; i < tiling.size(); i++){
        Plotter::plotZonotope(tiling[i],1,2,"b");
    }

    plt::show();
    // Plotter::plotpolyZonotope(polyZono,2,1,3);
    // plt::show();
    // cout << polyshape(V);
// polyZonotope<double> pZ1,pZ2;
// cout << polyZono << endl;
// splitLongestGen(polyZono,pZ1,pZ2);
// cout << pZ1 << endl;
// cout << pZ2 << endl;
// Matrix_t<double> poly;
    // Matrix_t<double> Vabs = V.cwiseAbs();
    // cout << Vabs << endl;

    // polyZonotope<double> polyZono(c,G,exp,Gr);
    // cout << polyZono << endl;
    // polyZonotope<double> res = project(polyZono, 1);
    // cout << res << endl;
}

