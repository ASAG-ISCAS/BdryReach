#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point.hpp>
#include <boost/geometry/geometries/polygon.hpp>
#include <cmath>
#include <vector>
#include <iostream>

#include <zonotope/polyZonotope.h>
#include <zonotope/Zonotope.h>

#include <eigen3/Eigen/Core>

#include <plotter/matplotlibcpp.h>
#include <plotter/plotter.h>

using namespace std;
using namespace reachSolver;
using namespace capd;
using capd::autodiff::Node;

namespace bg = boost::geometry;
namespace bgi = boost::geometry::index;

typedef bg::model::point<double, 2, bg::cs::cartesian> point;
typedef bg::model::polygon<point, true, true> polygon; // ccw, open polygon

int main()
{ 
    polygon p1, p2, p;
    Matrix_t<double> polyAll(6,2);
//     polyAll <<     5.2500,    7.2500,
//     8.0000,    8.0000,
//    10.0000,    8.0000,
//     8.7500,    5.7500,
//     7.6250,    4.6250,
//     6.3750,   2.3750,
//     4.0000,         0,
//     2.0000,         0,
//     0.7500,    0.7500,
//     0.8333,    0.9000,
//    -0.3750,    1.6250,
//    -0.1071,    2.2500,
//    -0.7500,    2.2500,
//          0,    4.0000,
//     2.3750,   6.3750,
//     5.1250,    7.1250;

polyAll <<0, 0, 
4 ,0 ,
7.5, 3.5, 
10, 8 ,
6, 8 ,
2.5, 4.5 ;
    Matrix_t<double> poly(6,2);
//     poly << -0.7500,    2.2500,
//          0 ,   4.0000,
//     2.3750,    6.3750,
//     5.1250 ,   7.1250,
//     7.1250,    7.1250,
//     6.3750 ,   5.3750,
//     6.0000 ,   5.0000,
//     5.2500 ,   3.2500,
//     2.8750 ,   0.8750,
//     0.8750 ,   0.8750,
//    -0.3750 ,   1.6250,
//    -0.1071 ,  2.2500;

poly << -1.5, 0.5 ,
2.5, 0.5, 
6, 4, 
7.5, 7.5, 
3.5, 7.5, 
0, 4 ;
    // for(int i = 0; i < polyAll.rows(); i++){
    //     p1.outer().push_back(point(polyAll(i,0),polyAll(i,1)));
    // }
    // for(int i = 0; i < poly.rows(); i++){
    //     p2.outer().push_back(point(poly(i,0),poly(i,1)));
    // }

    // std::vector<polygon> output;
    // boost::geometry::union_(p1,p2,output);

    // Matrix_t<double> res(output[0].outer().size() - 1,2);
    // for(int i = 0; i < output[0].outer().size() - 1; i++){
    //     res(i,0) = output[0].outer()[i].get<0>();
    //     res(i,1) = output[0].outer()[i].get<1>();
    // }
    // auto = output[0].outer();
    // cout << output.size() << endl;

    Matrix_t<double> res;
    Plotter::unionPoly(polyAll,poly,res);
    cout << res << endl;
}