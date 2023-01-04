#include <plotter/matplotlibcpp.h>
#include <plotter/plotter.h>
//#include <cmath>
#include <zonotope/Zonotope.h>
#include <eigen3/Eigen/Core>

namespace plt = matplotlibcpp;
using namespace std;
using namespace reachSolver;

int main() {
  std::vector<double> x = {1, 2, 3, 4};
  std::vector<double> y = {1, 4, 9, 16};

  plt::plot(x, y, {{"label", "f(x)"}});  // add the label f(x)
  plt::legend(); 
  plt::show();

  return 0;
}
// int main()
// {
//     Vector_t<double> c(2);
// 	Matrix_t<double> matrix(2,7);
// 	c << -4.5, 8.44;
// 	matrix << 2.523, 0, 0, 4.35, 2, 2, 1,
// 		 0, 0, 0, -1, 2, 0, 2;
// 	Zonotope<double> Z(c,matrix);
//     Plotter::plotZonotope(Z);
// 	/*
//     Matrix_t<double> p = convert2polygon(Z);
// 	cout << p;
//     int colNum = p.cols();
//     vector<double> p_1(colNum), p_2(colNum);
//     for(auto i = 0; i < colNum; i++){
//         p_1[i] = p(0,i);
//         p_2[i] = p(1,i);
//     }

//     // Set the size of output image to 1200x780 pixels
//     plt::figure_size(1200, 780);
//     // Plot line from given x and y data. Color is selected automatically.
//     plt::plot(p_1, p_2);

//     for(auto i = 0; i < colNum; i++){
//         p_1[i] = p(0,i)+1;
//         p_2[i] = p(1,i)+1;
//     }

//     plt::plot(p_1, p_2);
//     /*
//     // Plot a red dashed line from given x and y data.
//     plt::plot(x, w,"r--");
//     // Plot a line whose name will show up as "log(x)" in the legend.
//     plt::named_plot("log(x)", x, z);
//     // Set x-axis to interval [0,1000000]
//     plt::xlim(0, 1000*1000);
//     // Add graph title
//     plt::title("Sample figure");
//     // Enable legend.
//     plt::legend();
//     // Save the image (file format is determined by the extension)
//     plt::save("./basic.png"); 
//     */
//     plt::show();
   
// }
