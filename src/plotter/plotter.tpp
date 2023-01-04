/**
 * @author Dejin Ren
 * @date  Jan 2022
 * @version 1.0
 *
 */

#include <plotter/plotter.h>

namespace reachSolver
{
	template <typename Number>
	static void plotZonotope(const Zonotope<Number> & zono){
		Matrix_t<Number> p = convert2polygon(zono);
		int colNum = p.cols();
		vector<Number> p_1(colNum), p_2(colNum);
		for(auto i = 0; i < colNum; i++){
			p_1[i] = p(0,i);
			p_2[i] = p(1,i);
		}
		// Set the size of output image to 1200x780 pixels
		plt::figure_size(1200, 780);
		// Plot line from given x and y data. Color is selected automatically.
		plt::plot(p_1, p_2);
	}
} // namespace reachSolver