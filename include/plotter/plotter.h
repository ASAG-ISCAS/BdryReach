/**
 * @file   plotter.h
 * @brief  Zonotope plotter
 * @author Dejin Ren
 * @date  Jan 2022
 * @version 1.0
 *
 */

#ifndef PLOTTER_H
#define PLOTTER_H

#include <zonotope/Zonotope.h>
#include <zonotope/polyZonotope.h>
#include <eigen3/Eigen/Core>
#include <plotter/matplotlibcpp.h>
#include <zonotope/globalfunctions.h>
#include <dynamics/NonlinearSys.h>

#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point.hpp>
#include <boost/geometry/geometries/polygon.hpp>

namespace bg = boost::geometry;
namespace bgi = boost::geometry::index;

namespace plt = matplotlibcpp;
using namespace std;

namespace reachSolver
{
	class Plotter
	{
	public:
		/**
		 * @brief plot the Zonotope
		 */
		template <typename Number>
		static void plotZonotope(const Zonotope<Number> & zono, int dim1, int dim2, string color){
			Vector_t<Number> c(2);
			Matrix_t<Number> G(2, zono.generators().cols());

			c(0) = zono.center()(dim1 - 1);
			c(1) = zono.center()(dim2 - 1);

			G.row(0) = zono.generators().row(dim1 - 1);
			G.row(1) = zono.generators().row(dim2 - 1);

			Zonotope<Number> Zono = Zonotope<Number>(c, G);

			Matrix_t<Number> p = convert2polygon(Zono);
			int colNum = p.cols();
			vector<Number> p_1(colNum), p_2(colNum);
			for(auto i = 0; i < colNum; i++){
				p_1[i] = p(0,i);
				p_2[i] = p(1,i);
			}
			// Set the size of output image to 1200x780 pixels
			//plt::figure_size(1200, 780);
			// Plot line from given x and y data. Color is selected automatically.
			// plt::plot(p_1, p_2, color);

			plt::plot(p_1, p_2, {{"lw", "1"}, {"color", color}});  // add the label f(x)
  			plt::legend();
		}

		template <typename Number>
		static void plotZonotope(const Zonotope<Number> & zono, int dim1, int dim2, string color, string lw){
			Vector_t<Number> c(2);
			Matrix_t<Number> G(2, zono.generators().cols());

			c(0) = zono.center()(dim1 - 1);
			c(1) = zono.center()(dim2 - 1);

			G.row(0) = zono.generators().row(dim1 - 1);
			G.row(1) = zono.generators().row(dim2 - 1);

			Zonotope<Number> Zono = Zonotope<Number>(c, G);

			Matrix_t<Number> p = convert2polygon(Zono);
			int colNum = p.cols();
			vector<Number> p_1(colNum), p_2(colNum);
			for(auto i = 0; i < colNum; i++){
				p_1[i] = p(0,i);
				p_2[i] = p(1,i);
			}
			// Set the size of output image to 1200x780 pixels
			//plt::figure_size(1200, 780);
			// Plot line from given x and y data. Color is selected automatically.
			// plt::plot(p_1, p_2, color);

			plt::plot(p_1, p_2, {{"lw", lw}, {"color", color}});  // add the label f(x)
  			plt::legend();
		}

		template <typename Number>
		static void plotpolyZonotope(const polyZonotope<Number> & polyZono, int splits, int dim1, int dim2, string color){
			
			//Returns a polyZonotope which is projected onto the specified

			Vector_t<Number> c(2);
			Matrix_t<Number> G(2,polyZono.generators().cols());
			Matrix_t<int> expMat = polyZono.expMat();
			Matrix_t<Number> Grest(2,polyZono.Grest().cols());
			
			c(0) = polyZono.center()(dim1 - 1);
			c(1) = polyZono.center()(dim2 - 1);

			G.row(0) = polyZono.generators().row(dim1 - 1);
			G.row(1) = polyZono.generators().row(dim2 - 1);

			Grest.row(0) =polyZono.Grest().row(dim1 - 1);
			Grest.row(1) =polyZono.Grest().row(dim2 - 1);

			polyZonotope<Number> pZ = polyZonotope<Number>(c, G, expMat, Grest);
			
			//default values for the optional input arguments
			// int splits = 2;

			//delete all zero-generators
			deleteZeros(pZ);

			//split the polynomial zonotope multiple times to obtain a better over-approximation of the real shape
			vector<polyZonotope<Number>> pZsplit;
			pZsplit.push_back(pZ);

			Matrix_t<Number> polyPrev, V;
			getPolygon(pZ, V, polyPrev);
			vector<Matrix_t<Number>> Vlist;
			Vlist.push_back(V);
			
			// cout << polyPrev << endl;
			// cout << V << endl;

			Matrix_t<Number> polyAll;
			
			for(int i = 0; i < splits; i++){

				polyAll.resize(0,0);
				vector<polyZonotope<Number>> pZnew;
				vector<Matrix_t<Number>> Vnew;
				Matrix_t<Number> poly;

				for(int j = 0; j < pZsplit.size(); j++){
					// cout << "i:" << i <<", " << "j:" << j << " ";
					// int a = ismembertol(Vlist[j], polyPrev);
					// cout << a << endl;
					// cout << Vlist[j] << endl;
					// cout << polyPrev << endl;
					if(ismembertol(Vlist[j], polyPrev)){
						//split the polynomial zonotope
						polyZonotope<Number> res1, res2;
						splitLongestGen(pZsplit[j], res1, res2);

						// cout << "Good" << endl;

						//compute the corresponding polygons
						Matrix_t<Number> poly1, V1, poly2, V2;
						getPolygon(res1, V1, poly1);
						pZnew.push_back(res1);
						Vnew.push_back(V1);

						getPolygon(res2, V2, poly2);
						pZnew.push_back(res2);
						Vnew.push_back(V2);
						
						// cout << "Good" << endl;
						// cout << V1 << endl;
						// cout << poly1 << endl;

						// cout << V2 << endl;
						// cout << poly2 << endl;
						//unite the polygons
						unionPoly(poly1, poly2, poly);

						// cout << poly1 << endl;
						// cout << poly2 << endl;
						// cout << poly << endl;
												
					}else{
						Matrix_t<Number> V;
						getPolygon(pZsplit[j], V, poly);
						pZnew.push_back(pZsplit[j]);
						Vnew.push_back(V);
					}

					// cout << poly << endl;

					//unite with previous polygons
					if(polyAll.rows() == 0){
						polyAll = poly;
					}else{
						Matrix_t<Number> polyAllTemp = polyAll;
						unionPoly(polyAllTemp, poly, polyAll);
					}
				}

				// cout << polyAll << endl;

				//update lists
				pZsplit = pZnew;
				Vlist = Vnew;
				polyPrev = polyAll;
			}

			// cout << polyAll << endl;
			//plot the polygon
			int rowNum = polyAll.rows();
			vector<Number> p_1(rowNum + 1), p_2(rowNum + 1);
			for(auto i = 0; i < rowNum; i++){
				p_1[i] = polyAll(i,0);
				p_2[i] = polyAll(i,1);
			}
			p_1[rowNum] = polyAll(0,0);
			p_2[rowNum] = polyAll(0,1);
			// Set the size of output image to 1200x780 pixels
			//plt::figure_size(1200, 780);
			// Plot line from given x and y data. Color is selected automatically.
			plt::plot(p_1, p_2, color);
		}

		template <typename Number>
		static void plotReach(const polyReachableSet<Number>& R, int splits, int dim1, int dim2, string color){

        	// plt::figure_size(1200, 780);       

			polyTimeInt<Number> Rtime_interval = R.time_Rinterval();

			for(int i = 0; i < Rtime_interval.set().size(); i++){
				polyZonotope<Number> polyZ = Rtime_interval.set()[i][0].rs();
				Plotter::plotpolyZonotope(polyZ, splits, dim1, dim2, color);
			}

			// plt::show();
		}

		template <typename Number>
		static void plotReach(const polyReachableSet<Number>& R, int dim1, int dim2, string color){

        	// plt::figure_size(1200, 780);       

			polyTimeInt<Number> Rtime_interval = R.time_Rinterval();

			for(int i = 0; i < Rtime_interval.set().size(); i++){
				polyZonotope<Number> polyZ = Rtime_interval.set()[i][0].rs();
				Plotter::plotpolyZonotope(polyZ, 2, dim1, dim2, color);
			}

			// plt::show();
		}

		template <typename Number>
		static void plotReach(const ReachableSet<Number>& R, int dim1, int dim2, string color){

        	// plt::figure_size(1200, 780);       

			//polyTimeInt<Number> Rtime_interval = R.time_Rinterval();

			for(int i = 0; i < R.time_Rinterval().set().size(); i++){
				Zonotope<Number> Z = R.time_Rinterval().set()[i][0].rs();
				Plotter::plotZonotope(Z, dim1, dim2, color);
			}

			// plt::show();
		}

		template <typename Number>
		static void plotVecReach(const vector<ReachableSet<Number>>& R, int dim1, int dim2, string color){

        	// plt::figure_size(1200, 780);       

			//polyTimeInt<Number> Rtime_interval = R.time_Rinterval();

			for(int k = 0 ; k < R.size(); k++){
				for(int i = 0; i < R[k].time_Rinterval().set().size(); i++){
					Zonotope<Number> Z = R[k].time_Rinterval().set()[i][0].rs();
					Plotter::plotZonotope(Z, dim1, dim2, color);
				}
			}

			plt::show();
		}

		template <typename Number>
		static void plotOvertime(const polyReachableSet<Number>& R, int dim, string color){

        	plt::figure_size(1200, 780);       

			polyTimeInt<Number> Rtime_interval = R.time_Rinterval();

			for(int i = 0; i < Rtime_interval.set().size(); i++){
				polyZonotope<Number> Rset = Rtime_interval.set()[i][0].rs();
				const Interval_t<Number> Rtime = Rtime_interval.time(i);

				//get intervals
				IntervalMatrix<Number> intX = project(Rset, dim).interval();

				vector<Number> p_1(5), p_2(5);

				p_1[0] = Rtime.leftBound();
				p_2[0] = intX(0,0).leftBound();

				p_1[1] = Rtime.leftBound();
				p_2[1] = intX(0,0).rightBound();

				p_1[2] = Rtime.rightBound();
				p_2[2] = intX(0,0).rightBound();

				p_1[3] = Rtime.rightBound();
				p_2[3] = intX(0,0).leftBound();

				p_1[4] = p_1[0];
				p_2[4] = p_2[0];

				plt::plot(p_1, p_2, color);
			}

			plt::show();
		}

		template <typename Number>
		static void plotOvertime(const ReachableSet<Number>& R, int dim, string color){

        	plt::figure_size(1200, 780);       

			TimeInt<Number> Rtime_interval = R.time_Rinterval();

			for(int i = 0; i < Rtime_interval.set().size(); i++){
				Zonotope<Number> Rset = Rtime_interval.set()[i][0].rs();
				const Interval_t<Number> Rtime = Rtime_interval.time(i);
				
				//get intervals
				IntervalMatrix<Number> intX = project(Rset, dim).interval();

				vector<Number> p_1(5), p_2(5);

				p_1[0] = Rtime.leftBound();
				p_2[0] = intX(0,0).leftBound();

				p_1[1] = Rtime.leftBound();
				p_2[1] = intX(0,0).rightBound();

				p_1[2] = Rtime.rightBound();
				p_2[2] = intX(0,0).rightBound();

				p_1[3] = Rtime.rightBound();
				p_2[3] = intX(0,0).leftBound();

				p_1[4] = p_1[0];
				p_2[4] = p_2[0];

				plt::plot(p_1, p_2, color);
			}

			plt::show();
		}

		template <typename Number>
		static void unionPoly(const Matrix_t<Number>& poly1, 
					   const Matrix_t<Number>& poly2,
					   Matrix_t<Number>& poly)
		{
			typedef bg::model::point<double, 2, bg::cs::cartesian> point;
			typedef bg::model::polygon<point, false, true> polygon; // ccw, open polygon

			polygon p1, p2, p;
			for(int i = 0; i < poly1.rows(); i++){
				p1.outer().push_back(point(poly1(i,0),poly1(i,1)));
			}
			p1.outer().push_back(point(poly1(0,0),poly1(0,1)));
			for(int i = 0; i < poly2.rows(); i++){
				p2.outer().push_back(point(poly2(i,0),poly2(i,1)));
			}
			p2.outer().push_back(point(poly2(0,0),poly2(0,1)));
			std::vector<polygon> output;
			boost::geometry::union_(p1,p2,output);

			// cout << "hahaha good" << endl;
			// cout << output.size() << endl;

			poly.conservativeResize(output[0].outer().size() - 1,2);
			for(int i = 0; i < output[0].outer().size() - 1; i++){
				poly(i,0) = output[0].outer()[i].get<0>();
				poly(i,1) = output[0].outer()[i].get<1>();
			}

			// cout << poly << endl;
		}
	};
	//static void plotZonotope(const Zonotope<double> & zono)
} // namespace reachSolver
#ifdef SEPARATE_TEMPLATE
#include "../../src/plotter/plotter.tpp"
#endif
#endif // PLOTTER