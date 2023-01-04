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
#include <zonotope/polyZonotope.h>
#include <dynamics/LinearSys.h>
#include <dynamics/NonlinearSys.h>
#include <eigen3/Eigen/Core>
#include <functional>
#include <iostream>
#include <cmath>

#include <plotter/matplotlibcpp.h>
#include <plotter/plotter.h>

#include<underApprox/underApprox.h>

using namespace std;
using namespace capd;
using capd::autodiff::Node;


namespace reachSolver
{
	class OverApprox
	{
	public:
		/**
		 * @brief get the over-approximation using the boundary theory with zonotope as set representation.
		 */
		template <typename Number>
		static vector<ReachableSet<Number>> BdReach(NonlinearSys<Number> mysys, ReachOptions<Number> options, Zonotope<Number> R0){

			vector<Zonotope<Number>> R0_border = Zono_border(R0);

			vector<ReachableSet<Number>> BdReachset(R0_border.size());
				for(int i = 0; i < R0_border.size(); i++){

					options.set_R0(R0_border[i]);
					mysys.reach(options, BdReachset[i]);
				}

			return BdReachset;
		}


		/**
		 * @brief get the over-approximation using the boundary theory with zonotope as set representation(Split the boundary of the R0)
		 */
		template <typename Number>
		static vector<ReachableSet<Number>> BdReachSplit(NonlinearSys<Number> mysys, ReachOptions<Number> options, Zonotope<Number> R0, double radius){

			vector<Zonotope<Number>> R0_border = Zono_border(R0);

			vector<Zonotope<Number>> R0_border_split = UnderApprox::boundSplit(R0_border, radius);

			cout << R0_border_split.size() << endl;

			vector<ReachableSet<Number>> BdReachset(R0_border_split.size());
				for(int i = 0; i < R0_border_split.size(); i++){

					options.set_R0(R0_border_split[i]);
					mysys.reach(options, BdReachset[i]);
				}

			return BdReachset;
		}
				/**
		 * @brief get the over-approximation using the boundary theory with zonotope as set representation.
		 */
		template <typename Number>
		static vector<polyReachableSet<Number>> BdReach(NonlinearSys<Number> mysys, polyReachOptions<Number> options, Zonotope<Number> R0){

			vector<Zonotope<Number>> R0_border = Zono_border(R0);

			vector<polyReachableSet<Number>> BdReachset(R0_border.size());
				for(int i = 0; i < R0_border.size(); i++){

					options.set_R0(polyZonotope<Number>(R0_border[i]));
					mysys.reach(options, BdReachset[i]);
				}

			return BdReachset;
		}

	
	};
} // namespace reachSolver
