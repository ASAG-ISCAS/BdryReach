/**
 * @file   NonlinearSys.h
 * @brief  NonlinearSys class representation
 * @author Yunjun Bai
 * @date April 2021
 * @version 1.0
 *
 * Reference:
 * CORA ../contDynamics/nonlinearSys/all
 */

/**
 * @author Changyuan Zhao
 * @date  Nov 2021
 * @version 1.1
 *
 */

#ifndef NONLINEARSYS_H
#define NONLINEARSYS_H

#include <capd/capdlib.h>
#include <dynamics/ContDynamics.h>
#include <dynamics/LinearSys.h>
#include <iostream>
#include <zonotope/commonType.h>
#include <zonotope/globalfunctions.h>

//#include <ctime>
// #include "derivatives.h"
namespace reachSolver
{

	/**
 * @class      NonlinearSys
 * @brief      Class for NonlinearSys.
 * @ingroup    structure
 * @{
 */
	template <typename Number>
	class NonlinearSys : private ContDynamics<Number>
	{
	private:
		// typedef Vector_t<Number> (*function_type)(Vector_t<Number> param1,
		// Vector_t<Number> param2);

		Vector_t<Number> mFile_(Vector_t<Number> vector1, Vector_t<Number> vector2);
		// autodiff::dual2nd (*mFile_f[6])(autodiff::ArrayXdual2nd& vector1,
		// autodiff::ArrayXdual2nd& vector2);
		// Number (*mFile_f_[6])(Vector_t<Number> vector1, Vector_t<Number> vector2);

		void jacobian_(IntervalMatrix<Number>  im1,
					   IntervalMatrix<Number>  im2,
					   IntervalMatrix<Number>& matrix1,
					   IntervalMatrix<Number>& matrix2);

		void jacobian_(Vector_t<Number> vector1, Vector_t<Number> vector2, Matrix_t<Number>& A, Matrix_t<Number>& B);

		std::vector<IntervalMatrix<Number>> hessian_(IntervalMatrix<Number> im1, IntervalMatrix<Number> im2);
		std::vector<Matrix_t<Number>>		hessian_(Vector_t<Number> im1, Vector_t<Number> im2);

		Matrix_t<IntervalMatrix<Number>> thirdOrderTensor(IntervalMatrix<Number>		 im1,
														  IntervalMatrix<Number>		 im2,
														  Vector_t<std::vector<size_t>>& ind);

		// FUNCTION myFun;
		// typedef void (*function_type)(capd::autodiff::Node/* t*/,
		// capd::autodiff::Node in[], int /*dimIn*/, capd::autodiff::Node out[], int/*
		// dimOut*/, capd::autodiff::Node params[], int noParams); void
		// (*_f)(capd::autodiff::Node/* t*/, capd::autodiff::Node in[], int /*dimIn*/,
		// capd::autodiff::Node out[], int/* dimOut*/, capd::autodiff::Node params[],
		// int noParams);
		capd::IMap myFun;

		// function_type thirdOrderTensor_;
		// function_type tensors_;

		linerror_type linerror_;

		/**
   * @brief: the dimesion of the output
   */
		size_t out_dimesion;

		size_t state_num;

		size_t input_num;

		bool noinput;

		// used in constructor
		// std::vector<size_t> numberOfInputs(function_type fun_handle, size_t
		// inpArgs);

		// used in reach: create an object of class reachSet that stores the reachable
		// set
		ReachableSet<Number>	 createReachSetObject(TimeInt<Number>& time_int, TimePoint<Number>& time_point);
		std::vector<ReachableSet<Number>>	 createReachSetObjectSplit(TimeInt<Number>& time_int, TimePoint<Number>& time_point);
		polyReachableSet<Number> createReachSetObject(polyTimeInt<Number>& time_int, polyTimePoint<Number>& time_point);

		// used in createReachSetObject etc. change input to an uniform
		// std::vector<ReachableSetElement<Number>> UniformOutput(TimeInt<Number>&
		// time_int); std::vector<ReachableSetElement<Number>>
		// UniformOutput(TimePoint<Number>& time_point);

		// used in initReach
		ReachableSet<Number> initReach_linRem(std::vector<ReachableSetElement<Number>>& Rinit,
											  ReachOptions<Number>&						options);

		// Splits a zonotope bundle into two zonotope bundles. This is done for one or
		// every generator resulting in n possible splits where n is the system
		// dimension; it is also possible to use a splitting hyperplane
		std::vector<Zonotope<Number>> split(Zonotope<Number> input, int number);

		// used in linReach
		LinearSys<Number> linearize(ReachOptions<Number>& options,
									Zonotope<Number>&	  R,
									ReachOptions<Number>& linOptions);
		LinearSys<Number> linearize(polyReachOptions<Number>& options,
									polyZonotope<Number>&	  R,
									polyReachOptions<Number>& linOptions);

		double linReach_linRem(LinearReachableSet<Number>& R,
							   Zonotope<Number>&		   Rinit,
							   Zonotope<Number>&		   Rdelta,
							   ReachOptions<Number>&	   options,
							   Zonotope<Number>&		   Rti,
							   Zonotope<Number>&		   Rtp);

	public:
		/****************************************************************************
   *                                                                           *
   *                           Constructors and Destructors                    *
   *                                                                           *
   *****************************************************************************/

		/**
   * @brief Constructor with no params
   */
		// NonlinearSys();

		/**
   * @brief Constructor with function handle to the dynamic equation
   * @param fun_handle function handle to the dynamic equation
   */
		// NonlinearSys(size_t dimension);

		/**
   * @brief Constructor with two params
   * @param name name of dynamics
   * @param fun_handle function handle to the dynamic equation
   */
		// NonlinearSys(std::string name, function_type fun_handle);

		/**
   * @brief Constructor with three params
   * @param fun_handle function handle to the dynamic equation
   * @param num_states number of states
   * @param num_inputs number of inputs
   * @param out_dimension dimension of output
   */
		NonlinearSys(capd::IMap fun_handle, size_t num_states, size_t num_inputs, size_t dimension);

		/**
   * @brief Constructor with four params
   * @param name name of dynamics
   * @param fun_handle function handle to the dynamic equation
   * @param num_states number of states
   * @param num_inputs number of inputs
   * @param out_dimension dimension of output
   */
		NonlinearSys(std::string name, capd::IMap fun_handle, size_t num_states, size_t num_inputs, size_t dimension);

		/**
   * @brief Copy Constructor - constructs a NonlinearSys from an existing one.
   * @param other Another NonlinearSys, from which a new NonlinearSys is
   * constructed
   */
		NonlinearSys(const NonlinearSys& other) = default;

		virtual ~NonlinearSys();

		/*****************************************************************************
   *                                                                           *
   *                       Public Functions on Properties                      *
   *                                                                           *
   *****************************************************************************/

		/**
   * @brief Get the mFile handle
   * @return mFile handle
   */
		// const function_type mFile() const;

		/**
   * @brief Get the jacobian
   * @return jacobian handle
   */
		// const function_type jacobian() const;

		/**
   * @brief Get the hessian
   * @return hessian handle
   */
		// const function_type hessian() const;

		/**
   * @brief Get the thirdOrderTensor
   * @return thirdOrderTensor handle
   */
		// const function_type thirdOrderTensor() const;

		/**
   * @brief Get the tensors
   * @return tensors handle
   */
		// const function_type tensors() const;

		/**
   * @brief Get the linerror
   * @return linerror
   */
		// const linerror_type linerror() const;

		/**
   * @brief Set the linerror
   * @param linerror
   */
		// void set_linerror(linerror_type linerror) const;

		/*****************************************************************************
   *                                                                           *
   *                       Compute Reachable Set                               *
   *                                                                           *
   *****************************************************************************/

		/**
   * @brief computes the reachable continuous set for the entire time horizon of
   * a continuous system
   * @param options options for the computation of reachable sets
   * @param spec object of class specification
   * @param R object of class reachSet storing the computed reachable set
   * @return 1 if specifications are satisfied, 0 if not
   */
		// int reach(ReachOptions<Number>& options, ReachSpecification& spec,
		// ReachableSet<Number> & R);
		int reach(ReachOptions<Number>& options, ReachableSet<Number>& R);
		int reachSplit(ReachOptions<Number>& options, vector<ReachableSet<Number>>& R);
		int reach(polyReachOptions<Number>& options, polyReachableSet<Number>& R);

		/**
   * @brief checks if all necessary options are there and have valid values
   * @param options options for nonlinearSys
   * @param hyb called from hybrid Automaton (0/1)
   * @return options for nonlinearSys
   */
		void checkOptionsReach(ReachOptions<Number>& options, int hyb);
		void checkOptionsReach(polyReachOptions<Number>& options, int hyb);

		/**
   * @brief computes multivariate derivatives (jacobians, hessians, etc.) of
   * nonlinear systems in a symbolic way; the result is stored in m-files and
   * passed by a handle
   * @param options options struct
   */
		void derivatives(ReachOptions<Number>& options);

		/**
   * @brief computes the reachable continuous set for the first time step
   * @param Rinit initial reachable set
   * @param options struct containing the algorithm settings
   * @return first reachable set
   */
		ReachableSet<Number> initReach(std::vector<ReachableSetElement<Number>> Rinit, ReachOptions<Number>& options);
		polyReachableSet<Number> initReach(std::vector<polyReachableSetElement<Number>> Rinit,
										   polyReachOptions<Number>&					options);

		/**
   * @brief computes the reachable continuous set for one time step of a
   * nonlinear system by overapproximative linearization
   * @param R reachable set of the previous time step
   * @param options options for the computation of the reachable set
   * @return reachable set of the next time step
   */
		ReachableSet<Number>	 post(ReachableSet<Number>& R, ReachOptions<Number>& options);
		polyReachableSet<Number> post(polyReachableSet<Number>& R, polyReachOptions<Number>& options);

		/**
   * @brief computes the reachable set after linearization
   * @param options struct with algorithm settings
   * @param Rstart initial reachable set
   * @param Rti reachable set for time interval
   * @param Rtp reachable set for time point
   * @return dimForSplit - dimension that is split to reduce the lin. error
   */
		int linReach(ReachOptions<Number>&		  options,
					 ReachableSetElement<Number>& Rstart,
					 ReachableSetElement<Number>& Rti,
					 ReachableSetElement<Number>& Rtp);
		int linReach(polyReachOptions<Number>&		  options,
					 polyReachableSetElement<Number>& Rstart,
					 polyReachableSetElement<Number>& Rti,
					 polyReachableSetElement<Number>& Rtp);

		/**
   * @brief computes the solution due to the linearization error
   * @param options - options struct
   * @param R - reachable set (time-interval solution from linearized system +
   * estimated set of abstraction errors)
   * @param Verrordyn - abstraction error (zonotope)
   * @return trueError - abstraction error (interval)
   */
		Vector_t<Number> abstrerr_lin(ReachOptions<Number>& options, Zonotope<Number> R, Zonotope<Number>& VerrorDyn);

		/**
   * @brief: standardized console output if options verbose = true
   * @param step current step
   * @param t start time of current step
   * @param options options struct
   */
		void verboseLog(int step, double t, ReachOptions<Number>& options);
		void verboseLog(int step, double t, polyReachOptions<Number>& options);

		/**
   * @brief: precompute thr second order static error along with hessian matrix
   * @param Rdelta shifted reachable set at the beginning of the time step
   * @param H hessian matrix
   * @param Zdelta zonotope over-approximating the reachable set at the
   * beginning of the time step extended by the input set
   * @param errorStat static linearization error
   * @param T third-order tensor
   * @param ind3 indices at which the third-order tensor is not zero
   * @param Zdelta3 set Zdelta reduced to the zonotope order for the evaluation
   * of the third-order tensor
   */
		void precompStatError(polyZonotope<Number>				Rdelta,
							  polyReachOptions<Number>&			options,
							  std::vector<Matrix_t<Number>>&	H,
							  Zonotope<Number>&					Zdelta,
							  polyZonotope<Number>&				errorStat,
							  Matrix_t<IntervalMatrix<Number>>& T,
							  std::vector<size_t>&				ind3,
							  Zonotope<Number>&					Zdelta3);

		/**
   * @brief: computes the abstacion error for the polymiolization approach
   * @param options options struct
   * @param Rall time-interval reachable set
   * @param Rdiff difference between the reachable set at the beginning of the
   * time interval and the time-interval reachable set
   * @param H Hessian matrix
   * @param Zdelta zonotope over-approximating the reachable set at the
   * beginning o the time step extended by the input errors
   * @param VerrorStat set of static linearization errors
   * @param T third-order tensor
   * @param ind3 indices of non-zero entries in the third-order tensor
   * @param Zdelta3 set Zdelta reduced to the zonotope order for the evaluation
   * of the third-order tensor
   * @param VerrorDyn zonotope overapproximating the dynamic linearizatiion
   * error
   * @param VerrorStat2 zonoyope overapproximationg the static linearization
   * error
   * @return trueError interval overapproxiamting the overall linearization
   * error
   */
		Vector_t<Number> abstrerr_poly(polyReachOptions<Number>&		options,
									   polyZonotope<Number>				Rall,
									   polyZonotope<Number>				Rdiff,
									   std::vector<Matrix_t<Number>>	H,
									   Zonotope<Number>					Zdelta,
									   polyZonotope<Number>				VerrorStat,
									   Matrix_t<IntervalMatrix<Number>> T,
									   std::vector<size_t>				ind3,
									   Zonotope<Number>					Zdelta3,
									   Zonotope<Number>&				VerrorDyn,
									   polyZonotope<Number>&			VerrorStat2);

		/**
   * @brief: selects the split strategy of the reachable set causing the least
   * linearization error
   * @param  options struct containin the algortithm setting
   * @param Rinit initial reachables set
   * @return dimForSplit dimension that is split to reduce the linearization
   * error
   */
		int select(ReachOptions<Number> options, ReachableSetElement<Number>& Rinit);
	};

	    /**
   * @brief: Calculate the approximate ratio of the volumes 
   * between the dependent generator and the independent 
   * generator part of the polynomial zonotope 
   * ratio = (V_ind/V_dep)^(1/n)
   * @param polyZonotope polyZono
   * @return ratio
   */
    template <typename Number>
	Number approxVolumeRatio(const polyZonotope<Number>& pZ, const polyReachOptions<Number>& options);

	/** @} */
} // namespace reachSolver

#ifdef SEPARATE_TEMPLATE
#include "../../src/dynamics/NonlinearSys.tpp"
#endif
#endif // NONLINEARSYS_H