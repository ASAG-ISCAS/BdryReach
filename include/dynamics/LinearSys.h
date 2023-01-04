/**
 * @file   LinearSys.h
 * @brief  LinearSys class representation
 * @author Yunjun Bai
 * @date April 2021
 * @version 1.0
 *
 * Reference:
 * CORA ../contDynamics/linearSys/all
 */

/**
 * @author Changyuan Zhao
 * @date  Nov 2021
 * @version 1.1
 *
 */

#ifndef LINEARSYS_H
#define LINEARSYS_H

#include <dynamics/ContDynamics.h>

namespace reachSolver
{

	/**
 * @class      LinearSys
 * @brief      Class for LinearSys.
 * @ingroup    structure
 * @{
 */
	template <typename Number>
	class LinearSys : private ContDynamics<Number>
	{
	private:
		/**
   * @brief: represent of continuous dynamics by taylor type
   */
		struct taylor_type
		{
			/**
     * @brief: time step
     */
			double timeStep;
			/**
     * @brief: eAt
     */
			Matrix_t<Number> eAt;

			/**
     * @brief: F
     */
			IntervalMatrix<Number> F;

			/**
     * @brief: RV
     */
			Zonotope<Number> RV;

			/**
     * @brief: inputCorr
     */
			Zonotope<Number> inputCorr;

			/**
     * @brief: Rtrans
     */
			Zonotope<Number> Rtrans;

			/**
     * @brief: powers
     */
			std::vector<Matrix_t<Number>> powers;

			/**
     * @brief: error interval matrix
     */
			IntervalMatrix<Number> error;

			/**
     * @brief: inputF interval matrix
     */
			IntervalMatrix<Number> inputF;

			/**
     * @brief: zonotope V
     */
			Zonotope<Number> V;

			/**
     * @brief: eAtInt
     */
			IntervalMatrix<Number> eAtInt;
		};

		Matrix_t<Number> A_; // system matrix
		Matrix_t<Number> B_; // input matrix
		Matrix_t<Number> c_; // constant input
		Matrix_t<Number> C_; // output matrix
		Matrix_t<Number> D_; // throughput matrix
		Matrix_t<Number> k_; // output offset
		taylor_type		 taylor_;
		Matrix_t<Number> krylov_;

		/**
   * @brief computes the overapproximation of the exponential of a system matrix
   * up to a certain accuracy
   * @param options struct containing the algorithm settings
   */
		void exponential(ReachOptions<Number>& options);
		void exponential(polyReachOptions<Number>& options);

		/**
   * @brief tie=time interval error; computes the error done by building the
   * convex hull of time point solutions
   * @param options struct containing the algorithm settings
   */
		void tie(ReachOptions<Number>& options);
		void tie(polyReachOptions<Number>& options);

		/**
   * @brief computes the bloating due to the input
   * @param options struct containing the algorithm settings
   */
		void inputSolution(ReachOptions<Number>& options);
		void inputSolution(polyReachOptions<Number>& options);

		/**
   * @brief computes the error done by the linear assumption of the constant
   * input solution
   * @param options struct containing the algorithm settings
   */
		void inputTie(ReachOptions<Number>& options);
		void inputTie(polyReachOptions<Number>& options);

	public:
		/****************************************************************************
   *                                                                           *
   *                           Constructors and Destructors                    *
   *                                                                           *
   *****************************************************************************/

		/**
   * @brief Constructor with no params
   */
		// LinearSys();

		/**
   * @brief Constructor with two params
   * @param A system matrix
   * @param B input matrix
   */
		LinearSys(Matrix_t<Number>& A, Matrix_t<Number>& B);

		/**
   * @brief Constructor with three params
   * @param A system matrix
   * @param B input matrix
   * @param c constant input
   */
		LinearSys(Matrix_t<Number>& A, Matrix_t<Number>& B, Matrix_t<Number>& c);

		/**
   * @brief Constructor with four params
   * @param A system matrix
   * @param B input matrix
   * @param c constant input
   * @param C output matrix
   */
		LinearSys(Matrix_t<Number>& A, Matrix_t<Number>& B, Matrix_t<Number>& c, Matrix_t<Number>& C);

		/**
   * @brief Constructor with five params
   * @param A system matrix
   * @param B input matrix
   * @param c constant input
   * @param C output matrix
   * @param D throughput matrix
   */
		LinearSys(
			Matrix_t<Number>& A, Matrix_t<Number>& B, Matrix_t<Number>& c, Matrix_t<Number>& C, Matrix_t<Number>& D);

		/**
   * @brief Constructor with six params
   * @param A system matrix
   * @param B input matrix
   * @param c constant input
   * @param C output matrix
   * @param D throughput matrix
   * @param k output offset
   */
		LinearSys(Matrix_t<Number>& A,
				  Matrix_t<Number>& B,
				  Matrix_t<Number>& c,
				  Matrix_t<Number>& C,
				  Matrix_t<Number>& D,
				  Matrix_t<Number>& k);

		/**
   * @brief Constructor with name and two params
   * @param name name of LinearSys
   * @param A system matrix
   * @param B input matrix
   */
		LinearSys(std::string name, Matrix_t<Number>& A, Matrix_t<Number>& B);

		/**
   * @brief Constructor with name and three params
   * @param name name of LinearSys
   * @param A system matrix
   * @param B input matrix
   * @param c constant input
   */
		LinearSys(std::string name, Matrix_t<Number>& A, Matrix_t<Number>& B, Matrix_t<Number>& c);

		/**
   * @brief Constructor with name and four params
   * @param name name of LinearSys
   * @param A system matrix
   * @param B input matrix
   * @param c constant input
   * @param C output matrix
   */
		LinearSys(std::string name, Matrix_t<Number>& A, Matrix_t<Number>& B, Matrix_t<Number>& c, Matrix_t<Number>& C);

		/**
   * @brief Constructor with name and five params
   * @param name name of LinearSys
   * @param A system matrix
   * @param B input matrix
   * @param c constant input
   * @param C output matrix
   * @param D throughput matrix
   */
		LinearSys(std::string		name,
				  Matrix_t<Number>& A,
				  Matrix_t<Number>& B,
				  Matrix_t<Number>& c,
				  Matrix_t<Number>& C,
				  Matrix_t<Number>& D);

		/**
   * @brief Constructor with name and six params
   * @param name name of LinearSys
   * @param A system matrix
   * @param B input matrix
   * @param c constant input
   * @param C output matrix
   * @param D throughput matrix
   * @param k output offset
   */
		LinearSys(std::string		name,
				  Matrix_t<Number>& A,
				  Matrix_t<Number>& B,
				  Matrix_t<Number>& c,
				  Matrix_t<Number>& C,
				  Matrix_t<Number>& D,
				  Matrix_t<Number>& k);

		/**
   * @brief Copy Constructor - constructs a NonlinearSys from an existing one.
   * @param other Another NonlinearSys, from which a new NonlinearSys is
   * constructed
   */
		LinearSys(const LinearSys& other) = default;

		virtual ~LinearSys();

		/*****************************************************************************
   *                                                                           *
   *                       Public Functions on Properties                      *
   *                                                                           *
   *****************************************************************************/

		/*****************************************************************************
   *                                                                           *
   *                       Compute Reachable Set                               *
   *                                                                           *
   *****************************************************************************/

		/**
   * @brief computes the reachable continuous set for the first time step
   * @param Rinit initial reachable set
   * @param options struct containing the algorithm settings
   * @return first reachable set
   */
		LinearReachableSet<Number>	   initReach(Zonotope<Number>& Rinit, ReachOptions<Number>& options);
		polyLinearReachableSet<Number> initReach(polyZonotope<Number>& Rinit, polyReachOptions<Number>& options);

		/**
   * @brief computes the reachable continuous set for the first time step using
   * Krylov subspace methods
   * @param Rinit initial reachable set
   * @param options struct containing the algorithm settings
   * @return first reachable set
   */
		LinearReachableSet<Number>	   initReach_Krylov(Zonotope<Number>& Rinit, ReachOptions<Number>& options);
		polyLinearReachableSet<Number> initReach_Krylov(polyZonotope<Number>& Rinit, polyReachOptions<Number>& options);

		/**
   * @brief computes the reachable continuous set for the first time step in the
   * untransformed space
   * @param Rinit initial reachable set
   * @param options struct containing the algorithm settings
   * @return first reachable set
   */
		LinearReachableSet<Number>	   initReach_Euclidean(Zonotope<Number>& Rinit, ReachOptions<Number>& options);
		polyLinearReachableSet<Number> initReach_Euclidean(polyZonotope<Number>&	 Rinit,
														   polyReachOptions<Number>& options);

		/**
   * @brief computes the solution due to the linearization error
   * @param options struct containing the algorithm settings
   * @param Vdyn set of admissible errors (dynamic)
   * @param Vstat - set of admissible errors (static) (optional)
   * @return  reachable set due to the linearization error
   */
		Zonotope<Number> error_solution(ReachOptions<Number>& options, Zonotope<Number> Vdyn, Zonotope<Number> Vstat);
		Zonotope<Number> error_solution(polyReachOptions<Number>& options,
										Zonotope<Number>		  Vdyn,
										Zonotope<Number>		  Vstat);
		polyZonotope<Number> error_solution(polyReachOptions<Number>& options,
											Zonotope<Number>		  Vdyn,
											polyZonotope<Number>	  Vstat);

		/**
   * @brief: computes the reachable continuous set of different to the initial
   * sstate for all initial states
   * @param options - options struct
   * @param Rinit - initial reachable set
   * @return Rdelta difference reachanle set
   */
		polyZonotope<Number> deltaReach(polyReachOptions<Number>& options, polyZonotope<Number> Rinit);

		void lin_print();
	};

	/** @} */
} // namespace reachSolver
#ifdef SEPARATE_TEMPLATE
#include "../../src/dynamics/LinearSys.tpp"
#endif
#endif // LINEARSYS_H