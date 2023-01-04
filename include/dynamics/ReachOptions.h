/**
 * @file   ReachOptions.h
 * @brief  ReachOptions class representation
 * @author Yunjun Bai
 * @date April 2021
 * @version 1.0
 *
 */

/**
 * @author Changyuan Zhao
 * @date  Nov 2021
 * @version 1.1
 *
 */

#ifndef REACHOPTIONS_H
#define REACHOPTIONS_H

#include <zonotope/BasicObject.h>
#include <zonotope/Zonotope.h>
#include <zonotope/polyZonotope.h>
namespace reachSolver
{

	/**
 * @class      ReachableOptions
 * @brief      Class for ReachableOptions
 * @ingroup    structure
 * @{
 */

	template <typename Number>
	class ReachOptions
	{
	private:
		/**
   * @brief: time
   */
		double t_;

		/**
   * @brief: start time
   */
		double tStart_;

		/**
   * @brief: final time
   */
		double tFinal_;

		Zonotope<Number> R0_;

		Zonotope<Number> U_;

		Zonotope<Number> uTrans_lin_;
		Vector_t<Number> uTrans_nonlin_;
		Zonotope<Number> Rtrans_;
		Zonotope<Number> Rhom_;
		Zonotope<Number> Rhom_tp_;
		Zonotope<Number> Raux_;
		Zonotope<Number> Rpar_;
		Matrix_t<Number> uTrans_Vec_;

		/**
   * @brief: time step
   */
		double time_step_;

		/**
   * @brief: taylor terms
   */
		size_t taylor_terms_;

		/**
   * @brief: zonotope's order
   */
		size_t zonotope_order_;

		/**
   * @brief: intermediate order
   */
		size_t intermediate_order_;

		/**
   * @brief: error order
   */
		size_t error_order_;

		/**
   * @brief: tensor order
   */
		size_t tensor_order_;

		/**
   * @brief: factor
   */
		Number* factor_;

		bool verbose_;

		std::string alg_;

		/**
   * @brief: max error
   */
		Vector_t<Number> max_error_;
		int				 originContained_;

		double reductionInterval_;
		int	   usekrylovError_;

		/**
   * @brief: poly parameter
   */
		int			maxDepGenOrder_;
		double		maxPolyZonoRatio_;
		std::string restructureTechnique_;

		bool isRV_;

	public:
		/*****************************************************************************
   *                                                                           *
   *                           Constructors and Destructors                    *
   *                                                                           *
   *****************************************************************************/

		/**
   * @brief Constructor with no params
   */
		ReachOptions();

		/*****************************************************************************
   *                                                                           *
   *                       Public Functions on Properties                      *
   *                                                                           *
   *****************************************************************************/

		/**
   * @brief Get the isRV
   * @return the isRV
   */
		const double isRV() const;

		/**
   * @brief Replaces the isRV with the parameter
   * @param isRV
   */
		void set_isRV(double isRV);

		/**
   * @brief Get the t
   * @return the t
   */
		const double t() const;

		/**
   * @brief Replaces the current t with the parameter
   * @param t
   */
		void set_t(double t);

		/**
   * @brief Get the tStart
   * @return the tStart
   */
		const double tStart() const;

		/**
   * @brief Replaces the current tStart with the parameter
   * @param tStart
   */
		void set_tStart(double tStart);

		/**
   * @brief Get the tFinal
   * @return the tFinal
   */
		const double tFinal() const;

		/**
   * @brief Replaces the current tFinal with the parameter
   * @param tFinal
   */
		void set_tFinal(double tFinal);

		/**
   * @brief Get the R0
   * @return the R0
   */
		const Zonotope<Number> R0() const;

		/**
   * @brief Replaces the R0 with the parameter
   * @param R0
   */
		void set_R0(Zonotope<Number> R0);

		/**
   * @brief Get the polyR0
   * @return the polyR0
   */
		const polyZonotope<Number> polyR0() const;

		/**
   * @brief Replaces the polyR0 with the parameter
   * @param polyR0
   */
		void set_polyR0(polyZonotope<Number> polyR0);

		/**
   * @brief Get the polyR0
   * @return the polyR0
   */
		const polyZonotope<Number> polyU() const;

		/**
   * @brief Replaces the polyR0 with the parameter
   * @param polyR0
   */
		void set_polyU(polyZonotope<Number> polyU);

		/**
   * @brief Get the U
   * @return the U
   */
		const Zonotope<Number> U() const;

		/**
   * @brief Replaces the U with the parameter
   * @param U
   */
		void set_U(Zonotope<Number> U);

		/**
   * @brief Get the time_step
   * @return the time_step
   */
		const double time_step() const;

		/**
   * @brief Replaces the time_step with the parameter
   * @param time_step
   */
		void set_time_step(double time_step);

		/**
   * @brief Get the taylor_terms
   * @return the taylor_terms
   */
		const size_t taylor_terms() const;

		/**
   * @brief Replaces the taylor_terms with the parameter
   * @param taylor_terms
   */
		void set_taylor_terms(size_t taylor_terms);

		/**
   * @brief Get the zonotope_order
   * @return the zonotope_order
   */
		const size_t zonotope_order() const;

		/**
   * @brief Replaces the zonotope_order with the parameter
   * @param zonotope_order
   */
		void set_zonotope_order(size_t zonotope_order);

		/**
   * @brief Get the intermediate_order
   * @return the intermediate_order
   */
		const size_t intermediate_order() const;

		/**
   * @brief Replaces the intermediate_order with the parameter
   * @param intermediate_order
   */
		void set_intermediate_order(size_t intermediate_order);

		/**
   * @brief Get the error_order
   * @return the error_order
   */
		const size_t error_order() const;

		/**
   * @brief Replaces the error_order with the parameter
   * @param error_order
   */
		void set_error_order(size_t error_order);

		/**
   * @brief Get the tensor_order
   * @return the tensor_order
   */
		const size_t tensor_order() const;

		/**
   * @brief Replaces the tensor_order with the parameter
   * @param tensor_order
   */
		void set_tensor_order(size_t tensor_order);

		/**
   * @brief Get the factor
   * @return the factor
   */
		Number* factor() const;

		/**
   * @brief Replaces the factor with the parameter
   * @param factor
   */
		void set_factor(Number* factor);

		/**
   * @brief Get the verbose
   * @return the verbose
   */
		const bool verbose() const;

		/**
   * @brief Replaces the verbose with the parameter
   * @param verbose
   */
		void set_verbose(bool verbose);

		/**
   * @brief Get the alg
   * @return the alg
   */
		const std::string alg() const;

		/**
   * @brief Replaces the alg with the parameter
   * @param alg
   */
		void set_alg(std::string alg);

		/**
   * @brief Get the max_error
   * @return the max_error
   */
		const Vector_t<Number> max_error() const;

		/**
   * @brief Replaces the max_error with the parameter
   * @param max_error
   */
		void set_max_error(Vector_t<Number> max_error);

		/**
   * @brief Get the uTrans
   * @return the uTrans
   */
		const Zonotope<Number> uTrans_lin() const;

		/**
   * @brief Replaces the uTrans with the parameter
   * @param uTrans
   */
		void set_uTrans_lin(Zonotope<Number> uTrans_lin);

		/**
   * @brief Get the uTrans
   * @return the uTrans
   */
		const Vector_t<Number> uTrans_nonlin() const;

		/**
   * @brief Replaces the uTrans with the parameter
   * @param uTrans
   */
		void set_uTrans_nonlin(Vector_t<Number> uTrans_nonlin);

		/**
   * @brief Get the uTrans
   * @return the uTrans
   */
		const Matrix_t<Number> uTrans_Vec() const;

		/**
   * @brief Replaces the uTrans with the parameter
   * @param uTrans_Vec
   */
		void set_uTrans_Vec(Matrix_t<Number> uTrans_Vec);

		/**
   * @brief Get the Rhom
   * @return the Rhom
   */
		const Zonotope<Number> Rhom() const;

		/**
   * @brief Replaces the Rhom with the parameter
   * @param Rhom
   */
		void set_Rhom(Zonotope<Number> Rhom);

		/**
   * @brief Get the Rhom_tp
   * @return the Rhom_tp
   */
		const Zonotope<Number> Rhom_tp() const;

		/**
   * @brief Replaces the Rhom_tp with the parameter
   * @param Rhom_tp
   */
		void set_Rhom_tp(Zonotope<Number> Rhom_tp);

		/**
   * @brief Get the Raux
   * @return the Raux
   */
		const Zonotope<Number> Raux() const;

		/**
   * @brief Replaces the Raux with the parameter
   * @param Raux
   */
		void set_Raux(Zonotope<Number> Raux);

		/**
   * @brief Get the Rpar
   * @return the Rpar
   */
		const Zonotope<Number> Rpar() const;

		/**
   * @brief Replaces the Rpar with the parameter
   * @param Rpar
   */
		void set_Rpar(Zonotope<Number> Rpar);

		/**
   * @brief Get the Rtrans
   * @return the Rtrans
   */
		const Zonotope<Number> Rtrans() const;

		/**
   * @brief Replaces the Rtrans with the parameter
   * @param Rtrans
   */
		void set_Rtrans(Zonotope<Number> Rtrans);

		/**
   * @brief Get the originContained
   * @return the originContained
   */
		const int originContained() const;

		/**
   * @brief Replaces the originContained with the parameter
   * @param originContained
   */
		void set_originContained(int originContained);

		/**
   * @brief Get the reductionInterval
   * @return the reductionInterval
   */
		const double reductionInterval() const;

		/**
   * @brief Replaces the reductionInterval with the parameter
   * @param reductionInterval
   */
		void set_reductionInterval(double reductionInterval);

		/**
   * @brief Get the usekrylovError
   * @return the usekrylovError
   */
		const int usekrylovError() const;

		/**
   * @brief Replaces the usekrylovError with the parameter
   * @param usekrylovError
   */
		void set_usekrylovError(int usekrylovError);

		/**
   * @brief Get the maxDepGenOrder
   * @return the maxDepGenOrder
   */
		const int maxDepGenOrder() const;

		/**
   * @brief Replaces the maxDepGenOrder with the parameter
   * @param usekrylovError
   */
		void set_maxDepGenOrder(int maxDepGenOrder);

		/**
   * @brief Get the maxDepGenOrder
   * @return the maxDepGenOrder
   */
		const double maxPolyZonoRatio() const;

		/**
   * @brief Replaces the maxDepGenOrder with the parameter
   * @param usekrylovError
   */
		void set_maxPolyZonoRatio(double maxDepGenOrder);

		/**
   * @brief Get the maxDepGenOrder
   * @return the maxDepGenOrder
   */
		const std::string restructureTechnique() const;

		/**
   * @brief Replaces the maxDepGenOrder with the parameter
   * @param usekrylovError
   */
		void set_restructureTechnique(std::string restructureTechnique);
	};
	/** @} */
} // namespace reachSolver
#ifdef SEPARATE_TEMPLATE
#include "../../src/dynamics/ReachOptions.tpp"
#endif
#endif // REACHOPTIONS_H