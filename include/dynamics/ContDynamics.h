/**
 * @file   ContDynamics.h
 * @brief  ContDynamics class representation
 * @author Yunjun Bai
 * @date April 2021
 * @version 1.0
 *
 * Reference:
 * CORA ../contDynamics/ContDynamics/all
 */

/**
 * @author Changyuan Zhao
 * @date  Nov 2021
 * @version 1.1
 *
 */
#ifndef CONTDYNAMICS_H
#define CONTDYNAMICS_H

#include <algorithm>
#include <cfloat>
#include <cmath>
#include <dynamics/ReachOptions.h>
#include <dynamics/ReachableSet.h>
#include <dynamics/TimeRelated.h>
#include <dynamics/polyReachOptions.h>
#include <exception>
#include <string>
#include <vector>
#include <zonotope/BasicObject.h>
#include <zonotope/globalfunctions.h>

namespace reachSolver
{

	struct SetExplosionException : public std::exception
	{
		const char* what() const throw() { return "Set Explosion Exception!"; }
	};

	/**
 * @class      ContDynamics
 * @brief      Class for ContDynamics.
 * @ingroup    structure
 * @{
 */
	template <typename Number>
	class ContDynamics
	{
	private:
		// typedef Number (*function_type)(Number param1, Number param2);

		/**
   * @brief: name
   */
		std::string name_;

		/**
   * @brief: dimension
   */
		size_t dim_;

		/**
   * @brief: input number
   */
		size_t num_inputs_;

		/**
   * @brief: output number
   */
		size_t num_outputs_;

		// used in constructor
		// std::vector<size_t> numberOfInputs(function_type fun_handle, size_t
		// inpArgs);

		// used in initReach
		// ReachableSet<Number>
		// initReach_linRem(std::vector<ReachableSetElement<Number>>& Rinit,
		// ReachOptions<Number>& options);

		// std::vector<Zonotope<Number>> split(Zonotope<Number> input, int number);

	public:
		/****************************************************************************
   *                                                                           *
   *                           Constructors and Destructors                    *
   *                                                                           *
   *****************************************************************************/

		/**
   * @brief Constructor with no params
   */
		// ContDynamics();

		/**
   * @brief Constructor with function handle to the dynamic equation
   * @param name name of dynamics
   */
		explicit ContDynamics(std::string name);

		/**
   * @brief Constructor with four params
   * @param name name of dynamics
   * @param dim number of states
   * @param num_inputs number of inputs
   * @param num_outputs number of outputs
   */
		ContDynamics(std::string name, size_t dim, size_t num_inputs, size_t num_outputs);

		virtual ~ContDynamics();

		/*****************************************************************************
   *                                                                           *
   *                       Public Functions on Properties                      *
   *                                                                           *
   *****************************************************************************/
		const std::string name() const;

		const size_t num_outputs() const;

		const size_t dim() const;

		const size_t num_inputs() const;

		void set_num_outputs(const size_t num_outputs);

		void set_dim(const size_t dim);

		void set_num_inputs(const size_t num_inputs);

		/*****************************************************************************
   *                                                                           *
   *                       Compute Reachable Set                               *
   *                                                                           *
   *****************************************************************************/
	};

	/** @} */
} // namespace reachSolver
#ifdef SEPARATE_TEMPLATE
#include "../../src/dynamics/ContDynamics.tpp"
#endif
#endif // CONTDYNAMICS_H