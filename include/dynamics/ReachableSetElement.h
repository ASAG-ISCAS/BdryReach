/*
 * @Author: Changyuan zhao
 * @Date: 2021-11-18 22:22:36
 */

#ifndef REACHABLESETELEMENT_H
#define REACHABLESETELEMENT_H

#include <dynamics/ReachOptions.h>
#include <zonotope/BasicObject.h>
#include <zonotope/Zonotope.h>
namespace reachSolver
{

	/**
 * @class      ReachableSetElement
 * @brief      Class for ReachableSetElement.
 * @ingroup    structure
 * @{
 */
	template <typename Number>
	class ReachableSetElement
	{
	private:
		/**
   * @brief: a set represented by a zonotope
   */
		Zonotope<Number> set_;

		/**
   * @brief: error vector
   */
		Vector_t<Number> error_;

		int prev_;
		int parent_;

	public:
		/****************************************************************************
   *                                                                           *
   *                           Constructors and Destructors                    *
   *                                                                           *
   *****************************************************************************/
		/**
   * @brief Constructor with no params
   */
		ReachableSetElement();

		/**
   * @brief Constructor with two params
   * @param set Zonotope<Number>
   * @param error error vector
   */
		ReachableSetElement(Zonotope<Number> set, Vector_t<Number> error);

		/*****************************************************************************
   *                                                                           *
   *                       Public Functions on Properties                      *
   *                                                                           *
   *****************************************************************************/
		/**
   * @brief Get the set
   * @return the set
   */
		Zonotope<Number> rs() const;

		/**
   * @brief Replaces the current set with the parameter
   * @param set
   */
		void set_rs(Zonotope<Number>& set);

		/**
   * @brief Get the error
   * @return the error
   */
		Vector_t<Number> error() const;

		/**
   * @brief Replaces the current error with the parameter
   * @param time
   */
		void set_error(Vector_t<Number> error);

		/**
   * @brief Get the prev
   * @return the prev
   */
		const int prev() const;

		/**
   * @brief Replaces the current prev with the parameter
   * @param prev
   */
		void set_prev(int prev);

		/**
   * @brief Get the parent
   * @return the parent
   */
		const int parent() const;

		/**
   * @brief Replaces the current parent with the parameter
   * @param parent
   */
		void set_parent(int parent);

		/**
   * @brief: get set
   * @return set
   */
	};
} // namespace reachSolver
#ifdef SEPARATE_TEMPLATE
#include "../../src/dynamics/ReachableSetElement.tpp"
#endif
#endif // REACHABLESETELEMENT_H