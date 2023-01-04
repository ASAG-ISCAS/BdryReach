/*
 * @Author: Changyuan zhao
 * @Date: 2021-11-18 22:22:36
 */

#ifndef POLYREACHABLESETELEMENT_H
#define POLYREACHABLESETELEMENT_H

#include <dynamics/ReachOptions.h>
#include <zonotope/BasicObject.h>
#include <zonotope/Zonotope.h>
#include <zonotope/polyZonotope.h>
namespace reachSolver
{

	/**
 * @class      polyReachableSetElement
 * @brief      Class for polyReachableSetElement.
 * @ingroup    structure
 * @{
 */
	template <typename Number>
	class polyReachableSetElement
	{
	private:
		/**
   * @brief: a set represented by a polyZonotope
   */
		polyZonotope<Number> set_;

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
		polyReachableSetElement();

		/**
   * @brief Constructor with two params
   * @param set polyZonotope<Number>
   * @param error error vector
   */
		polyReachableSetElement(polyZonotope<Number> set, Vector_t<Number> error);

		/*****************************************************************************
   *                                                                           *
   *                       Public Functions on Properties                      *
   *                                                                           *
   *****************************************************************************/
		/**
   * @brief Get the set
   * @return the set
   */
		polyZonotope<Number> rs() const;

		/**
   * @brief Replaces the current set with the parameter
   * @param set
   */
		void set_rs(polyZonotope<Number>& set);

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
#include "../../src/dynamics/polyReachableSetElement.tpp"
#endif
#endif // POLYREACHABLESETELEMENT_H