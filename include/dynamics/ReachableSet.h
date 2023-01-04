/**
 * @file   ReachableSet.h
 * @brief  ReachableSet class representation
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

#ifndef REACHABLESET_H
#define REACHABLESET_H

#include <dynamics/ReachOptions.h>
#include <dynamics/ReachableSetElement.h>
#include <dynamics/TimeRelated.h>
#include <dynamics/polyReachOptions.h>
#include <dynamics/polyReachableSetElement.h>
#include <zonotope/BasicObject.h>
#include <zonotope/Zonotope.h>
namespace reachSolver
{

	/** @} */

	/**
 * @class      ReachableSet
 * @brief      Class for ReachableSet.
 * @ingroup    structure
 * @{
 */
	template <typename Number>
	class ReachableSet
	{
	private:
		/**
   * @brief: time point
   */
		std::vector<ReachableSetElement<Number>> time_point_;

		/**
   * @brief: time interval
   */
		std::vector<ReachableSetElement<Number>> time_interval_;

		/**
   * @brief: result time interval
   */
		TimeInt<Number> Rtime_interval_;

		/**
   * @brief: result time point
   */
		TimePoint<Number> Rtime_point_;

		/**
   * @brief: reachable set
   */
		std::vector<ReachableSetElement<Number>> R0_;
		int										 parent_rs_; // index of parent reachable set
		int										 loc_;		 // index of the location

		int internalCount_;

	public:
		/****************************************************************************
   *                                                                           *
   *                           Constructors and Destructors                    *
   *                                                                           *
   *****************************************************************************/
		/**
   * @brief Constructor with no params
   */
		ReachableSet();

		/**
   * @brief Constructor with TimePoint
   * @param time_point struct with fields .set and .time storing the time point
   * reachable set
   */
		explicit ReachableSet(std::vector<ReachableSetElement<Number>>& time_point);

		/**
   * @brief Constructor with two params
   * @param time_point struct with fields .set and .time storing the time point
   * reachable set
   * @param parent index of the parent reachable set
   */
		ReachableSet(std::vector<ReachableSetElement<Number>>& time_point, int parent);

		/**
   * @brief Constructor with two params
   * @param time_point struct with fields .set and .time storing the time point
   * reachable set
   * @param time_interval struct with fields .set, .time, and .algebraic
   * (nonlinDAsys) storing the time interval reachable set
   */
		ReachableSet(std::vector<ReachableSetElement<Number>>& time_point,
					 std::vector<ReachableSetElement<Number>>& time_interval);

		/**
   * @brief Constructor with three params
   * @param time_point struct with fields .set and .time storing the time point
   * reachable set
   * @param parent index of the parent reachable set
   * @param loc index of the location
   */
		ReachableSet(std::vector<ReachableSetElement<Number>>& time_point, int parent, int loc);

		/**
   * @brief Constructor with three params
   * @param time_point struct with fields .set and .time storing the time point
   * reachable set
   * @param time_interval struct with fields .set, .time, and .algebraic
   * (nonlinDAsys) storing the time interval reachable set
   * @param parent index of the parent reachable set
   */
		ReachableSet(std::vector<ReachableSetElement<Number>>& time_point,
					 std::vector<ReachableSetElement<Number>>& time_interval,
					 int									   parent);

		/**
   * @brief Constructor with four params
   * @param time_point struct with fields .set and .time storing the time point
   * reachable set
   * @param time_interval struct with fields .set, .time, and .algebraic
   * (nonlinDAsys) storing the time interval reachable set
   * @param parent index of the parent reachable set
   * @param loc index of the location
   */
		ReachableSet(std::vector<ReachableSetElement<Number>>& time_point,
					 std::vector<ReachableSetElement<Number>>& time_interval,
					 int									   parent,
					 int									   loc);

		/**
   * @brief Copy Constructor - constructs a ReachableSet from an existing one.
   * @param other Another ReachableSet, from which a new ReachableSet is
   * constructed
   */
		ReachableSet(const ReachableSet<Number>& other) = default;

		virtual ~ReachableSet();

		/*****************************************************************************
   *                                                                           *
   *                       Public Functions on Properties                      *
   *                                                                           *
   *****************************************************************************/

		/**
   * @brief Get the time_point
   * @return time_point
   */
		const std::vector<ReachableSetElement<Number>> time_point() const;

		/**
   * @brief Replaces the current time_point with the parameter
   * @param time_point
   */
		void set_time_point(std::vector<ReachableSetElement<Number>>& time_point);

		/**
   * @brief Get the time_interval
   * @return time_interval
   */
		const std::vector<ReachableSetElement<Number>> time_interval() const;

		/**
   * @brief Replaces the current time_interval with the parameter
   * @param time_interval
   */
		void set_time_interval(std::vector<ReachableSetElement<Number>>& time_interval);

		/**
   * @brief Get the result time_point
   * @return time_point
   */
		TimePoint<Number> Rtime_point() const;

		/**
   * @brief Replaces the result time_point with the parameter
   * @param time_point
   */
		void set_Rtime_point(TimePoint<Number>& time_point);

		/**
   * @brief Get the result time_interval
   * @return time_interval
   */
		TimeInt<Number> time_Rinterval() const;

		/**
   * @brief Replaces the result time_interval with the parameter
   * @param time_interval
   */
		void set_Rtime_interval(TimeInt<Number>& time_interval);

		/**
   * @brief Get the time_interval
   * @return time_interval
   */
		const std::vector<ReachableSetElement<Number>> R0() const;

		/**
   * @brief Replaces the current R0 with the parameter
   * @param R0
   */
		void set_R0(std::vector<ReachableSetElement<Number>>& R0);

		/**
   * @brief Get the parent_rs
   * @return parent_rs
   */
		const int parent_rs() const;

		/**
   * @brief Replaces the current parent_rs with the parameter
   * @param parent_rs
   */
		void set_parent_rs(int parent_rs);

		/**
   * @brief Get the loc
   * @return loc
   */
		const int loc() const;

		/**
   * @brief Replaces the current loc with the parameter
   * @param loc
   */
		void set_loc(int loc);

		/**
   * @brief Get the internalCount
   * @return internalCount
   */
		const int internalCount() const;

		/**
   * @brief Replaces the current internalCount with the parameter
   * @param internalCount
   */
		void set_internalCount(int internalCount);

		/*****************************************************************************
   *                                                                           *
   *                       Other Functions                                     *
   *                                                                           *
   *****************************************************************************/

		/**
   * @brief delete reachable sets that are already covered by other sets
   * @param Rold reachable sets of previous time steps
   * @param options options for the computation of the reachable set
   */
		void deleteRedundantSets(ReachableSet<Number> Rold, ReachOptions<Number>& options);
	};

	/**
 * @class      polyReachableSet
 * @brief      Class for polyReachableSet.
 * @ingroup    structure
 * @{
 */
	template <typename Number>
	class polyReachableSet
	{
	private:
		/**
   * @brief: time point
   */
		std::vector<polyReachableSetElement<Number>> time_point_;

		/**
   * @brief: time interval
   */
		std::vector<polyReachableSetElement<Number>> time_interval_;

		/**
   * @brief: result time interval
   */
		polyTimeInt<Number> Rtime_interval_;

		/**
   * @brief: result time point
   */
		polyTimePoint<Number> Rtime_point_;

		/**
   * @brief: reachable set
   */
		std::vector<polyReachableSetElement<Number>> R0_;
		int											 parent_rs_; // index of parent reachable set
		int											 loc_;		 // index of the location

		int internalCount_;

	public:
		/****************************************************************************
   *                                                                           *
   *                           Constructors and Destructors                    *
   *                                                                           *
   *****************************************************************************/
		/**
   * @brief Constructor with no params
   */
		polyReachableSet();

		/**
   * @brief Constructor with TimePoint
   * @param time_point struct with fields .set and .time storing the time point
   * reachable set
   */
		explicit polyReachableSet(std::vector<polyReachableSetElement<Number>>& time_point);

		/**
   * @brief Constructor with two params
   * @param time_point struct with fields .set and .time storing the time point
   * reachable set
   * @param parent index of the parent reachable set
   */
		polyReachableSet(std::vector<polyReachableSetElement<Number>>& time_point, int parent);

		/**
   * @brief Constructor with two params
   * @param time_point struct with fields .set and .time storing the time point
   * reachable set
   * @param time_interval struct with fields .set, .time, and .algebraic
   * (nonlinDAsys) storing the time interval reachable set
   */
		polyReachableSet(std::vector<polyReachableSetElement<Number>>& time_point,
						 std::vector<polyReachableSetElement<Number>>& time_interval);

		/**
   * @brief Constructor with three params
   * @param time_point struct with fields .set and .time storing the time point
   * reachable set
   * @param parent index of the parent reachable set
   * @param loc index of the location
   */
		polyReachableSet(std::vector<polyReachableSetElement<Number>>& time_point, int parent, int loc);

		/**
   * @brief Constructor with three params
   * @param time_point struct with fields .set and .time storing the time point
   * reachable set
   * @param time_interval struct with fields .set, .time, and .algebraic
   * (nonlinDAsys) storing the time interval reachable set
   * @param parent index of the parent reachable set
   */
		polyReachableSet(std::vector<polyReachableSetElement<Number>>& time_point,
						 std::vector<polyReachableSetElement<Number>>& time_interval,
						 int										   parent);

		/**
   * @brief Constructor with four params
   * @param time_point struct with fields .set and .time storing the time point
   * reachable set
   * @param time_interval struct with fields .set, .time, and .algebraic
   * (nonlinDAsys) storing the time interval reachable set
   * @param parent index of the parent reachable set
   * @param loc index of the location
   */
		polyReachableSet(std::vector<polyReachableSetElement<Number>>& time_point,
						 std::vector<polyReachableSetElement<Number>>& time_interval,
						 int										   parent,
						 int										   loc);

		/**
   * @brief Copy Constructor - constructs a polyReachableSet from an existing
   * one.
   * @param other Another polyReachableSet, from which a new polyReachableSet is
   * constructed
   */
		polyReachableSet(const polyReachableSet<Number>& other) = default;

		virtual ~polyReachableSet();

		/*****************************************************************************
   *                                                                           *
   *                       Public Functions on Properties                      *
   *                                                                           *
   *****************************************************************************/

		/**
   * @brief Get the time_point
   * @return time_point
   */
		const std::vector<polyReachableSetElement<Number>> time_point() const;

		/**
   * @brief Replaces the current time_point with the parameter
   * @param time_point
   */
		void set_time_point(std::vector<polyReachableSetElement<Number>>& time_point);
      void set_time_point_rs(int i, polyZonotope<Number>& polyZono);
		/**
   * @brief Get the time_interval
   * @return time_interval
   */
		const std::vector<polyReachableSetElement<Number>> time_interval() const;

		/**
   * @brief Replaces the current time_interval with the parameter
   * @param time_interval
   */
		void set_time_interval(std::vector<polyReachableSetElement<Number>>& time_interval);

		/**
   * @brief Get the result time_point
   * @return time_point
   */
		polyTimePoint<Number> Rtime_point() const;

		/**
   * @brief Replaces the result time_point with the parameter
   * @param time_point
   */
		void set_Rtime_point(polyTimePoint<Number>& time_point);

		/**
   * @brief Get the result time_interval
   * @return time_interval
   */
		polyTimeInt<Number> time_Rinterval() const;

		/**
   * @brief Replaces the result time_interval with the parameter
   * @param time_interval
   */
		void set_Rtime_interval(polyTimeInt<Number>& time_interval);

		/**
   * @brief Get the time_interval
   * @return time_interval
   */
		const std::vector<polyReachableSetElement<Number>> R0() const;

		/**
   * @brief Replaces the current R0 with the parameter
   * @param R0
   */
		void set_R0(std::vector<polyReachableSetElement<Number>>& R0);

		/**
   * @brief Get the parent_rs
   * @return parent_rs
   */
		const int parent_rs() const;

		/**
   * @brief Replaces the current parent_rs with the parameter
   * @param parent_rs
   */
		void set_parent_rs(int parent_rs);

		/**
   * @brief Get the loc
   * @return loc
   */
		const int loc() const;

		/**
   * @brief Replaces the current loc with the parameter
   * @param loc
   */
		void set_loc(int loc);

		/**
   * @brief Get the internalCount
   * @return internalCount
   */
		const int internalCount() const;

		/**
   * @brief Replaces the current internalCount with the parameter
   * @param internalCount
   */
		void set_internalCount(int internalCount);

		/*****************************************************************************
   *                                                                           *
   *                       Other Functions                                     *
   *                                                                           *
   *****************************************************************************/

		/**
   * @brief delete reachable sets that are already covered by other sets
   * @param Rold reachable sets of previous time steps
   * @param options options for the computation of the reachable set
   */
		void deleteRedundantSets(polyReachableSet<Number> Rold, polyReachOptions<Number>& options);
	};

	/** @} */

	/**
 * @class      ReachableSet
 * @brief      Class for ReachableSet.
 * @ingroup    structure
 * @{
 */
	template <typename Number>
	class LinearReachableSet
	{
	private:
		Zonotope<Number> time_point_;
		Zonotope<Number> time_interval_;

	public:
		/****************************************************************************
   *                                                                           *
   *                           Constructors and Destructors                    *
   *                                                                           *
   *****************************************************************************/
		/**
   * @brief Constructor with no params
   */
		LinearReachableSet();

		/**
   * @brief Constructor with two params
   * @param time_point struct with fields .set and .time storing the time point
   * reachable set
   * @param time_interval struct with fields .set, .time, and .algebraic
   * (nonlinDAsys) storing the time interval reachable set
   */
		LinearReachableSet(Zonotope<Number>& time_point, Zonotope<Number>& time_interval);

		/**
   * @brief Copy Constructor - constructs a ReachableSet from an existing one.
   * @param other Another ReachableSet, from which a new ReachableSet is
   * constructed
   */
		LinearReachableSet(const LinearReachableSet<Number>& other) = default;

		virtual ~LinearReachableSet();

		/*****************************************************************************
   *                                                                           *
   *                       Public Functions on Properties                      *
   *                                                                           *
   *****************************************************************************/

		/**
   * @brief Get the time_point
   * @return time_point
   */
		const Zonotope<Number> time_point() const;

		/**
   * @brief Replaces the current time_point with the parameter
   * @param time_point
   */
		void set_time_point(Zonotope<Number>& time_point);

		/**
   * @brief Get the time_interval
   * @return time_interval
   */
		const Zonotope<Number> time_interval() const;

		/**
   * @brief Replaces the current time_interval with the parameter
   * @param time_interval
   */
		void set_time_interval(Zonotope<Number>& time_interval);

		/*****************************************************************************
   *                                                                           *
   *                       Other Functions                                     *
   *                                                                           *
   *****************************************************************************/
	};

	template <typename Number>
	class polyLinearReachableSet
	{
	private:
		polyZonotope<Number> time_point_;
		polyZonotope<Number> time_interval_;

	public:
		/****************************************************************************
   *                                                                           *
   *                           Constructors and Destructors                    *
   *                                                                           *
   *****************************************************************************/
		/**
   * @brief Constructor with no params
   */
		polyLinearReachableSet();

		/**
   * @brief Constructor with two params
   * @param time_point struct with fields .set and .time storing the time point
   * reachable set
   * @param time_interval struct with fields .set, .time, and .algebraic
   * (nonlinDAsys) storing the time interval reachable set
   */
		polyLinearReachableSet(polyZonotope<Number>& time_point, polyZonotope<Number>& time_interval);

		/**
   * @brief Copy Constructor - constructs a ReachableSet from an existing one.
   * @param other Another ReachableSet, from which a new ReachableSet is
   * constructed
   */
		polyLinearReachableSet(const polyLinearReachableSet<Number>& other) = default;

		virtual ~polyLinearReachableSet();

		/*****************************************************************************
   *                                                                           *
   *                       Public Functions on Properties                      *
   *                                                                           *
   *****************************************************************************/

		/**
   * @brief Get the time_point
   * @return time_point
   */
		const polyZonotope<Number> time_point() const;

		/**
   * @brief Replaces the current time_point with the parameter
   * @param time_point
   */
		void set_time_point(polyZonotope<Number>& time_point);

		/**
   * @brief Get the time_interval
   * @return time_interval
   */
		const polyZonotope<Number> time_interval() const;

		/**
   * @brief Replaces the current time_interval with the parameter
   * @param time_interval
   */
		void set_time_interval(polyZonotope<Number>& time_interval);

		/*****************************************************************************
   *                                                                           *
   *                       Other Functions                                     *
   *                                                                           *
   *****************************************************************************/
	};

	/** @} */
} // namespace reachSolver
#ifdef SEPARATE_TEMPLATE
#include "../../src/dynamics/ReachableSet.tpp"
#endif
#endif // REACHABLESET_H