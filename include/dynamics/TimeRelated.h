/**
 * @file   TimeRelated.h
 * @brief  TimeRelated class representation
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

#ifndef TIMERELATED_H
#define TIMERELATED_H

#include <dynamics/ReachableSetElement.h>
#include <dynamics/polyReachableSetElement.h>
#include <zonotope/BasicObject.h>
#include <zonotope/Zonotope.h>
#include <zonotope/polyZonotope.h>

namespace reachSolver
{

	/**
 * @class      TimePoint
 * @brief      Class for TimePoint.
 * @ingroup    structure
 * @{
 */
	template <typename Number>
	class TimePoint
	{
	private:
		std::vector<std::vector<ReachableSetElement<Number>>> set_;
		std::vector<double>									  time_;

	public:
		/*****************************************************************************
   *                                                                           *
   *                           Constructors and Destructors                    *
   *                                                                           *
   *****************************************************************************/

		/**
   * @brief Constructor with no params
   */
		TimePoint();

		/**
   * @brief Constructor with size
   */
		TimePoint(int size);

		/*****************************************************************************
   *                                                                           *
   *                       Public Functions on Properties                      *
   *                                                                           *
   *****************************************************************************/
		/**
   * @brief Get the set(i)
   * @return the set
   */
		std::vector<ReachableSetElement<Number>> rs(int i) const;

		/**
   * @brief Get the set
   * @return the set
   */
		std::vector<std::vector<ReachableSetElement<Number>>> set();

        /**
   * @brief Add set;
   */
		void addSet(std::vector<ReachableSetElement<Number>> set){
      set_.push_back(set);
    }

      /**
   * @brief Add time;
   */
		void addTime(double time){
      time_.push_back(time);
    }

		/**
   * @brief Replaces the current set with the parameter
   * @param set
   */
		void set_rs(int i, std::vector<ReachableSetElement<Number>> set);

		/**
   * @brief Get the time
   * @return the time
   */
		const double time(int i) const;

		/**
   * @brief Replaces the current time with the parameter
   * @param time
   */
		void set_time(int i, double time);
	};
	/** @} */

	/**
 * @class      TimeInt
 * @brief      Class for TimeInt.
 * @ingroup    structure
 * @{
 */
	template <typename Number>
	class TimeInt
	{
	private:
		std::vector<std::vector<ReachableSetElement<Number>>> set_;
		std::vector<Interval_t<Number>>						  time_;

	public:
		/*****************************************************************************
   *                                                                           *
   *                           Constructors and Destructors                    *
   *                                                                           *
   *****************************************************************************/

		/**
   * @brief Constructor with no params
   */
		TimeInt();

		/**
   * @brief Constructor with size
   */
		TimeInt(int size);

		/*****************************************************************************
   *                                                                           *
   *                       Public Functions on Properties                      *
   *                                                                           *
   *****************************************************************************/
		/**
   * @brief Get the set(i)
   * @return the set
   */
		const std::vector<ReachableSetElement<Number>> rs(int i) const;

    /**
   * @brief Get the set
   * @return the set
   */
		std::vector<std::vector<ReachableSetElement<Number>>> set();

    /**
   * @brief Add set;
   */
		void addSet(std::vector<ReachableSetElement<Number>> set){
      set_.push_back(set);
    }

      /**
   * @brief Add time;
   */
		void addTime(Interval_t<Number> time){
      time_.push_back(time);
    }


		/**
   * @brief Replaces the current set with the parameter
   * @param set
   */
		void set_rs(int i, std::vector<ReachableSetElement<Number>> set);

		/**
   * @brief Get the time
   * @return the time
   */
		const Interval_t<Number> time(int i) const;

		/**
   * @brief Replaces the current time with the parameter
   * @param time
   */
		void set_time(int i, Interval_t<Number> time);
	};
	/** @} */

	/**
 * @class      TimePoint
 * @brief      Class for TimePoint.
 * @ingroup    structure
 * @{
 */
	template <typename Number>
	class polyTimePoint
	{
	private:
		std::vector<std::vector<polyReachableSetElement<Number>>> set_;
		std::vector<double>										  time_;

	public:
		/*****************************************************************************
   *                                                                           *
   *                           Constructors and Destructors                    *
   *                                                                           *
   *****************************************************************************/

		/**
   * @brief Constructor with no params
   */
		polyTimePoint();

		/**
   * @brief Constructor with size
   */
		polyTimePoint(int size);

		/*****************************************************************************
   *                                                                           *
   *                       Public Functions on Properties                      *
   *                                                                           *
   *****************************************************************************/
		/**
   * @brief Get the set(i)
   * @return the set
   */
		std::vector<polyReachableSetElement<Number>> rs(int i) const;

		/**
   * @brief Get the set
   * @return the set
   */
		std::vector<std::vector<polyReachableSetElement<Number>>> set();

		/**
   * @brief Replaces the current set with the parameter
   * @param set
   */
		void set_rs(int i, std::vector<polyReachableSetElement<Number>> set);

		/**
   * @brief Get the time
   * @return the time
   */
		const double time(int i) const;

		/**
   * @brief Replaces the current time with the parameter
   * @param time
   */
		void set_time(int i, double time);
	};
	/** @} */

	/**
 * @class      polyTimeInt
 * @brief      Class for polyTimeInt.
 * @ingroup    structure
 * @{
 */
	template <typename Number>
	class polyTimeInt
	{
	private:
		std::vector<std::vector<polyReachableSetElement<Number>>> set_;
		std::vector<Interval_t<Number>>							  time_;

	public:
		/*****************************************************************************
   *                                                                           *
   *                           Constructors and Destructors                    *
   *                                                                           *
   *****************************************************************************/

		/**
   * @brief Constructor with no params
   */
		polyTimeInt();

		/**
   * @brief Constructor with size
   */
		polyTimeInt(int size);

		/*****************************************************************************
   *                                                                           *
   *                       Public Functions on Properties                      *
   *                                                                           *
   *****************************************************************************/
		/**
   * @brief Get the set(i)
   * @return the set
   */
		const std::vector<polyReachableSetElement<Number>> rs(int i) const;

    	/**
   * @brief Get the set
   * @return the set
   */
		std::vector<std::vector<polyReachableSetElement<Number>>> set();

		/**
   * @brief Replaces the current set with the parameter
   * @param set
   */
		void set_rs(int i, std::vector<polyReachableSetElement<Number>> set);

		/**
   * @brief Get the time
   * @return the time
   */
		const Interval_t<Number> time(int i) const;

		/**
   * @brief Replaces the current time with the parameter
   * @param time
   */
		void set_time(int i, Interval_t<Number> time);
	};

} // namespace reachSolver
#ifdef SEPARATE_TEMPLATE
#include "../../src/dynamics/TimeRelated.tpp"
#endif
#endif // TIMERELATED_H