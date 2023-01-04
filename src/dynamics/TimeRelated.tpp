/**
 * @file   TimeRelated.cpp
 * @brief  TimeRelated class
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

#include <dynamics/TimeRelated.h>

namespace reachSolver
{

	template class TimePoint<double>;
	template class TimeInt<double>;

	/*****************************************************************************
 *                                                                           *
 *                           TimePoint Class                                 *
 *                                                                           *
 *****************************************************************************/

	template <typename Number>
	TimePoint<Number>::TimePoint()
		: set_(std::vector<std::vector<ReachableSetElement<Number>>>())
		, time_(std::vector<double>())
	{}

	template <typename Number>
	TimePoint<Number>::TimePoint(int size)
		: set_(std::vector<std::vector<ReachableSetElement<Number>>>(size))
		, time_(std::vector<double>(size))
	{}

	template <typename Number>
	std::vector<ReachableSetElement<Number>> TimePoint<Number>::rs(int i) const
	{
		return set_[i];
	}

	template <typename Number>
	std::vector<std::vector<ReachableSetElement<Number>>> TimePoint<Number>::set()
	{
		return set_;
	}

	template <typename Number>
	void TimePoint<Number>::set_rs(int i, std::vector<ReachableSetElement<Number>> set)
	{
		set_[i] = set;
	}

	template <typename Number>
	const double TimePoint<Number>::time(int i) const
	{
		return time_[i];
	}

	template <typename Number>
	void TimePoint<Number>::set_time(int i, double time)
	{
		time_[i] = time;
	}

	/*****************************************************************************
 *                                                                           *
 *                           TimeInt Class                                 *
 *                                                                           *
 *****************************************************************************/

	template <typename Number>
	TimeInt<Number>::TimeInt()
		: set_(std::vector<std::vector<ReachableSetElement<Number>>>())
		, time_(std::vector<Interval_t<Number>>())
	{}

	template <typename Number>
	TimeInt<Number>::TimeInt(int size)
		: set_(std::vector<std::vector<ReachableSetElement<Number>>>(size))
		, time_(std::vector<Interval_t<Number>>(size))
	{}

	template <typename Number>
	const std::vector<ReachableSetElement<Number>> TimeInt<Number>::rs(int i) const
	{
		return set_[i];
	}

	template <typename Number>
	std::vector<std::vector<ReachableSetElement<Number>>> TimeInt<Number>::set()
	{
		return set_;
	}

	template <typename Number>
	void TimeInt<Number>::set_rs(int i, std::vector<ReachableSetElement<Number>> set)
	{
		set_[i] = set;
	}

	template <typename Number>
	const Interval_t<Number> TimeInt<Number>::time(int i) const
	{
		return time_[i];
	}

	template <typename Number>
	void TimeInt<Number>::set_time(int i, Interval_t<Number> time)
	{
		time_[i] = time;
	}

	/*****************************************************************************
 *                                                                           *
 *                           polyTimePoint Class                              *
 *                                                                           *
 *****************************************************************************/

	template <typename Number>
	polyTimePoint<Number>::polyTimePoint()
		: set_(std::vector<std::vector<polyReachableSetElement<Number>>>())
		, time_(std::vector<double>())
	{}

	template <typename Number>
	polyTimePoint<Number>::polyTimePoint(int size)
		: set_(std::vector<std::vector<polyReachableSetElement<Number>>>(size))
		, time_(std::vector<double>(size))
	{}

	template <typename Number>
	std::vector<polyReachableSetElement<Number>> polyTimePoint<Number>::rs(int i) const
	{
		return set_[i];
	}

	template <typename Number>
	std::vector<std::vector<polyReachableSetElement<Number>>> polyTimePoint<Number>::set()
	{
		return set_;
	}

	template <typename Number>
	void polyTimePoint<Number>::set_rs(int i, std::vector<polyReachableSetElement<Number>> set)
	{
		set_[i] = set;
	}

	template <typename Number>
	const double polyTimePoint<Number>::time(int i) const
	{
		return time_[i];
	}

	template <typename Number>
	void polyTimePoint<Number>::set_time(int i, double time)
	{
		time_[i] = time;
	}

	/*****************************************************************************
 *                                                                           *
 *                           polyTimeInt Class                                 *
 *                                                                           *
 *****************************************************************************/

	template <typename Number>
	polyTimeInt<Number>::polyTimeInt()
		: set_(std::vector<std::vector<polyReachableSetElement<Number>>>())
		, time_(std::vector<Interval_t<Number>>())
	{}

	template <typename Number>
	polyTimeInt<Number>::polyTimeInt(int size)
		: set_(std::vector<std::vector<polyReachableSetElement<Number>>>(size))
		, time_(std::vector<Interval_t<Number>>(size))
	{}

	template <typename Number>
	const std::vector<polyReachableSetElement<Number>> polyTimeInt<Number>::rs(int i) const
	{
		return set_[i];
	}

	template <typename Number>
	std::vector<std::vector<polyReachableSetElement<Number>>> polyTimeInt<Number>::set()
	{
		return set_;
	}

	template <typename Number>
	void polyTimeInt<Number>::set_rs(int i, std::vector<polyReachableSetElement<Number>> set)
	{
		set_[i] = set;
	}

	template <typename Number>
	const Interval_t<Number> polyTimeInt<Number>::time(int i) const
	{
		return time_[i];
	}

	template <typename Number>
	void polyTimeInt<Number>::set_time(int i, Interval_t<Number> time)
	{
		time_[i] = time;
	}

} // namespace reachSolver