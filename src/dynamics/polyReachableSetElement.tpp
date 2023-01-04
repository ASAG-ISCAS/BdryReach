/*
 * @Author: changyuan zhao
 * @Date: 2021-11-18 22:24:53
 */
#include <dynamics/polyReachableSetElement.h>

namespace reachSolver
{

	template class polyReachableSetElement<double>;

	/*************************Class polyReachableSetElement************************/
	template <typename Number>
	polyReachableSetElement<Number>::polyReachableSetElement()
		: set_(polyZonotope<Number>())
		, error_(Vector_t<Number>())
		, prev_(0)
		, parent_(0)
	{}

	template <typename Number>
	polyReachableSetElement<Number>::polyReachableSetElement(polyZonotope<Number> set, Vector_t<Number> error)
		: set_(set)
		, error_(error)
		, prev_(0)
		, parent_(0)
	{}

	template <typename Number>
	polyZonotope<Number> polyReachableSetElement<Number>::rs() const
	{
		return set_;
	}

	template <typename Number>
	void polyReachableSetElement<Number>::set_rs(polyZonotope<Number>& set)
	{
		set_ = set;
	}

	template <typename Number>
	Vector_t<Number> polyReachableSetElement<Number>::error() const
	{
		return error_;
	}

	template <typename Number>
	void polyReachableSetElement<Number>::set_error(Vector_t<Number> error)
	{
		error_ = error;
	}

	template <typename Number>
	const int polyReachableSetElement<Number>::prev() const
	{
		return prev_;
	}

	template <typename Number>
	void polyReachableSetElement<Number>::set_prev(int prev)
	{
		prev_ = prev;
	}

	template <typename Number>
	const int polyReachableSetElement<Number>::parent() const
	{
		return parent_;
	}

	template <typename Number>
	void polyReachableSetElement<Number>::set_parent(int parent)
	{
		parent_ = parent;
	}

} // namespace reachSolver