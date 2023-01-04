/*
 * @Author: changyuan zhao
 * @Date: 2021-11-18 22:24:53
 */

#include <dynamics/ReachableSetElement.h>

namespace reachSolver
{

	template class ReachableSetElement<double>;

	/*************************Class  ReachableSetElement************************/
	template <typename Number>
	ReachableSetElement<Number>::ReachableSetElement()
		: set_(Zonotope<Number>())
		, error_(Vector_t<Number>())
		, prev_(0)
		, parent_(-1)
	{}

	template <typename Number>
	ReachableSetElement<Number>::ReachableSetElement(Zonotope<Number> set, Vector_t<Number> error)
		: set_(set)
		, error_(error)
		, prev_(0)
		, parent_(-1)
	{}

	template <typename Number>
	Zonotope<Number> ReachableSetElement<Number>::rs() const
	{
		return set_;
	}

	template <typename Number>
	void ReachableSetElement<Number>::set_rs(Zonotope<Number>& set)
	{
		set_ = set;
	}

	template <typename Number>
	Vector_t<Number> ReachableSetElement<Number>::error() const
	{
		return error_;
	}

	template <typename Number>
	void ReachableSetElement<Number>::set_error(Vector_t<Number> error)
	{
		error_ = error;
	}

	template <typename Number>
	const int ReachableSetElement<Number>::prev() const
	{
		return prev_;
	}

	template <typename Number>
	void ReachableSetElement<Number>::set_prev(int prev)
	{
		prev_ = prev;
	}

	template <typename Number>
	const int ReachableSetElement<Number>::parent() const
	{
		return parent_;
	}

	template <typename Number>
	void ReachableSetElement<Number>::set_parent(int parent)
	{
		parent_ = parent;
	}

} // namespace reachSolver