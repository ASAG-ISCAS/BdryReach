/*
 * @Author: your name
 * @Date: 2021-11-11 20:10:50
 * @LastEditTime: 2021-12-22 22:06:38
 * @LastEditors: Please set LastEditors
 * @Description: 打开koroFileHeader查看配置 进行设置:
 * https://github.com/OBKoro1/koro1FileHeader/wiki/%E9%85%8D%E7%BD%AE
 * @FilePath: /Solver/src/dynamics/ContDynamics.cpp
 */
/**
 * @file   ContDynamics.cpp
 * @brief  ContDynamics class
 * @author Yunjun Bai
 * @date April 2021
 * @version 1.0
 *
 * Reference:
 * CORA ../contDynamics/contDynamics/all
 */

/**
 * @author Changyuan Zhao
 * @date  Nov 2021
 * @version 1.1
 *
 */

#include <dynamics/ContDynamics.h>

namespace reachSolver
{

	template class ContDynamics<double>;

	template <typename Number>
	ContDynamics<Number>::ContDynamics(std::string name)
	{
		name_ = name;
	}

	template <typename Number>
	ContDynamics<Number>::ContDynamics(std::string name, size_t dim, size_t num_inputs, size_t num_outputs)
	{
		name_		 = name;
		dim_		 = dim;
		num_inputs_	 = num_inputs;
		num_outputs_ = num_outputs;
	}

	template <typename Number>
	ContDynamics<Number>::~ContDynamics()
	{}

	template <typename Number>
	const std::string ContDynamics<Number>::name() const
	{
		return name_;
	}

	template <typename Number>
	const size_t ContDynamics<Number>::num_outputs() const
	{
		return num_outputs_;
	}

	template <typename Number>
	const size_t ContDynamics<Number>::dim() const
	{
		return dim_;
	}

	template <typename Number>
	const size_t ContDynamics<Number>::num_inputs() const
	{
		return num_inputs_;
	}

	template <typename Number>
	void ContDynamics<Number>::set_num_outputs(const size_t num_outputs)
	{
		num_outputs_ = num_outputs;
	}

	template <typename Number>
	void ContDynamics<Number>::set_dim(const size_t dim)
	{
		dim_ = dim;
	}

	template <typename Number>
	void ContDynamics<Number>::set_num_inputs(const size_t num_inputs)
	{
		num_inputs_ = num_inputs;
	}

} // namespace reachSolver
