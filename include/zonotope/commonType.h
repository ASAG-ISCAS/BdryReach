/*
 * @Author: your name
 * @Date: 2021-11-12 16:24:21
 * @LastEditTime: 2021-12-24 18:33:38
 * @LastEditors: Please set LastEditors
 * @Description: 打开koroFileHeader查看配置 进行设置:
 * https://github.com/OBKoro1/koro1FileHeader/wiki/%E9%85%8D%E7%BD%AE
 * @FilePath: /Solver/include/zonotope/commonType.h
 */
/**
 * @file   commonType.h
 * @brief
 * @author Yunjun Bai
 * @date April 2021
 * @version 1.0
 *
 */

/**
 * @author Changyuan Zhao
 * @date  Nov 2021
 * @version 1.1
 */

#ifndef COMMONTYPE_H
#define COMMONTYPE_H

#include <capd/capdlib.h>
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/StdVector>
#include <eigen3/unsupported/Eigen/MatrixFunctions>
#include <iostream>
#include <set>

namespace reachSolver
{

	/// Wrapper for an Eigen::Matrix type with only one column.
	template <typename Number>
	using Vector_t = Eigen::Matrix<Number, Eigen::Dynamic, 1>;

	/// Wrapper for an Eigen::Matrix type.
	template <typename Number>
	using Matrix_t = Eigen::Matrix<Number, Eigen::Dynamic, Eigen::Dynamic>;

	template <typename Number>
	using RowVector_t = Eigen::Matrix<Number, 1, Eigen::Dynamic>;

	/// Wrapper Eigen::Matrix type with only one column as a set
	template <typename Number>
	using Vector_set_t = std::set<Vector_t<Number>>;

	/// Wrapper capd::interval
	template <typename Number>
	using Interval_t = capd::intervals::Interval<Number>;

	template <typename Number>
	// using IntervalMatrix = Eigen::Matrix<capd::intervals::Interval<Number>,
	// Eigen::Dynamic,Eigen::Dynamic>; using IntervalMatrix =
	// Matrix_t<capd::intervals::Interval<Number>>;
	using IntervalMatrix = Matrix_t<Interval_t<Number>>;

	/**
 * @brief: a struct of the interbal matrix
 */
	// template <typename Number>
	// struct IntervalMatrix
	//{
	/**
 * @brief: the inf of an interbal matrix
 */
	// Matrix_t<Interval<Number>> inf;

	/**
 * @brief: the sup of an interbal matrix
 */
	// Matrix_t<Number> sup;
	//};

	/**
 * @brief: a struct of the interbal
 */
	/*struct Interval
{
    /**
     * @brief: an array with interval's inf and sup
     */
	// double interval[2];
	//};

	/**
 * @brief:
 */
	struct linerror_p_type
	{
		Vector_t<double> x;
		Vector_t<double> u;
	};

	/**
 * @brief:
 */
	struct linerror_type
	{
		Vector_t<double> f0;
		linerror_p_type	 p;
	};

	/**
 * @brief: reachability
 */
	enum class REACHABILITY_RESULT
	{
		SAFE,
		UNKNOWN
	};

} // namespace reachSolver
#endif // COMMONTYPE_H