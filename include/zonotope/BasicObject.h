/**
 * @file   BasicObject.h
 * @brief  Basic class of these set representations
 * @author Yunjun Bai
 * @date April 2021
 * @version 1.0
 *
 * @defgroup structure
 * Reference:
 *   CORA ../contSet/contSet/contSet.m
 */

/**
 * @author Changyuan Zhao
 * @date  Nov 2021
 * @version 1.1
 */

#ifndef BASICOBJECT_H
#define BASICOBJECT_H

#include <cstddef>
#include <iostream>
#include <zonotope/commonType.h>

namespace reachSolver {

    /**
 * @class   BasicObject
 * @brief  basic objects
 * @ingroup structure
 */
    class BasicObject {

    protected:
        /**
   * @brief Constructor with no params
   */
        BasicObject() = default;

        /**
   * @brief: Constructor with a const input
   * @param in a const input
   */
        BasicObject(const BasicObject &in) = default;

        /**
   * @brief: Constructor with a input
   * @param in a input
   */
        BasicObject(BasicObject &in) = default;

        /**
   * @brief assign the specified input
   * @param in a const input
   * @return a basicobject with given input
   */
        BasicObject &operator=(const BasicObject &in) = default;

        /**
   * @brief assign the specified input
   * @param in a input
   * @return a basicobject with given input
   */
        BasicObject &operator=(BasicObject &in) = default;
    };
    /** @} */

}// namespace reachSolver
#endif// BASICOBJECT_h