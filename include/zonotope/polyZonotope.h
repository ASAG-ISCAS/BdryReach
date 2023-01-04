/**
 * @file   polyZonotope.h
 * @brief  polynomial Zonotope class representation
 * @author Changyuan Zhao
 * @date  Nov 2021
 * @version 1.1
 *
 * Reference:
 *   CORA ../contSet/polyZonotope/polyZonotope.m
 */

#ifndef POLYZONOTOPE_H
#define POLYZONOTOPE_H

#include <algorithm>
#include <iostream>
#include <zonotope/BasicObject.h>
#include <zonotope/Zonotope.h>

namespace reachSolver
{

	/**
 * @class      polyZonotope
 * @brief      Class for polynomial Zonotopes.
 * @tparam     Number    The used number type.
 * @ingroup    structure
 * @{
 */
	template <typename Number>
	class polyZonotope : private BasicObject
	{
	private:
		/**
   * @brief: dimension of the polynomial zonotope
   */
		size_t dimension_;

		/**
   * @brief: center vector of the polynomial zonotope
   */
		Vector_t<Number> center_;

		/**
   * @brief: generators of the polynomial zonotope
   */
		Matrix_t<Number> generators_;

		/**
   * @brief: expMat of the polynomial zonotope
   */
		Matrix_t<int> expMat_;

		/**
   * @brief: Grest of the polynomial zonotope
   */
		Matrix_t<Number> Grest_;

		/**
   * @brief: id of the polynomial zonotope
   */
		std::vector<size_t> id_;

	public:
		/*****************************************************************************
   *                                                                           *
   *                           Constructors and Destructors                    *
   *                                                                           *
   *****************************************************************************/

		/**
   * @brief Constructor with no params
   */
		polyZonotope();

		/**
   * @brief Constructor with dimension
   * @param dimension Dimensionality of Zonotope
   */
		explicit polyZonotope(size_t dimension);

		/**
   * @brief Constructor with center and generators.
   * @param center A vector
   * @param generators A matrix
   * @param expMat A matrix
   * @param Grest A matrix
   */
		polyZonotope(const Vector_t<Number>& center,
					 const Matrix_t<Number>& generators,
					 const Matrix_t<int>&	 expMat,
					 const Matrix_t<Number>& Grest);

		/**
   * @brief Constructor with center and generators.
   * @param center A vector
   * @param generators A matrix
   * @param expMat A matrix
   * @param Grest A matrix
   * @param id A vector
   */
		polyZonotope(const Vector_t<Number>&   center,
					 const Matrix_t<Number>&   generators,
					 const Matrix_t<int>&	   expMat,
					 const Matrix_t<Number>&   Grest,
					 const std::vector<size_t> id);

		/**
   * @brief Copy Constructor - constructs a zonotope from an existing one.
   * @param other Another Zonotope, from which a new zonotope is constructed
   */
		polyZonotope(const polyZonotope& other) = default;

		/**
   * @brief Constructor with a zonotope.
   * @param zono A zonotope
   * */
		polyZonotope(const Zonotope<Number> zono);

		/**
   * @brief: destructor
   */
		virtual ~polyZonotope();

		/*****************************************************************************
   *                                                                           *
   *                       Public Functions on Properties                      *
   *                                                                           *
   *****************************************************************************/

		/**
   * @brief Get the dimension of polynomial Zonotope
   * @return the dimension
   */
		size_t dimension() const;

		/**
   * @brief Number of generators
   * @return number of generators
   */
		size_t numGenerators() const;

		/**
   * @brief Get the current center
   * @return center a nx1 matrix
   */
		const Vector_t<Number>& center() const;

		/**
   * @brief Replaces the current center with the parameter center
   * @param center a nx1 matrix
   */
		void set_center(const Vector_t<Number>& center);

		/**
   * @brief Get the current generators
   * @return center a nxm matrix
   */
		const Matrix_t<Number>& generators() const;

		/**
   * @brief Replaces the current matrix of generators with the parameter
   * generators
   * @param generators a nxm matrix
   */
		void set_generators(const Matrix_t<Number>& generators);

		/**
   * @brief Get the current expMat
   * @return center a nxm matrix
   */
		const Matrix_t<int>& expMat() const;

		/**
   * @brief Replaces the current matrix of expMat with the parameter expMat
   * @param generators a nxm matrix
   */
		void set_expMat(const Matrix_t<int>& expMat);

		/**
   * @brief Get the current expMat
   * @return center a nxm matrix
   */
		const Matrix_t<Number>& Grest() const;

		/**
   * @brief Replaces the current matrix of expMat with the parameter expMat
   * @param generators a nxm matrix
   */
		void set_Grest(const Matrix_t<Number>& Grest);

		/**
   * @brief Get the current expMat
   * @return center a nxm matrix
   */
		const std::vector<size_t>& id() const;

		/**
   * @brief Replaces the current matrix of expMat with the parameter expMat
   * @param generators a nxm matrix
   */
		void set_id(const std::vector<size_t>& id);

		/**
   * @brief Removes zero generators in generator matrix
   */
		void DeleteZeroGenerators();

		/**
   * @brief Changes the dimension of a Zonotope. if new_dim > old dim, new rows
   * are initialized with null
   * @param new_dim The new dimension of the Zonotope
   * @return True, if change in dimension was successful
   */
		bool ChangeDimension(size_t new_dim);

		/**
   * @brief: remove specified column Grest
   * @param colomn a given column
   * @return the zonotope without the specified column Grest
   */
		void RemoveGrest(unsigned int colomn);

		/**
   * @brief: remove specified column generators
   * @param colomn a given column
   * @return the zonotope without the specified column generators
   */
		void Removegenerators(unsigned int colomn);

		/**
   * @brief: remove specified column expMat
   * @param colomn a given column
   * @return the zonotope without the specified column expMat
   */
		void RemoveexpMat(unsigned int colomn);

		/**
   * @brief Reduces the order of a zonotope
   * @param limitOrder order of reduced zonotope
   */
		void Reduce(unsigned int limitOrder);

		/**
   * @brief computes an enclosing zonotope of the polynomial zonotope
   * @return zonotope
   */
		Zonotope<Number> toZonotope() const;

		/**
   * @brief Whether the Zonotope is empty
   * @return 0 if empty, 1 if not empty
   */
		int Empty();

		/**
   * @brief Clears the generators and center of the Zonotope and sets
   * dimensionality to zero
   */
		void Clear();

		/**
   * @brief display the zonotope
   */
		void Display() const;

		/**
   * @brief Judge whether two zonotope is equal
   * @param another_zonotope
   * @return True, if equal
   */
		bool operator==(const Zonotope<Number>& another_zonotope) const
		{
			if (this->dimension_ != another_zonotope.dimension())
			{
				return false;
			}
			if (this->center_ != another_zonotope.center())
			{
				return false;
			}
			if (this->generators_ != another_zonotope.generators())
			{
				return false;
			}
			return true;
		}

		/*****************************************************************************
   *                                                                           *
   *                          Basic Set Operations                             *
   *                                                                           *
   *****************************************************************************/

		/**
   * @brief implement the linear maps of a set, i.e., "*" operator
   * @param num a number
   * @return a zonotope = num * this zonotope
   */
		polyZonotope operator*(Number num) const;

		/**
   * @brief implement the linear maps of a set, i.e., "*" operator
   * @param matrix a matrix
   * @return a  zonotope = matrix * this zonotope
   */
		polyZonotope Times(const Matrix_t<Number>& matrix) const;

		/**
   * @brief implement the linear maps of a set, i.e., "*" operator
   * @param matrix a matrix
   * @return a  zonotope = matrix * this zonotope
   */
		polyZonotope operator*(const Matrix_t<Number>& matrix) const;

		/**
   * @brief implement the times of zonotope and intervalmatrix
   * @param int_matrix an IntervalMatrix
   * @return a  zonotope = int_matrix * this zonotope
   */
		polyZonotope Times(const IntervalMatrix<Number> int_matrix) const;

		/**
   * @brief implement the times of zonotope and intervalmatrix
   * @param int_matrix an IntervalMatrix
   * @return a  zonotope = int_matrix * this zonotope
   */
		polyZonotope operator*(const IntervalMatrix<Number> int_matrix) const;

		/**
   * @brief Get the minkowski addition of two zonotope,i.e., "+" operator
   * @param another_zonotope
   * @return a  zonotope = zonotope1 + zonotope2
   */
		polyZonotope Plus(const Zonotope<Number>& another_zonotope) const;

		/**
   * @brief Get the addition of two polyzonotopes
   * @param another_polyzonotope
   * @return a  polyzonotope = polyzonotope1 + polyzonotope2
   */
		polyZonotope Plus(const polyZonotope& another_polyzonotope) const;

		/**
   * @brief Get the addition of a polyzonotope and an intervalMatrix
   * @param inval
   * @return the sum
   */
		polyZonotope operator+(const IntervalMatrix<Number> inval) const;
		/**
   * @brief Get the minkowski addition of two zonotope,i.e., "+" operator
   * @param another_zonotope
   * @return the sum
   */
		polyZonotope operator+(const Zonotope<Number>& another_zonotope) const;

		/**
   * @brief Get the addition of two polyzonotopes
   * @param another_polyzonotope
   * @return a  polyzonotope = polyzonotope1 + polyzonotope2
   */
		polyZonotope operator+(const polyZonotope& another_polyzonotope) const;

		/**
   * @brief Get the addition of a zonotope and a vector
   * @param vector
   * @return the sum
   */
		polyZonotope operator+(const Vector_t<Number>& vector) const;

		/**
   * @brief Get the addition of a zonotope and a vector
   * @param num
   * @return the sum
   */
		polyZonotope operator+(const Number num) const;

		/**
   * @brief Get the subtraction of a zonotope'center and a num
   * @param num
   * @return the sum
   */
		polyZonotope operator-(const Number num) const;

		/**
   * @brief Get the subtraction of a zonotope'center and a num
   * @param vector
   * @return the sum
   */
		polyZonotope operator-(const Vector_t<Number>& vector) const;
		// /**
		//  * @brief overloads & operator, computes the intersection of two zonotopes
		//  * @param another_zonotope
		//  * @return the and
		//  */
		// Zonotope operator&(const Zonotope<Number>& another_zonotope) const;

		/**
   * @brief Get the enclosure for the convex hull of two zonotope
   * @param another_zonotope
   * @return a zonotope enclosing the convex hull
   */
		// Zonotope Convexhull(const Zonotope<Number>& another_zonotope) const;

		/**
   * @brief Generates a polyzonotope that encloses a polyzonotopes and its
   * linear transformation
   * @param another_zonotope second polyzonotope object, satisfying Z2 = (M *
   * Z1) + Zplus
   * @return zonotope, that encloses Z1 and Z2
   */
		polyZonotope enclose(const polyZonotope<Number>& another_polyzonotope) const;

		/**
   * @brief: return the cartesian produce of two polyZonotope object
   * @param pZ2 contSet object
   * @return polyzonotope
   */
		polyZonotope cartProd(const polyZonotope<Number> pZ2);

		polyZonotope cartProd(const Zonotope<Number> pZ2);

		/**
   * @brief: compute the quadratic map of a polyzonotope
   * @param  H  quadratic coefficients as a cell of interval matrices
   * @return resulting set as a ployZonotope object
   */
		polyZonotope quadMap(std::vector<Matrix_t<Number>> H);

		/**
   * @brief Overapproximates a zonotope by an interval hull
   * @return interval object
   */
		IntervalMatrix<Number> interval();

		/**
   * @brief: compute the addition of two sets while preserving the dependencies
   * between the two sets
   * @param polyZonotope pZ2
   * @return result of exact addition
   */
		polyZonotope exactPlus(polyZonotope pZ2);

		// friend std::ostream & operator << (std::ostream &,const Zonotope<Number>
		// &);
	};

    /**
   * @brief: deleteZeros - deletes all generators of length 0
   */
    template <typename Number>
	void deleteZeros(polyZonotope<Number> & polyZono);

    /**
   * @brief: computes an enclosing zonotope of the polynomial zonotope
   * @param polyZonotope polyZono
   * @return zonotope object 
   */
    template <typename Number>
	Zonotope<Number> zonotope(polyZonotope<Number> & polyZono);

    /**
   * @brief: Calculate a new over-approxmiating representation of a
   * polynomial zonotope in such a way that there remain no
   * independent generators
   * @param polyZonotope polyZono
   * @return resuling polyZonotope object which over-approximates pZ 
   */
    template <typename Number>
	polyZonotope<Number> restructure(polyZonotope<Number> polyZono, unsigned int limitOrder);

    /**
   * @brief: polyshape objects are created from the x- and 
   * y-coordinates that describe the polygon's boundaries.  
   * @param x- and y-coordinates
   * @return the polygon's boundaries
   */
    template <typename Number>
	Matrix_t<Number> polyshape(Matrix_t<Number>& matrix);

    /**
   * @brief: enclose polynomial zonotope with a polygon 
   */
    template <typename Number>
	void getPolygon(polyZonotope<Number> const & pZ, Matrix_t<Number>& V, Matrix_t<Number>& poly);

    /**
   * @brief: Splits the longest generator dependent generator with a 
   * polynomial order of 1 for a polynomial zonotope 
   */
    template <typename Number>
	void splitLongestGen(polyZonotope<Number> const & pZ, polyZonotope<Number>& pZ1, polyZonotope<Number>& pZ2);

    /**
   * @brief: check if a polynomial zonotope is full-dimensional
   */
    template <typename Number>
	bool isFullDim(polyZonotope<Number> const & pZ);

    /**
   * @brief: Returns a polyZonotope which is projected onto the specified dimensions
   */
    template <typename Number>
	polyZonotope<Number> project(polyZonotope<Number> const & pZ, int dim);

} // namespace reachSolver
#ifdef SEPARATE_TEMPLATE
#include "../../src/zonotope/polyZonotope.tpp"
#endif
#endif // POLYZONOTOPE_H