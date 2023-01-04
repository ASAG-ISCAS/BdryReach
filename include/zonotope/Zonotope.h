/**
 * @file   Zonotope.h
 * @brief  Zonotope class representation
 * @author Yunjun Bai
 * @date April 2021
 * @version 1.0
 *
 * Reference:
 *   CORA ../contSet/zonotope/zonotope.m
 */

/**
 * @author Changyuan Zhao
 * @date  Nov 2021
 * @version 1.1
 */
#ifndef ZONOTOPE_H
#define ZONOTOPE_H

#include <iostream>
#include <zonotope/BasicObject.h>
#include <zonotope/globalfunctions.h>
#include <cmath>
#include <eigen3/Eigen/Core>

using namespace std;

namespace reachSolver
{

	/**
 * @class      Zonotope
 * @brief      Class for Zonotopes.
 * @tparam     Number     The used number type.
 * @ingroup    structure
 * @{
 */
	template <typename Number>
	class Zonotope : private BasicObject
	{
	private:
		/**
   * @brief: dimension of the zonotope
   */
		size_t dimension_;

		/**
   * @brief: center vector of the zonotope
   */
		Vector_t<Number> center_;

		/**
   * @brief: generators of the zonotope
   */
		Matrix_t<Number> generators_;

		/**
   * @brief Remove a generator
   * @param colomn column of generator
   */
		void RemoveGenerator(unsigned int colomn);

	public:
		/*****************************************************************************
   *                                                                           *
   *                           Constructors and Destructors                    *
   *                                                                           *
   *****************************************************************************/

		/**
   * @brief Constructor with no params
   */
		Zonotope();

		/**
   * @brief Constructor with dimension
   * @param dimension Dimensionality of Zonotope
   */
		explicit Zonotope(size_t dimension);

		/**
   * @brief Constructor with interval_m
   * @param interval_m IntervalMatrix
   */
		explicit Zonotope(IntervalMatrix<Number> interval_m);

		/**
   * @brief Constructor with center and generators.
   * @param center A vector
   * @param generators A matrix
   */
		Zonotope(const Vector_t<Number>& center, const Matrix_t<Number>& generators);

		/**
   * @brief Constructor with center and generators.
   * @param center A vector
   */
		Zonotope(const Vector_t<Number>& center);

		/**
   * @brief Copy Constructor - constructs a zonotope from an existing one.
   * @param other Another Zonotope, from which a new zonotope is constructed
   */
		Zonotope(const Zonotope& other) = default;

		/**
   * @brief: destructor
   */
		virtual ~Zonotope();

		/*****************************************************************************
   *                                                                           *
   *                       Public Functions on Properties                      *
   *                                                                           *
   *****************************************************************************/

		/**
   * @brief Get the dimension of Zonotope
   * @return the dimension
   */
		size_t dimension() const;

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
   * @brief Add generators to Zonotope. Simply performs setGenerators if
   * generators was previously not initialized.
   * @param generators a nxm matrix
   * @return true if able to add generators
   */
		bool AddGenerators(const Matrix_t<Number>& generators);

		/**
   * @brief Get the order
   * @return zonotope order
   */
		Number order() const;

		/**
   * @brief Number of generators
   * @return number of generators
   */
		size_t numGenerators() const;

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
   * @brief Reduces the order of a zonotope
   * @param limitOrder order of reduced zonotope
   */
		void Reduce(unsigned limitOrder);

        /**
    * @brief Reduces the order of a zonotope using PCA
    * @param limitOrder order of reduced zonotope
    */
	    void ReducePCA(unsigned limitOrder);    

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
		Zonotope operator*(Number num) const;

		/**
   * @brief implement the linear maps of a set, i.e., "*" operator
   * @param matrix a matrix
   * @return a  zonotope = matrix * this zonotope
   */
		Zonotope Times(const Matrix_t<Number>& matrix) const;

		/**
   * @brief implement the linear maps of a set, i.e., "*" operator
   * @param matrix a matrix
   * @return a  zonotope = matrix * this zonotope
   */
		Zonotope operator*(const Matrix_t<Number>& matrix) const;

		/**
   * @brief implement the times of zonotope and intervalmatrix
   * @param int_matrix an IntervalMatrix
   * @return a  zonotope = int_matrix * this zonotope
   */
		Zonotope Times(const IntervalMatrix<Number> int_matrix) const;

		/**
   * @brief implement the times of zonotope and intervalmatrix
   * @param int_matrix an IntervalMatrix
   * @return a  zonotope = int_matrix * this zonotope
   */
		Zonotope operator*(const IntervalMatrix<Number> int_matrix) const;

		/**
   * @brief Get the minkowski addition of two zonotope,i.e., "+" operator
   * @param another_zonotope
   * @return a  zonotope = zonotope1 + zonotope2
   */
		Zonotope Plus(const Zonotope<Number>& another_zonotope) const;

		/**
   * @brief Get the addition of a zonotope and a vector
   * @param vector
   * @return the sum
   */
		// Zonotope Plus(const Vector_t<Number>& vector) const;

		/**
   * @brief Get the minkowski addition of two zonotope,i.e., "+" operator
   * @param another_zonotope
   * @return the sum
   */
		Zonotope operator+(const Zonotope<Number>& another_zonotope) const;

		/**
   * @brief Get the addition of a zonotope and a vector
   * @param vector
   * @return the sum
   */
		Zonotope operator+(const Vector_t<Number>& vector) const;

		/**
   * @brief Get the addition of a zonotope and a vector
   * @param num
   * @return the sum
   */
		Zonotope operator+(const Number num) const;

		/**
   * @brief Get the subtraction of a zonotope'center and a num
   * @param num
   * @return the sum
   */
		Zonotope operator-(const Number num) const;

		/**
   * @brief Get the subtraction of a zonotope'center and a num
   * @param vector
   * @return the sum
   */
		Zonotope operator-(const Vector_t<Number>& vector) const;
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
		Zonotope Convexhull(const Zonotope& another_zonotope) const;

		/**
   * @brief Generates a zonotope that encloses a zonotopes and its linear
   * transformation
   * @param another_zonotope second zonotope object, satisfying Z2 = (M * Z1) +
   * Zplus
   * @return zonotope, that encloses Z1 and Z2
   */
		Zonotope enclose(const Zonotope& another_zonotope) const;

		/**
   * @brief Overapproximates a zonotope by an interval hull
   * @return interval object
   */
		IntervalMatrix<Number> interval() const;

		/**
   * @brief: return the cartesian produce of two Zonotope object
   * @param Z2 zonotope object
   * @return polyzonotope
   */
		Zonotope cartProd(const Zonotope<Number> Z2);

		/**
   * @brief: compute the quadratic map of a zonotope
   * @param  H  quadratic coefficients as a cell of interval matrices
   * @return resulting set as a ployZonotope object
   */
		Zonotope quadMap(std::vector<Matrix_t<Number>> Q);

		/**
   * @brief: compute the quadratic map of a zonotope
   * @param  H  quadratic coefficients as a cell of interval matrices
   * @param Z2 a Zonotope
   * @return resulting set as a Zonotope object
   */
		Zonotope quadMapMixed(Zonotope<Number> Z2, std::vector<Matrix_t<Number>> Q);

		/**
   * @brief: center vector and generator matrix Z = [c,g]
   * @return center vector and generator matrix Z = [c,g]
   */
		Matrix_t<Number> Z();

		/**
   * @brief: splits azonotope into two or more enclosinb zonotope
   * @return cell array of parallelpipeds represented as zonotopes
   */
		std::vector<std::vector<Zonotope<Number>>> split();
		std::vector<Zonotope<Number>>			   split(int dim);
		std::vector<Zonotope<Number>>			   splitOneDim(Zonotope<Number> Z, int dim);

		// friend std::ostream & operator << (std::ostream &,const Zonotope<Number>
		// &);
	};

	template <typename Number>
	std::ostream& operator<<(std::ostream&, const Zonotope<Number>&);
	/** @} */

    //删除矩阵某一列
    template <typename Number>
	void RemoveColumn(Matrix_t<Number> & matrix, unsigned int colToRemove);

    //删除矩阵某一行
    template <typename Number>
	void RemoveRow(Matrix_t<Number> & matrix, unsigned int rowToRemove);

    /**
   * @brief: delete zero columns in matrix
   * @return matrix with no zero column
   */
    template <typename Number>
	void deleteZeros(Matrix_t<Number> & matrix);

    /**
   * @brief: converts a two-dimensional zonotope into a polygon and returns
   * @return ordered set of points of a polygon
   */
    template <typename Number>
    Matrix_t<Number> convert2polygon(const Zonotope<Number>&);

        /**
   * @brief: Selects generators to be reduced
   * @return c - center of reduced zonotope
   * Gunred - generators that are not reduced 
   * Gred - generators that are reduced
   */
    template <typename Number>
	void pickedGenerators(Zonotope<Number>& Z, 
                          unsigned order, 
                          Vector_t<Number>& c, 
                          Matrix_t<Number>& Gunred,
                          Matrix_t<Number>& Gred);

    /**
   * @brief: Bounds of zonotope
   */
    template <typename Number>
	std::vector<Zonotope<Number>> Zono_border(const Zonotope<Number>& Zono);

    /**
   * @brief: Matrix form of zonotope's bounds
   */
    template <typename Number>
	Matrix_t<Number> Zono_border_matrix(const Zonotope<Number>& Zono);

    /**
   * @brief: Matrix form of zonotope's tiling (will change generators' order)
   */
    template <typename Number>
	Matrix_t<Number> Zono_tiling_matrix(Zonotope<Number>& Zono);

    /**
   * @brief: Tiling of zonotope
   */
    template <typename Number>
	vector<Zonotope<Number>> Zono_tiling(const Zonotope<Number>& Zono);

    /**
   * @brief: Convert Matrix form of zonotope to G-representation of zonotope
   */
    template <typename Number>
	vector<Zonotope<Number>> Matrix2Zono(const Zonotope<Number>& Zono, const Matrix_t<Number>& Matrix);

    /**
   * @brief: Returns a Zonotope which is projected onto the specified dimensions
   */
    template <typename Number>
	Zonotope<Number> project(Zonotope<Number> const & Z, int dim);

} // namespace reachSolver
#ifdef SEPARATE_TEMPLATE
#include "../../src/zonotope/Zonotope.tpp"
#endif
#endif // ZONOTOPE_H