/**
 * @file   polyZonotope.cpp
 * @brief  polynomial Zonotope class representation
 * @author Changyuan Zhao
 * @date  Nov 2021
 * @version 1.0
 *
 * Reference:
 *   CORA ../contSet/polyZonotope/polyZonotope.m
 */

#include <zonotope/polyZonotope.h>
#include <zonotope/Zonotope.h>
using namespace std;

namespace reachSolver
{

	// template class polyZonotope<double>;
	/*****************************************************************************
 *                                                                           *
 *                           Constructors and Destructors                    *
 *                                                                           *
 *****************************************************************************/
	/**
 * @brief Constructor with no params
 */
	template <typename Number>
	polyZonotope<Number>::polyZonotope()
		: dimension_(0)
		, center_(0)
		, generators_(0, 0)
		, expMat_(0, 0)
		, Grest_(0, 0)
	{}

	/**
 * @brief Constructor with a given dimension
 * @param dimension Dimensionality of Zonotope
 */
	template <typename Number>
	polyZonotope<Number>::polyZonotope(size_t dimension)
		: dimension_(dimension)
		, center_(Vector_t<Number>::Zero(dimension))
		, generators_(Matrix_t<Number>::Zero(dimension, 1))
		, expMat_(Matrix_t<int>::Zero(1, 1))
		, Grest_(Matrix_t<Number>::Zero(dimension, 1))
	{
		assert(dimension != 0);
		id_ = std::vector<size_t>(expMat_.rows());
		std::iota(id_.begin(), id_.end(), 1);
	}

	/**
 * @brief Constructor with center and generators.
 * @param center A  vector
 * @param generators A  matrix
 */
	template <typename Number>
	polyZonotope<Number>::polyZonotope(const Vector_t<Number>& center,
									   const Matrix_t<Number>& generators,
									   const Matrix_t<int>&	   expMat,
									   const Matrix_t<Number>& Grest)
	{
		dimension_	= center.rows();
		center_		= center;
		// generators_ = generators;

		if(expMat.cols() != 0){
			Matrix_t<int>	 ExpNew;
			Matrix_t<Number> Gnew;
			removeRedundantExponents(expMat, generators, ExpNew, Gnew);
			generators_ = Gnew;
			expMat_ = ExpNew;
			// generators_ = generators;
			// expMat_ = expMat;
		}else{
			generators_ = generators;
			expMat_ = expMat;
		}

		if (generators_.size() == 0)
		{
			generators_.resize(dimension_, 0);
		}
		// expMat_ = expMat;
		Grest_	= Grest;
		if (Grest_.size() == 0)
		{
			Grest_.resize(dimension_, 0);
		}
		id_ = std::vector<size_t>(expMat_.rows());
		std::iota(id_.begin(), id_.end(), 1);
		assert(center_.rows() == generators_.rows());
	}

	template <typename Number>
	polyZonotope<Number>::polyZonotope(const Vector_t<Number>&	 center,
									   const Matrix_t<Number>&	 generators,
									   const Matrix_t<int>&		 expMat,
									   const Matrix_t<Number>&	 Grest,
									   const std::vector<size_t> id)
	{
		dimension_	= center.rows();
		center_		= center;
		generators_ = generators;
		if (generators_.size() == 0)
		{
			generators_.resize(dimension_, 0);
		}
		expMat_ = expMat;
		Grest_	= Grest;
		if (Grest_.size() == 0)
		{
			Grest_.resize(dimension_, 0);
		}
		id_ = id;
		assert(center_.rows() == generators_.rows());
		assert(id_.size() == Grest_.cols());
	}

	template <typename Number>
	polyZonotope<Number>::polyZonotope(const Zonotope<Number> zono)
	{
		Vector_t<Number> c = zono.center();
		Matrix_t<Number> G = zono.generators();
		Matrix_t<int>	 expMat;
		expMat.setIdentity(G.cols(), G.cols());
		dimension_	= c.rows();
		center_		= c;
		generators_ = G;
		expMat_		= expMat;
		Grest_.resize(dimension_, 0);
		id_ = std::vector<size_t>(expMat_.rows());
		std::iota(id_.begin(), id_.end(), 1);
	}

	/**
 * @brief: destructor
 */
	template <typename Number>
	polyZonotope<Number>::~polyZonotope()
	{}

	/*****************************************************************************
 *                                                                           *
 *                            Public Functions on Set Properties             *
 *                                                                           *
 *****************************************************************************/

	/**
 * @brief Get the dimension of Zonotope
 * @return the dimension
 */
	template <typename Number>
	size_t polyZonotope<Number>::dimension() const
	{
		return dimension_;
	}

	/**
 * @brief Get the current center
 * @return center a nx1 matrix
 */
	template <typename Number>
	const Vector_t<Number>& polyZonotope<Number>::center() const
	{
		return center_;
	}

	/**
 * @brief Replaces the current center with the parameter center
 * @param center a nx1 matrix
 */
	template <typename Number>
	void polyZonotope<Number>::set_center(const Vector_t<Number>& center)
	{
		if (dimension_ == 0)
		{
			dimension_	= center.rows();
			generators_ = Matrix_t<Number>::Zero(dimension_, 0);
			Grest_		= Matrix_t<Number>::Zero(dimension_, 0);
		}
		assert((std::size_t) center.rows() == dimension_ && "Center has to have same dimensionality as zonotope.");
		center_ = center;
		// uniteEqualVectors();
		// DeleteZeroGenerators();
	}

	/**
 * @brief Get the current generators
 * @return center a nxm matrix
 */
	template <typename Number>
	const Matrix_t<Number>& polyZonotope<Number>::generators() const
	{
		return generators_;
	}

	/**
 * @brief Replaces the current matrix of generators with the parameter
 * generators
 * @param generators a nxm matrix
 */
	template <typename Number>
	void polyZonotope<Number>::set_generators(const Matrix_t<Number>& generators)
	{
		if (dimension_ == 0)
		{
			dimension_ = generators.rows();
			center_	   = Vector_t<Number>::Zero(dimension_);
			Grest_	   = Matrix_t<Number>::Zero(dimension_, 0);
		}
		// assert((std::size_t)generators.cols() == Grest_.cols() && "Generators have
		// to have same colunmn number as Grest" );
		generators_ = generators;
		if (generators_.size() == 0)
		{
			generators_.resize(dimension_, 0);
		}
		assert((std::size_t) generators_.rows() == dimension_
			   && "Generators have to have same dimensionality as zonotope");
		// uniteEqualVectors();
		// DeleteZeroGenerators();
	}

	/**
 * @brief Get the current expMat
 * @return center a nxm matrix
 */
	template <typename Number>
	const Matrix_t<int>& polyZonotope<Number>::expMat() const
	{
		return expMat_;
	}

	/**
 * @brief Replaces the current matrix of gexpMat with the parameter expMat
 * @param generators a nxm matrix
 */
	template <typename Number>
	void polyZonotope<Number>::set_expMat(const Matrix_t<int>& expMat)
	{
		if (dimension_ == 0)
		{
			dimension_	= expMat.rows();
			center_		= Vector_t<Number>::Zero(dimension_);
			generators_ = Matrix_t<Number>::Zero(dimension_, expMat.cols());
			Grest_		= Matrix_t<Number>::Zero(dimension_, 0);
		}
		// assert((std::size_t) expMat.rows() > 0 && "expMat have to have at least one row");
		assert((std::size_t) expMat.cols() == (std::size_t) generators_.cols()
			   && "expMat have to have same colunmn number as Generators");
		expMat_ = expMat;
		id_		= std::vector<size_t>(expMat_.rows());
		std::iota(id_.begin(), id_.end(), 1);
		// uniteEqualVectors();
		// DeleteZeroGenerators();
	}

	template <typename Number>
	const std::vector<size_t>& polyZonotope<Number>::id() const
	{
		return id_;
	}

	template <typename Number>
	void polyZonotope<Number>::set_id(const std::vector<size_t>& id)
	{
		assert(dimension_ != 0);
		id_ = id;
	}

	/**
 * @brief Get the current Grest
 * @return center a nxm matrix
 */
	template <typename Number>
	const Matrix_t<Number>& polyZonotope<Number>::Grest() const
	{
		return Grest_;
	}

	/**
 * @brief Replaces the current matrix of Grest with the parameter Grest
 * @param generators a nxm matrix
 */
	template <typename Number>
	void polyZonotope<Number>::set_Grest(const Matrix_t<Number>& Grest)
	{
		if (dimension_ == 0)
		{
			dimension_	= Grest.rows();
			center_		= Vector_t<Number>::Zero(dimension_);
			generators_ = Matrix_t<Number>::Zero(dimension_, 0);
		}
		Grest_ = Grest;
		if (Grest_.size() == 0)
		{
			Grest_.resize(dimension_, 0);
		}
		assert((std::size_t) Grest_.rows() == dimension_ && "Grest have to have same dimensionality as zonotope");
		// uniteEqualVectors();
		// DeleteZeroGenerators();
	}

	/**
 * @brief Number of generators
 * @return number of generators
 */
	template <typename Number>
	size_t polyZonotope<Number>::numGenerators() const
	{
		return generators_.cols();
	}

	/**
 * @brief: remove specified column Grests
 * @param colomn a given column
 * @return the zonotope without the specified column Grests
 */
	template <typename Number>
	void polyZonotope<Number>::RemoveGrest(unsigned int colomn)
	{
		Eigen::Index numRows = Grest_.rows();
		Eigen::Index numCols = Grest_.cols() - 1;

		if (colomn < numCols)
		{
			Grest_.block(0, colomn, numRows, numCols - colomn) = Grest_.block(0, colomn + 1, numRows, numCols - colomn);
		}

		Grest_.conservativeResize(numRows, numCols);
	}

	template <typename Number>
	void polyZonotope<Number>::Removegenerators(unsigned int colomn)
	{
		Eigen::Index numRows = generators_.rows();
		Eigen::Index numCols = generators_.cols() - 1;

		if (colomn < numCols)
		{
			generators_.block(0, colomn, numRows, numCols - colomn)
				= generators_.block(0, colomn + 1, numRows, numCols - colomn);
		}

		generators_.conservativeResize(numRows, numCols);
	}

	template <typename Number>
	void polyZonotope<Number>::RemoveexpMat(unsigned int colomn)
	{
		Eigen::Index numRows = expMat_.rows();
		Eigen::Index numCols = expMat_.cols() - 1;

		if (colomn < numCols)
		{
			expMat_.block(0, colomn, numRows, numCols - colomn)
				= expMat_.block(0, colomn + 1, numRows, numCols - colomn);
		}

		expMat_.conservativeResize(numRows, numCols);
	}

	/**
 * @brief Removes zero generators in generator matrix
 */
	template <typename Number>
	void polyZonotope<Number>::DeleteZeroGenerators()
	{
		Vector_t<Number> zero_vector = Vector_t<Number>::Zero(dimension_);

		std::vector<unsigned> zeroIndex;
		for (unsigned i = 0; i < generators_.cols(); i++)
		{
			if (generators_.col(i) == zero_vector)
			{
				zeroIndex.push_back(i);
			}
		}

		for (std::vector<unsigned>::reverse_iterator r_it = zeroIndex.rbegin(); r_it != zeroIndex.rend(); ++r_it)
		{
			RemoveGrest(*r_it);
		}
	}

	/**
 * @brief Changes the dimension of a polyZonotope. if new_dim > old dim, new
 * rows are initialized with null
 * @param new_dim The new dimension of the Zonotope
 * @return True, if change in dimension was successful
 */
	template <typename Number>
	bool polyZonotope<Number>::ChangeDimension(size_t new_dim)
	{
		assert(new_dim != 0 && "Cannot change dimensionality of zonotope to zero");
		if (new_dim == dimension_)
		{
			return false;
		}
		else
		{
			center_.conservativeResize(new_dim, Eigen::NoChange);
			generators_.conservativeResize(new_dim, Eigen::NoChange);
			Grest_.conservativeResize(new_dim, Eigen::NoChange);
			expMat_.conservativeResize(new_dim, Eigen::NoChange);

			// If new dim > old dim, initialize all new rows to zero
			for (unsigned i = dimension_; i < new_dim; i++)
			{
				center_.row(i).setZero();
				generators_.row(i).setZero();
				Grest_.row(i).setZero();
				expMat_.row(i).setZero();
			}

			dimension_ = new_dim;
			return true;
		}
	}

	/**
 * @brief: compare two vectors' difference between 1-norm and infty-norm
 * @param v1 vector1
 * @param v2 vector2
 * @return True if v1 less than v2
 */
	/*template <typename Number>
bool compareVectors( const Vector_t<Number>& v1, const Vector_t<Number>& v2 ) {
        Number v1_sum = v1.array().abs().matrix().sum();
        Number v2_sum = v2.array().abs().matrix().sum();
        Number v1_inf = v1.array().abs().maxCoeff();
        Number v2_inf = v2.array().abs().maxCoeff();

        return ( v1_sum - v1_inf ) < ( v2_sum - v2_inf );
}*/
	template <typename Number>
	std::vector<size_t> sort_indexes(const std::vector<Number> v)
	{
		std::vector<size_t> idx(v.size());
		std::iota(idx.begin(), idx.end(), 0);
		std::sort(idx.begin(), idx.end(), [&v](size_t i1, size_t i2) { return v[i1] > v[i2]; });
		return idx;
	}

	template <typename Number>
	void polyZonotope<Number>::Reduce(unsigned int limitOrder)
	{

		clock_t start, end;
		// start = clock();

		// cout << endl << "USE Reduce !!!" << endl << endl;
		// extract dimensions
		size_t N = dimension_;
		size_t P = generators_.cols();
		size_t Q = Grest_.cols();

		// number of generators that stay unreduced (N generators are added again
		// after reduction)
		size_t tempk = std::floor(N * limitOrder - N);
		size_t K	 = std::max((size_t) 0, tempk);

		// cout << P << ' ' << Q << ' ' << N << ' ' << limitOrder << endl;

		// check if itis necessary to reduce the order
		if (P + Q > N * limitOrder && K >= 0)
		{
			// concatenate all generators
			Matrix_t<Number> G(N, generators_.cols() + Grest_.cols());
			G << generators_, Grest_;
			// half the generator length for eponents that are all even

			Matrix_t<int> tempexpMat = expMat_;
			for (auto i = 0; i < tempexpMat.rows(); i++)
			{
				for (auto j = 0; j < tempexpMat.cols(); j++)
				{
					tempexpMat(i, j) = tempexpMat(i, j) % 2;
				}
			}
			Matrix_t<int> temp1;
			temp1.setOnes(expMat_.rows(), expMat_.cols());
			temp1			   = temp1 - tempexpMat;
			Matrix_t<int> temp = temp1.colwise().prod();

			// Matrix_t<int> temp =
			// (Eigen::MatrixXd::Ones(expMat_.rows(),expMat_.cols())-
			// tempexpMat).colwise().prod();
			for (auto i = 0; i < temp.cols(); i++)
			{
				if (temp(0, i) == 1)
				{
					G.col(i) = 0.5 * G.col(i);
				}
			}
			// calculate the length of the generator vectors with a special metric
			// Matrix_t<Number> lens = (generators_.cwiseProduct(generators_)).colwise().sum();
			Matrix_t<Number> lens = (G.cwiseProduct(G)).colwise().sum();

			// cout << lens << endl;

			// determine the smallest generators (= generators that are removed)
			std::vector<Number> sortedGenerators;
			for (auto i = 0; i < lens.cols(); i++)
			{
				sortedGenerators.push_back(lens(0, i));
			}

			// cout << sortedGenerators << endl;

			std::vector<size_t> idx2 = sort_indexes(sortedGenerators);

			// cout << idx2 << endl;
			// cout << K << endl;

			std::vector<size_t> idx(idx2.begin() + K, idx2.end());

			// cout << idx.size() << endl;
			// cout << endl << "HERE GOOD !!!" << endl << endl;

			// split the indizes into ones for dependent and independent
			// generators
			std::vector<size_t> indDep;
			std::vector<size_t> indInd;
			for (auto i = 0; i < idx.size(); i++)
			{
				if (idx[i] <= P - 1)
				{
					indDep.push_back(idx[i]);
				}
				else
				{
					indInd.push_back(idx[i]);
				}
			}
			for (auto & i : indInd)
			{
				i -= P;
			}

			// cout << endl << "HERE GOOD !!!" << endl << endl;
			// end = clock();
			// start = clock();
			//这部分运行比较慢

			// constract a zonotope from the generators that are removed
			Matrix_t<Number> Grem(N, indDep.size());

			// cout << indDep << endl;
			// cout << generators_ << endl;
			// cout << expMat_ << endl;
			for (auto i = 0; i < indDep.size(); i++)
			{
				Grem.col(i) = generators_.col(indDep[i]);
			}
	
			Matrix_t<int> Erem(expMat_.rows(), indDep.size());
			for (auto i = 0; i < indDep.size(); i++)
			{
				// cout << indDep[i] << endl;
				Erem.col(i) = expMat_.col(indDep[i]);
			}

			// cout << Grem << endl;
			// cout << Erem << endl;

			Matrix_t<Number> GrestRem(N, indInd.size());
			for (auto i = 0; i < indInd.size(); i++)
			{
				GrestRem.col(i) = Grest_.col(indInd[i]);
			}

			Vector_t<Number>	 tempc = Eigen::MatrixXd::Zero(N, 1);
			polyZonotope<Number> pZtemp(tempc, Grem, Erem, GrestRem);

			Zonotope<Number> zono = pZtemp.toZonotope();

			// end = clock();
			// start = clock();

			zono.Reduce(1);

			// end = clock();

			set_center(center_ + zono.center());

			// cout << endl << "HERE GOOD !!!" << endl << endl;

			// Matrix_t<Number> Grest(N, Grest_.cols() + zono.generators().cols());
			// Grest << Grest_, zono.generators();

			// set_Grest(Grest);

			// cout << endl << "HERE GOOD !!!" << endl << endl;

			sort(indDep.begin(),indDep.end());
			sort(indInd.begin(),indInd.end());

			// start = clock();

			RemoveColumn_Vec(generators_, indDep);
			RemoveColumn_Vec(expMat_, indDep);

			RemoveColumn_Vec(Grest_, indInd);

			// for (auto i = 0; i < indDep.size(); i++)
			// {
			// 	Removegenerators(indDep[i] - i);
			// 	RemoveexpMat(indDep[i] - i);
			// }
			// for (auto i = 0; i < indInd.size(); i++)
			// {
			// 	RemoveGrest(indInd[i] - i);
			// }
			// cout << indDep.size() << endl;
			// cout << indInd.size() << endl;

			// end = clock();

			// for(auto i = 0; i < zeroColnum; i++){
			// 	RemoveColumn(matrix,zeroCol[i]-i);
			// }

			// cout << Grest_ << endl;
			// cout << indInd << endl;

			Matrix_t<Number> Grest(N, Grest_.cols() + zono.generators().cols());
			Grest << Grest_, zono.generators();

			set_Grest(Grest);

			// cout << Grest << endl;
		}

		// cout << endl << "HERE GOOD !!!" << endl << endl;

		//remove all exponent vector dimensions that have no entries
		// Matrix_t<int> expMat = expMat_;
		// std::vector<size_t> id = id_;

		// Matrix_t<int> row_min = expMat.rowwise().minCoeff();
		// Matrix_t<int> row_max = expMat.rowwise().maxCoeff();
		// vector<unsigned int> zeroRow;
		// unsigned int numRows = expMat.rows();
		// for(auto i = 0; i < numRows; i++){
		// 	if(row_min(i,0) == 0 && row_max(i,0) == 0){
		// 		zeroRow.push_back(i);
		// 	}
		// }
		// unsigned int zeroRownum = zeroRow.size();
		// for(auto i = 0; i < zeroRownum; i++){
		// 	// RemoveRow(expMat,zeroRow[i] - i);
		// 	id.erase(id.begin() + zeroRow[i] - i);
		// }
		Matrix_t<int> zero_vector = Matrix_t<int>::Zero(1, expMat_.cols());

		std::vector<size_t> zeroIndex;
		for (Eigen::Index i = 0; i < expMat_.rows(); i++)
		{
			if (expMat_.row(i) == zero_vector)
			{
				zeroIndex.push_back(i);
			}
		}

		for (std::vector<size_t>::reverse_iterator r_it = zeroIndex.rbegin(); r_it != zeroIndex.rend(); ++r_it){
			id_.erase(id_.begin() + *r_it);
		}
		// RemoveRow_Vec(expMat, zeroRow);
		RemoveRow_Vec(expMat_, zeroIndex);
		// set_expMat(expMat);
		// set_id(id);

		// cout << endl << "GOOD" << endl << endl;
		// end = clock();
		// cout << "Reduce time cost: " << (double) (end - start) / CLOCKS_PER_SEC << endl;
	}

	/**
 * @brief computes an enclosing zonotope of the polynomial zonotope
 * @return zonotope
 */
	template <typename Number>
	Zonotope<Number> polyZonotope<Number>::toZonotope() const
	{
		Zonotope<Number> result;
		if (generators_.size() != 0)
		{
			Matrix_t<int> tempexpMat = expMat_;
			for (auto i = 0; i < tempexpMat.rows(); i++)
			{
				for (auto j = 0; j < tempexpMat.cols(); j++)
				{
					tempexpMat(i, j) = tempexpMat(i, j) % 2;
				}
			}
			Matrix_t<int> temp1;
			temp1.setOnes(expMat_.rows(), expMat_.cols());
			temp1			   = temp1 - tempexpMat;
			Matrix_t<int> temp = temp1.colwise().prod();

			// Matrix_t<Number> temp =
			// (Eigen::MatrixXd::Ones(expMat_.rows(),expMat_.cols())-tempexpMat).colwise().prod();
			Matrix_t<Number> Gquad(dimension_, 0);
			Matrix_t<Number> Gone(dimension_, 0);
			// size_t dimexp = generators_.rows();
			for (auto i = 0, j = 0, k = 0; i < temp.cols(); i++)
			{
				if (temp(0, i) == 1)
				{
					Gquad.conservativeResize(dimension_, j + 1);
					Gquad.col(j) = generators_.col(i);
					j++;
				}
				else
				{
					Gone.conservativeResize(dimension_, k + 1);
					Gone.col(k) = generators_.col(i);
					k++;
				}
			}

			Matrix_t<Number> c = center_ + 0.5 * Gquad.rowwise().sum();
			Matrix_t<Number> G(dimension_, Gone.cols() + Gquad.cols() + Grest_.cols());
			G << Gone, 0.5 * Gquad, Grest_;
			result.set_center(c);
			result.set_generators(G);
		}
		else
		{
			result.set_center(center_);
			result.set_generators(Grest_);
		}
		return result;
	}

	/**
 * @brief judge wether the zonotope is empty
 * @return True if it's empty
 */
	template <typename Number>
	int polyZonotope<Number>::Empty()
	{
		if (this->dimension_ == 0 && this->center_.size() == 0 && this->generators_.size() == 0
			&& this->Grest_.size() == 0 && this->expMat_.size() == 0)
		{
			return 1;
		}
		else
		{
			return 0;
		}
	}

	/**
 * @brief Clears the generators and center of the Zonotope and sets
 * dimensionality to zero
 */
	template <typename Number>
	void polyZonotope<Number>::Clear()
	{
		generators_.resize(0, 0);
		Grest_.resize(0, 0);
		expMat_.resize(0, 0);
		center_.resize(0, 1);
		dimension_ = 0;
	}

	/**
 * @brief display the zonotope
 */
	template <typename Number>
	void polyZonotope<Number>::Display() const
	{
		std::cout << "This Zonotope's dimension is " << dimension_ << std::endl;
		std::cout << "This Zonotope's center is \n" << center_ << std::endl;
		std::cout << "This Zonotope's generators are \n" << generators_ << std::endl;
	}

	/*****************************************************************************
 *                                                                           *
 *                          Basic Set Operations                               *
 *                                                                           *
 *****************************************************************************/

	template <typename Number>
	polyZonotope<Number> polyZonotope<Number>::operator*(Number num) const
	{
		polyZonotope<Number> result;
		result.set_center(num * center_);
		result.set_generators(num * generators_);
		result.set_Grest(num * Grest_);
		result.set_expMat(expMat_);
		return result;
	}

	/**
 * @brief implement the linear maps of a set, i.e., "*" operator
 *        since zonotope could only be times by right
 * @param matrix a matrix
 * @return a  zonotope = matrix * this zonotope
 */
	template <typename Number>
	polyZonotope<Number> polyZonotope<Number>::Times(const Matrix_t<Number>& matrix) const
	{
		assert(matrix.cols() == center_.rows());
		assert(matrix.cols() == center_.rows());

		polyZonotope<Number> result;
		result.set_center(matrix * center_);
		Matrix_t<Number> G = matrix * generators_;
		result.set_generators(G);
		Matrix_t<Number> Grest;
		if (Grest_.size() != 0)
		{
			Grest = matrix * Grest_;
			result.set_Grest(Grest);
		}
		result.set_expMat(expMat_);
		return result;
	}

	template <typename Number>
	polyZonotope<Number> polyZonotope<Number>::operator*(const Matrix_t<Number>& matrix) const
	{
		return Times(matrix);
	}

	template <typename Number>
	polyZonotope<Number> polyZonotope<Number>::Times(const IntervalMatrix<Number> int_matrix) const
	{
		Matrix_t<Number> M_min = inf(int_matrix);
		Matrix_t<Number> M_max = sup(int_matrix);
		// get center of interval matrix
		Matrix_t<Number> T = 0.5 * (M_max + M_min);
		// get symmetric interval matrix
		Matrix_t<Number> S = 0.5 * (M_max - M_min);

		// Matrix_t<Number> tempInf = T - S;
		// Matrix_t<Number> tempSup = T + S;

		// IntervalMatrix<Number> temp = setIntervalMatrix(tempInf, tempSup);
		// M_min = inf(temp);
		// M_max = sup(temp);
		// T = 0.5 * (M_max + M_min);
		// S = 0.5 * (M_max - M_min);

		// cout << int_matrix << endl;


		Matrix_t<Number> Z = Matrix_t<Number>(center_.rows(), center_.cols() + generators_.cols());
		Z << center_, generators_;
		Matrix_t<Number> Zabssum = Z.cwiseAbs().rowwise().sum();

		Zonotope<Number>	   pz = toZonotope();
		IntervalMatrix<Number> I  = pz.interval();
		Matrix_t<Number>	   centerI(I.rows(), I.cols());
		centerI << 0.5 * (inf(I) + sup(I));
		Matrix_t<Number> radI(I.rows(), I.cols());
		radI << 0.5 * (sup(I) - inf(I));
		Matrix_t<Number> s(I.rows(), I.cols());
		s << centerI.cwiseAbs() + radI;
		// compute new zonotope
		polyZonotope<Number> result;
		result.set_center(T * center_);
		if (generators_.size() != 0)
		{
			result.set_generators(T * generators_);
		}
		if (Grest_.size() != 0)
		{

			// cout << S << endl;
			// cout << s << endl;

			Matrix_t<Number> diagrs = (S * s).asDiagonal();
			Matrix_t<Number> mGrest = T * Grest_;
			Matrix_t<Number> rGrest(mGrest.rows(), diagrs.cols() + mGrest.cols());
			rGrest << mGrest, diagrs;

			// cout << mGrest << endl;
			// cout << diagrs << endl;

			result.set_Grest(rGrest);
		}
		else
		{
			Matrix_t<Number> diagrs = (S * s).asDiagonal();
			result.set_Grest(diagrs);
		}
		result.set_expMat(expMat_);
		return result;
	}

	template <typename Number>
	polyZonotope<Number> polyZonotope<Number>::operator*(const IntervalMatrix<Number> int_matrix) const
	{
		return Times(int_matrix);
	}

	/**
 * @brief Get the minkowski addition of two zonotope,i.e., "+" operator
 * @param another_zonotope
 * @return a  zonotope = zonotope1 + zonotope2
 */
	template <typename Number>
	polyZonotope<Number> polyZonotope<Number>::Plus(const Zonotope<Number>& another_zonotope) const
	{
		assert(dimension_ == another_zonotope.dimension());
		polyZonotope<Number> sum(*this);
		sum.set_center(this->center_ + another_zonotope.center());
		Matrix_t<Number> sum_generators;
		sum_generators.resize(dimension_, Grest_.cols() + another_zonotope.generators().cols());
		sum_generators << Grest_, another_zonotope.generators();
		sum.set_Grest(sum_generators);
		// sum.uniteEqualVectors();
		sum.DeleteZeroGenerators();
		// sum.reduce();
		return sum;
	}

	template <typename Number>
	polyZonotope<Number> polyZonotope<Number>::Plus(const polyZonotope<Number>& another_polyzonotope) const
	{
		assert(dimension_ == another_polyzonotope.dimension());
		polyZonotope<Number> sum;
		sum.set_center(center_ + another_polyzonotope.center());
		Matrix_t<Number> sum_Grest;
		sum_Grest.resize(dimension_, Grest_.cols() + another_polyzonotope.Grest().cols());
		sum_Grest << Grest_, another_polyzonotope.Grest();
		sum.set_Grest(sum_Grest);

		Matrix_t<Number> sum_generators;
		sum_generators.resize(dimension_, generators_.cols() + another_polyzonotope.generators().cols());
		sum_generators << generators_, another_polyzonotope.generators();
		sum.set_generators(sum_generators);

		Matrix_t<int> sum_expMat;
		sum_expMat.setZero(expMat_.rows() + another_polyzonotope.expMat().rows(),
						   expMat_.cols() + another_polyzonotope.expMat().cols());

		sum_expMat.block(0, 0, expMat_.rows(), expMat_.cols()) = expMat_;
		sum_expMat.block(
			expMat_.rows(), expMat_.cols(), another_polyzonotope.expMat().rows(), another_polyzonotope.expMat().cols())
			= another_polyzonotope.expMat();
		sum.set_expMat(sum_expMat);

		// sum.uniteEqualVectors();
		sum.DeleteZeroGenerators();
		// sum.reduce();
		return sum;
	}

	template <typename Number>
	polyZonotope<Number> polyZonotope<Number>::operator+(const Zonotope<Number>& another_zonotope) const
	{
		return Plus(another_zonotope);
	}

	template <typename Number>
	polyZonotope<Number> polyZonotope<Number>::operator+(const polyZonotope<Number>& another_polyzonotope) const
	{
		return Plus(another_polyzonotope);
	}

	template <typename Number>
	polyZonotope<Number> polyZonotope<Number>::operator+(const Vector_t<Number>& vector) const
	{
		polyZonotope<Number> newZ(this->center_, this->generators_, this->expMat_, this->Grest_);
		newZ.set_center(this->center_ + vector);
		return newZ;
	}

	template <typename Number>
	polyZonotope<Number> polyZonotope<Number>::operator+(const Number num) const
	{
		polyZonotope<Number> newZ(this->center_, this->generators_, this->expMat_, this->Grest_);
		size_t				 dim	= newZ.dimension();
		Vector_t<Number>	 vector = Vector_t<Number>(dim);
		for (size_t i = 0; i < dim; i++)
		{
			vector[i] = num;
		}
		newZ.set_center(this->center_ + vector);
		return newZ;
	}

	template <typename Number>
	polyZonotope<Number> polyZonotope<Number>::operator-(const Vector_t<Number>& vector) const
	{
		polyZonotope<Number> newZ(this->center_, this->generators_, this->expMat_, this->Grest_);
		newZ.set_center(this->center_ - vector);

		return newZ;
	}

	template <typename Number>
	polyZonotope<Number> polyZonotope<Number>::operator-(const Number num) const
	{
		polyZonotope<Number> newZ(this->center_, this->generators_, this->expMat_, this->Grest_);
		size_t				 dim	= newZ.dimension();
		Vector_t<Number>	 vector = Vector_t<Number>(dim);
		for (size_t i = 0; i < dim; i++)
		{
			vector[i] = num;
		}
		newZ.set_center(this->center_ - vector);
		return newZ;
	}

	template <typename Number>
	polyZonotope<Number> polyZonotope<Number>::operator+(const IntervalMatrix<Number> inval) const
	{
		Zonotope<Number> newZ(inval);
		return Plus(newZ);
	}

	/**
 * @brief Generates a polyzonotope that encloses a polyzonotopes and its linear
 * transformation
 * @param another_zonotope second polyzonotope object, satisfying Z2 = (M * Z1)
 * + Zplus
 * @return zonotope, that encloses Z1 and Z2
 */
	template <typename Number>
	polyZonotope<Number> polyZonotope<Number>::enclose(const polyZonotope<Number>& another_polyzonotope) const
	{
		Matrix_t<Number> generators(dimension_, generators_.cols() * 2 + 1);
		generators << 0.5 * generators_ + 0.5 * another_polyzonotope.generators(),
			0.5 * generators_ - 0.5 * another_polyzonotope.generators(),
			0.5 * center_ - 0.5 * another_polyzonotope.center();

		// cout << generators << endl;

		Vector_t<Number> center = 0.5 * center_ + 0.5 * another_polyzonotope.center();
		Matrix_t<int>	 temp;
		temp.setOnes(1, expMat_.cols());
		Matrix_t<int> expMat(expMat_.rows() + 1, expMat_.cols() + another_polyzonotope.expMat().cols());
		expMat << expMat_, expMat_, 0 * temp, temp;
		Matrix_t<int> expMat2(expMat.rows(), expMat.cols() + 1);
		Matrix_t<int> temp2;
		temp2.setZero(expMat.rows(), 1);
		temp2(expMat.rows() - 1, 0) = 1;
		expMat2 << expMat, temp2;
		Vector_t<Number> temp3 = Eigen::MatrixXd::Zero(dimension_, 1);
		Zonotope<Number> zono1(temp3, Grest_);
		Zonotope<Number> zono2(temp3, another_polyzonotope.Grest());
		zono1					   = zono1.enclose(zono2);
		Matrix_t<Number>	 Grest = zono1.generators();
		polyZonotope<Number> result(center, generators, expMat2, Grest);
		return result;
	}

	template <typename Number>
	polyZonotope<Number> polyZonotope<Number>::cartProd(const polyZonotope<Number> pZ2)
	{
		if (pZ2.dimension() == 0)
		{
			return this;
		}
		if (dimension_ == 0)
		{
			return pZ2;
		}

		// get dimensions
		int dim1 = dimension_;
		int dim2 = pZ2.dimension();

		Matrix_t<Number> G;
		Matrix_t<int>	 expMat;

		Vector_t<Number> c(dim1 + dim2);
		c << center_, pZ2.center();
		if (generators_.size() == 0)
		{
			if (pZ2.generators().size() != 0)
			{
				G.resize(dim1 + dim2, 2 * pZ2.generators().cols());
				G << Eigen::MatrixXd::Zero(dim1, pZ2.generators().cols()), pZ2.generators();
				expMat = pZ2.expMat();
			}
		}
		else
		{
			if (pZ2.generators().size() == 0)
			{
				G.resize(dim1 + dim2, 2 * generators_.cols());
				G << generators_, Eigen::MatrixXd::Zero(dim2, generators_.cols());
				expMat = expMat_;
			}
			else
			{
				G.resize(dim1 + dim2, pZ2.generators().cols() + generators_.cols());
				G << generators_, Eigen::MatrixXd::Zero(dim1, pZ2.generators().cols()),
					Eigen::MatrixXd::Zero(dim2, generators_.cols()), pZ2.generators();
				expMat.resize(expMat_.rows() + pZ2.expMat().rows(), expMat_.cols() + pZ2.expMat().cols());
				Matrix_t<int> temp12;
				temp12.setZero(expMat_.rows(), pZ2.expMat().cols());
				Matrix_t<int> temp21;
				temp21.setZero(pZ2.expMat().rows(), expMat_.cols());

				expMat << expMat_, temp12, temp21, pZ2.expMat();
			}
		}

		// matrix of independent generators
		Matrix_t<Number> Grest;
		if (Grest_.size() == 0)
		{
			if (pZ2.Grest() != 0)
			{
				Grest.resize(dim1 + dim2, 2 * pZ2.Grest().cols());
				Grest << Eigen::MatrixXd::Zero(dim1, pZ2.Grest().cols()), pZ2.Grest();
			}
		}
		else
		{
			if (pZ2.Grest() == 0)
			{
				Grest.resize(dim1 + dim2, 2 * Grest_.cols());
				Grest << Grest_, Eigen::MatrixXd::Zero(dim2, Grest_.cols());
			}
			else
			{
				Grest.resize(dim1 + dim2, pZ2.Grest().cols() + Grest_.cols());
				Grest << Grest_, Eigen::MatrixXd::Zero(dim1, pZ2.Grest().cols()),
					Eigen::MatrixXd::Zero(dim2, Grest_.cols()), pZ2.Grest();
			}
		}
		polyZonotope<Number> result(c, G, Grest, expMat);
		return result;
	}

	template <typename Number>
	polyZonotope<Number> polyZonotope<Number>::cartProd(const Zonotope<Number> pZ2)
	{

		// cout << endl << "HAHAHA" << endl << endl;
		// cout << *this << endl;
		// cout << pZ2 << endl;

		if (pZ2.dimension() == 0)
		{
			return polyZonotope<Number>(*this);
		}

		// cout << *this << endl;
		// cout << pZ2 << endl;

		Matrix_t<Number> rG;
		Matrix_t<int>	 rexpMat;

		polyZonotope<Number> pZ(pZ2.center(), rG, rexpMat, pZ2.generators());
		if (dimension_ == 0)
		{
			return pZ;
		}

		// get dimensions
		int dim1 = dimension_;
		int dim2 = pZ.dimension();

		Matrix_t<Number> G;
		Matrix_t<int>	 expMat;

		Vector_t<Number> c(dim1 + dim2);
		c << center_, pZ.center();
		if (generators_.size() == 0)
		{
			if (pZ.generators().size() != 0)
			{
				G.resize(dim1 + dim2, 2 * pZ.generators().cols());
				G << Eigen::MatrixXd::Zero(dim1, pZ.generators().cols()), pZ.generators();
				expMat = pZ.expMat();
			}
		}
		else
		{
			if (pZ.generators().size() == 0)
			{
				G.resize(dim1 + dim2, generators_.cols());
				G << generators_, Eigen::MatrixXd::Zero(dim2, generators_.cols());
				expMat = expMat_;
			}
			else
			{
				G.resize(dim1 + dim2, pZ.generators().cols() + generators_.cols());
				G << generators_, Eigen::MatrixXd::Zero(dim1, pZ.generators().cols()),
					Eigen::MatrixXd::Zero(dim2, generators_.cols()), pZ.generators();
				expMat.resize(expMat_.rows() + pZ.expMat().rows(), expMat_.cols() + pZ.expMat().cols());
				Matrix_t<int> temp12;
				temp12.setZero(expMat_.rows(), pZ.expMat().cols());
				Matrix_t<int> temp21;
				temp21.setZero(pZ.expMat().rows(), expMat_.cols());

				expMat << expMat_, temp12, temp21, pZ.expMat();
			}
		}

		// cout << endl << "HAHAHA" << endl << endl;

		// matrix of independent generators
		Matrix_t<Number> Grest;
		if (Grest_.size() == 0)
		{
			if (pZ.Grest().size() != 0)
			{
				Grest.resize(dim1 + dim2, pZ.Grest().cols());
				Grest << Eigen::MatrixXd::Zero(dim1, pZ.Grest().cols()), pZ.Grest();
			}
		}
		else
		{
			if (pZ.Grest().size() == 0)
			{
				// cout << endl << "HAHAHA" << endl << endl;
				Grest.resize(dim1 + dim2, Grest_.cols());
				Grest << Grest_, Eigen::MatrixXd::Zero(dim2, Grest_.cols());
				// cout << endl << "HAHAHA" << endl << endl;
			}
			else
			{
				// cout << endl << "HAHAHA" << endl << endl;
				Grest.resize(dim1 + dim2, pZ.Grest().cols() + Grest_.cols());
				Grest << Grest_, Eigen::MatrixXd::Zero(dim1, pZ.Grest().cols()),
					Eigen::MatrixXd::Zero(dim2, Grest_.cols()), pZ.Grest();
				// cout << endl << "HAHAHA" << endl << endl;
			}
		}
		polyZonotope<Number> result(c, G, expMat, Grest);
		return result;
	}

	template <typename Number>
	polyZonotope<Number> polyZonotope<Number>::quadMap(std::vector<Matrix_t<Number>> H)
	{
		// clock_t start,end;
		// start = clock();

		// compute an over-approximation of the quadratic map
		//{x_i = x^T Q{i} x | x \in pZ}
		// split into a zonotope Z that overapproximates the dependent generators and
		// a zonotope Zrem that caontains the independent generators

		polyZonotope<Number> pZtemp(*this);
		Matrix_t<Number>	 Grest2;
		pZtemp.set_Grest(Grest2);
		Zonotope<Number> Z = pZtemp.toZonotope();
		Zonotope<Number> Zrem(0. * center_, Grest_);

		// extend generator and exponent matrix by center
		Matrix_t<Number> Gext(dimension_, center_.cols() + generators_.cols());
		Gext << center_, generators_;
		Matrix_t<int> Eext(expMat_.rows(), 1 + expMat_.cols());
		Matrix_t<int> tempEext;
		tempEext.setZero(expMat_.rows(), 1);
		Eext << tempEext, expMat_;

		// initialize the resulting generator and exponent matrix
		int N	= Gext.cols();
		int dim = H.size();
		int M	= N * (N + 1) / 2;

		Matrix_t<int> Equad;
		Equad.setZero(expMat_.rows(), M);
		Matrix_t<Number> Gquad = Eigen::MatrixXd::Zero(dim, M);

		// create the exponent matrix that corresponds to the quadratic map
		int counter = 0;
		for (auto j = 0; j < N; j++)
		{
			Matrix_t<int> Eexttemp;
			Eexttemp.setOnes(1, N - j);
			Eexttemp									 = Eext.col(j) * Eexttemp;
			Equad.block(0, counter, Equad.rows(), N - j) = Eext.block(0, j, Eext.rows(), N - j) + Eexttemp;
			counter										 = counter + N - j;
		}
		// loop over all dimensions
		for (auto i = 0; i < dim; i++)
		{
			// quadratic evaluation

			// interval Matrix will be wrong
			Matrix_t<Number> quadMat;
			quadMat = Gext.adjoint() * H[i] * Gext;
			// copy the result into the generator matrix
			counter = 0;
			for (auto j = 0; j < N - 1; j++)
			{
				Gquad.block(i, counter, 1, N - j)
					<< quadMat.block(j, j, 1, 1) , quadMat.block(j, j + 1, 1, N - j - 1) + quadMat.block(j + 1, j, N - j - 1, 1).adjoint();
				counter = counter + N - j;
			}
			Gquad.block(i, counter, 1, 1)
					<< quadMat.block(N - 1, N - 1, 1, 1) ;
		}

		// add up all generators that belong to identical exponents
		Matrix_t<int>	 ExpNew;
		Matrix_t<Number> Gnew;
		removeRedundantExponents(Equad, Gquad, ExpNew, Gnew);

		// cout << Equad << endl;
		// cout << Gquad << endl;
		// cout << ExpNew << endl;
		// cout << Gnew << endl;

		// quadratic and mixed multiplication of remaining generators
		Vector_t<Number> center;
		Matrix_t<Number> Grest;
		if (Grest_.size() != 0)
		{
			// quadratic multiplication
			Zonotope<Number> Ztemp1 = Zrem.quadMap(H);
			Zonotope<Number> Ztemp2 = Z.quadMapMixed(Zrem, H);
			Zonotope<Number> Ztemp3 = Zrem.quadMapMixed(Z, H);

			center = Ztemp1.center() + Ztemp2.center() + Ztemp3.center();
			Matrix_t<Number> tempGrest(
				Ztemp1.generators().rows(),
				Ztemp1.generators().cols() + Ztemp2.generators().cols() + Ztemp3.generators().cols());
			tempGrest << Ztemp1.generators(), Ztemp2.generators(), Ztemp3.generators();

			// delete generators of length zero
			Grest = Matrix_t<Number>(tempGrest.rows(), tempGrest.cols());
			int j = 0;
			for (auto i = 0; i < tempGrest.cols(); i++)
			{
				if (tempGrest.col(i).any())
				{
					Grest.col(j) = tempGrest.col(i);
					j++;
				}
			}
			Grest.conservativeResize(tempGrest.rows(), j);
		}
		else
		{
			center = Eigen::MatrixXd::Zero(H.size(), 1);
		}
		Matrix_t<Number> G;
		Matrix_t<int>	 expMat;
		// assemble the properties of the resulting polynomial zonotope

		if (ExpNew.col(0).sum() == 0)
		{
			center = center + Gnew.col(0);
			G	   = Gnew.block(0, 1, Gnew.rows(), Gnew.cols() - 1);
			expMat = ExpNew.block(0, 1, ExpNew.rows(), ExpNew.cols() - 1);
		}
		else
		{
			G	   = Gnew;
			expMat = ExpNew;
		}
		polyZonotope<Number> result(center, G, expMat, Grest);
		// polyZonotope<Number> result;

		// end = clock();
		// cout << "quadMap time cost: " << (double) (end - start) / CLOCKS_PER_SEC << endl;

		return result;
	}

	template <typename Number>
	IntervalMatrix<Number> polyZonotope<Number>::interval()
	{
		return toZonotope().interval();
	}

	template <typename Number>
	polyZonotope<Number> polyZonotope<Number>::exactPlus(polyZonotope<Number> pZ2)
	{
		// bring the exponent matrices to a commmon representation
		Matrix_t<int>		expMat1;
		Matrix_t<int>		expMat2;
		std::vector<size_t> id = mergeExpMatrix<double>(id_, pZ2.id(), expMat_, pZ2.expMat(), expMat1, expMat2);

		// cout << expMat1 << endl;
		// cout << expMat2 << endl;
		// cout << id << endl;

		Matrix_t<int> ExpNew;
		Matrix_t<Number> Gnew;

		Matrix_t<int> ExpMat(expMat1.rows(),expMat2.cols() + expMat1.cols());
		ExpMat << expMat1, expMat2;
		Matrix_t<Number> G(generators_.rows(), generators_.cols() + pZ2.generators().cols()); 
		G << generators_,  pZ2.generators();

		// cout << ExpMat <<endl;

		removeRedundantExponents(ExpMat, G, ExpNew,Gnew);
		
		// cout << ExpMat <<endl;
		// cout << ExpNew <<endl;

		polyZonotope<Number> result(*this);
		result.set_generators(Gnew);
		result.set_expMat(ExpNew);
		result.set_center(center_ + pZ2.center());

		Matrix_t<Number> Grest(Grest_.rows(),Grest_.cols() + pZ2.Grest().cols());
		Grest << Grest_, pZ2.Grest();

		result.set_Grest(Grest);
		result.set_id(id);

		// polyZonotope<Number> result(*this);

		//cout << result << endl;

		return result;
	}

	template <typename Number>
	std::ostream& operator<<(std::ostream& os, const polyZonotope<Number>& zono)
	{
		os.precision(5);
		if (zono.dimension() != 0)
		{
			os << "Dimension： " << zono.dimension() << std::endl;
			os << "Center:" << std::endl;
			for (auto i = 0; i < zono.center().rows(); i++)
			{
				os << zono.center()(i, 0) << " ";
			}
			os << std::endl;
			os << "Genertors:"
			   << "Number " << zono.generators().cols() << std::endl;
			if(zono.generators().cols() != 0){
				os << "[";
				for (auto i = 0; i < zono.center().rows(); i++)
				{
					for (auto j = 0; j < zono.generators().cols() - 1; j++)
					{
						os << zono.generators()(i, j) << ", ";
					}
				os << zono.generators()(i, zono.generators().cols()-1) << ";" << std::endl;
				}
				os << "]";
				os << std::endl;
				os << "expMat:"
				<< "Number " << zono.expMat().cols() << std::endl;
				os << "[";
				for (auto i = 0; i < zono.expMat().rows(); i++)
				{
					for (auto j = 0; j < zono.expMat().cols() - 1; j++)
					{
						os << zono.expMat()(i, j) << ", ";
					}
					os << zono.expMat()(i, zono.expMat().cols() - 1) << ";" << std::endl;
				}
				os << "]";
				os << std::endl;
			}
			os << "Grest:"
			   << "Number " << zono.Grest().cols() << std::endl;
			if(zono.Grest().cols() != 0){
				os << "[";
				for (auto i = 0; i < zono.Grest().rows(); i++)
				{
					for (auto j = 0; j < zono.Grest().cols() - 1; j++)
					{
						os << zono.Grest()(i, j) << ", ";
					}
						os << zono.Grest()(i, zono.Grest().cols()-1) << ";" << std::endl;
				}
				os << "]";
				os << std::endl;
			}
			os << "id:"
			   << "Number " << zono.id().size() << std::endl;
			for (auto i : zono.id())
			{
				os << i << " ";
				os << std::endl;
			}
		}
		return os;
	}

	template <typename Number>
	void deleteZeros(polyZonotope<Number> & polyZono){
		Matrix_t<Number> generators = polyZono.generators();
		Matrix_t<Number> Grest = polyZono.Grest();
		Matrix_t<int> expMat = polyZono.expMat();
		std::vector<size_t> id = polyZono.id();

		//delete zero generators
		deleteZeros(Grest);

		Matrix_t<double> col_min = generators.colwise().minCoeff();
		Matrix_t<double> col_max = generators.colwise().maxCoeff();
		vector<int> zeroCol;
		unsigned int numCols = generators.cols();
		for(auto i = 0; i < numCols; i++){
			if(col_min(0,i) == 0 && col_max(0,i) == 0){
				zeroCol.push_back(i);
			}
		}
		// unsigned int zeroColnum = zeroCol.size();
		// for(auto i = 0; i < zeroColnum; i++){
		// 	RemoveColumn(generators,zeroCol[i]-i);
		// 	RemoveColumn(expMat,zeroCol[i]-i);
		// }

		RemoveColumn_Vec(generators, zeroCol);
		RemoveColumn_Vec(expMat, zeroCol);
		
		//delete zero exponents
		Matrix_t<int> row_min = expMat.rowwise().minCoeff();
		Matrix_t<int> row_max = expMat.rowwise().maxCoeff();
		vector<unsigned int> zeroRow;
		unsigned int numRows = expMat.rows();
		for(auto i = 0; i < numRows; i++){
			if(row_min(i,0) == 0 && row_max(i,0) == 0){
				zeroRow.push_back(i);
			}
		}
		// unsigned int zeroRownum = zeroRow.size();
		// for(auto i = 0; i < zeroRownum; i++){
		// 	// RemoveRow(expMat,zeroRow[i] - i);
		// 	id.erase(id.begin() + zeroRow[i] - i);
		// }
		for (auto r_it = zeroRow.rbegin(); r_it != zeroRow.rend(); ++r_it){
			id.erase(id.begin() + *r_it);
		}
		RemoveRow_Vec(expMat, zeroRow);
		polyZono.set_generators(generators);
		polyZono.set_Grest(Grest);
		polyZono.set_expMat(expMat);
		polyZono.set_id(id);
	}

	template <typename Number>
	Zonotope<Number> zonotope(polyZonotope<Number> & polyZono){
		Vector_t<Number> center = polyZono.center();
		Matrix_t<Number> generators = polyZono.generators();
		Matrix_t<Number> Grest = polyZono.Grest();
		Matrix_t<int> expMat = polyZono.expMat();
		int cols = expMat.cols();
		int rows = expMat.rows();

		if(generators.rows() == 0 || generators.cols() == 0){
			Zonotope<Number> zono(center,generators);
			return zono;
		}else{
			//determine dependent generators with exponents that are all even
			vector<int> temp(cols,1);
			for(int i = 0; i < cols; i++){
				for(int j = 0; j < rows; j++){
					if(expMat(j,i) % 2 == 1){
						temp[i] = 0;
						break;
					}
				}
			}
			Matrix_t<Number> Gquad(generators.rows(),0);
			Matrix_t<Number> Geven(generators.rows(),0);
			for(int i = 0; i < cols; i++){
				if(temp[i] == 1){
					Gquad.conservativeResize(generators.rows(),Gquad.cols()+1);
					Gquad.col(Gquad.cols()-1) = generators.col(i);
				}else{
					Geven.conservativeResize(generators.rows(),Geven.cols()+1);
					Geven.col(Geven.cols()-1) = generators.col(i);
				}
			}

			//compute zonotope parameter
			Vector_t<Number> newCenter = center + 0.5 * Gquad.rowwise().sum(); 
			Matrix_t<Number> newGenerators(generators.rows(),Geven.cols() + Gquad.cols() + Grest.cols());
			newGenerators << Geven, 0.5 * Gquad, Grest;
			Zonotope<Number> zono(newCenter,newGenerators);
			return zono;
		}	
	}

	template <typename Number>
	polyZonotope<Number> restructure(polyZonotope<Number> polyZono, unsigned int limitOrder){
		Vector_t<Number> center = polyZono.center();
		Matrix_t<Number> generators = polyZono.generators();
		Matrix_t<Number> Grest = polyZono.Grest();
		Matrix_t<int> expMat = polyZono.expMat();

		//check if the maximum zonotope order is exceeded
		Number dim_x = center.rows();
		Number o = (polyZono.id().size() + dim_x)/dim_x;

		//max order satisfied
		if(o <= limitOrder){
			Zonotope<Number> zono;

			// //check if additional generators need to be removed
			// Number o_ = (generators.cols() + dim_x) - limitOrder * dim_x;

			// if(o_ > 0){

			// 	//half the generator length for exponents that are all even
			// 	Matrix_t<Number> Gtemp = generators;
			// 	Matrix_t<Number> expMatOnes;
			// 	expMatOnes.setOnes(expMat.rows(), expMat.cols());

			// 	Matrix_t<Number> expMatMod;
			// 	for(int i = 0; i < expMat.rows(); i++){
			// 		for(int j = 0; j < expMat.cols(); i++){
			// 			expMatMod(i,j) = expMat(i,j) % 2;
			// 		}
			// 	}

			// 	Matrix_t<Number> temp = (expMatOnes - expMatMod).colwise().prod();
			// 	vector<Number> ind = findOnes(temp);
			// 	for(int i = 0; i < ind.size(); i++){
			// 		Gtemp.col(ind[i]) = 0.5 * Gtemp.col(ind[i]);
			// 	}

			// 	//determine length of the generators
			// 	Matrix_t<Number> lens = (Gtemp.cwiseProduct(Gtemp)).colwise().sum();
			// 	vector<Number> sortedGenerators;
			// 	for (auto i = 0; i < lens.cols(); i++)
			// 	{
			// 		sortedGenerators.push_back(lens(0, i));
			// 	}
			// 	vector<size_t> idx2 = sort_indexes(sortedGenerators);

			// 	//reduce the smallest generators
			// 	vector<size_t> idx(idx2.begin(), idx2.end() + o_);

			// 	sort(idx.begin(),idx.end());

			// 	Matrix_t<Number> Grem(generators.rows(), idx.size());
			// 	Matrix_t<int> expMatRem(expMat.rows(), idx.size());

			// 	for(int i = 0; i < idx.size(); i++){
			// 		Grem.col(i) = generators.col(idx[i]);
			// 		expMatRem.col(i) = expMat.col(idx[i]);
			// 	}

			// 	for (auto i = 0; i < idx.size(); i++)
			// 	{
			// 		RemoveColumn(generators, idx[i] - i);
			// 		RemoveColumn(expMat, idx[i] - i);
			// 	}

			// 	// reduce the polynomial zonotope that corresponds to the generators that are removed	
			// 	Vector_t<Number>	 tempc = Eigen::MatrixXd::Zero(dim_x, 1);
			// 	polyZonotope<Number> pZtemp(tempc, Grem, expMat, Grest);
			// 	zono = pZtemp.toZonotope();
			// 	zono.Reduce(1);
			// }else{

			//reduce the zonotope that corresponds to the independent generators
			Vector_t<Number>	 tempc = Eigen::MatrixXd::Zero(dim_x, 1);
			zono = Zonotope<Number>(tempc, Grest);

			// cout << zono << endl;

			zono.ReducePCA(1);
			// zono.Reduce(1);
			// cout << zono << endl;
			//construct the restructured polynomial zonotope
			Matrix_t<Number> Gzono = zono.generators();
			Matrix_t<Number> G(generators.rows(),generators.cols() + Gzono.cols());
			G << generators, Gzono;

			Matrix_t<int> expMatNew(expMat.rows() + Gzono.cols(), expMat.cols() + Gzono.cols());
			Matrix_t<int> zeros1;
			Matrix_t<int> zeros2;
			Matrix_t<int> eye;
			zeros1.setZero(expMat.rows(), Gzono.cols());
			zeros2.setZero(Gzono.cols(), expMat.cols());
			eye.setIdentity(Gzono.cols(), Gzono.cols());

			expMatNew << expMat, zeros1,
						 zeros2, eye;

			Matrix_t<Number> GrestNew;
			polyZonotope<Number> res = polyZonotope<Number>(center + zono.center(), G, expMatNew, GrestNew);
			// res.set_center(center + zono.center());
			// res.set_generators(G);
			// res.set_expMat(expMatNew);
			// res.set_Grest(GrestNew);

			return(res);
		}else{
			//number of dependent generators that need to be removed
			int n = polyZono.id().size() + dim_x -dim_x * limitOrder;

			// cout << n << endl;

			//calculate reference zonotope that is added to the generators in order to compare the volumes
			IntervalMatrix<Number> inter = polyZono.interval();
			Vector_t<Number> znonZeros;
			znonZeros.setZero(dim_x, 1);
			Matrix_t<Number> M_min = inf(inter);
			Matrix_t<Number> M_max = sup(inter);
			Matrix_t<Number> rad = 0.5 * (M_max - M_min) / 100;
			Matrix_t<Number> diag = rad.asDiagonal();
			Zonotope<Number> zonoRef = Zonotope<Number>(znonZeros, diag);

			// cout << zonoRef << endl;

			//calculate the volume for all dependent generators
			vector<Number> Vdep(polyZono.id().size(),0);
			vector<vector<int>> indicesDep(polyZono.id().size());

			for(int i = 0; i < Vdep.size(); i++){

				//find all generators that that depend on the current factor
				Matrix_t<int> sortTemp(1,expMat.cols());
				sortTemp.row(0) = expMat.row(i);
				vector<int> ind = findPositive(sortTemp);

				// cout << ind << endl;

				indicesDep[i] = ind;

				Matrix_t<Number> Gtemp(generators.rows(),ind.size());
				Matrix_t<int> expMattemp(expMat.rows(),ind.size());
				for(int i = 0; i < ind.size(); i++){
					Gtemp.col(i) = generators.col(ind[i]);
					expMattemp.col(i) = expMat.col(ind[i]);
				}

				// cout << Gtemp << endl;

				polyZonotope<Number> pZ_ = polyZonotope<Number>(znonZeros, Gtemp, expMattemp, zonoRef.generators());

				// cout << pZ_ << endl;

				Zonotope<Number> zono_ = pZ_.toZonotope();

				// cout << zono_ << endl;

				//calculate volume of the zonotope over-approximation
				IntervalMatrix<Number> zonoInt = zono_.interval();
				Matrix_t<Number> M_minzonoInt = inf(zonoInt);
				Matrix_t<Number> M_maxzonoInt = sup(zonoInt);
				Matrix_t<Number> radzonoIntTwice = M_maxzonoInt - M_minzonoInt;
				Vdep[i] = radzonoIntTwice.prod();
			}

			// find generators with the smallest volume => smallest
       		// over-approximation by removal
			vector<size_t> idx = sort_indexes_ascend(Vdep);

			vector<Number> ind1(n);
			for(int i = 0; i < n; i++){
				ind1[i] = idx[i];
			}

			//determine the indices of all generators that are removed
			vector<vector<int>> indicesDep_(ind1.size());
			for(int i = 0; i < ind1.size(); i++){
				indicesDep_[i] = indicesDep[ind1[i]];
			}
			vector<int> indDep = indicesDep_[0];
			for(int i = 0; i < indicesDep_.size() - 1; i++){
				indDep = Unique(indDep, indicesDep_[i + 1]);
			}

			Matrix_t<Number> Grem(generators.rows(), indDep.size());
			Matrix_t<int> expMatRem(expMat.rows(), indDep.size());
			for(int i = 0; i < indDep.size(); i++){
				Grem.col(i) = generators.col(indDep[i]);
				expMatRem.col(i) = expMat.col(indDep[i]);
				// RemoveColumn(generators,indDep[i]);
				// RemoveColumn(expMat,indDep[i]);
			}

			sort(indDep.begin(), indDep.end());

			RemoveColumn_Vec(generators, indDep);
			RemoveColumn_Vec(expMat, indDep);

			
			//reduce the zonotope that corresponds to the independent generators
			polyZonotope<Number> PZ_ = polyZonotope<Number>(znonZeros, Grem, expMatRem, Grest);
			Zonotope<Number> zono = PZ_.toZonotope();

			zono.ReducePCA(1);

			//construct the restructured polynomial zonotope
			Matrix_t<Number> Gzono = zono.generators();
			Matrix_t<Number> G(generators.rows(),generators.cols() + Gzono.cols());
			G << generators, Gzono;

			Matrix_t<int> expMatNew(expMat.rows() + Gzono.cols(), expMat.cols() + Gzono.cols());
			Matrix_t<int> zeros1;
			Matrix_t<int> zeros2;
			Matrix_t<int> eye;
			zeros1.setZero(expMat.rows(), Gzono.cols());
			zeros2.setZero(Gzono.cols(), expMat.cols());
			eye.setIdentity(Gzono.cols(), Gzono.cols());

			expMatNew << expMat, zeros1,
						 zeros2, eye;

			Matrix_t<Number> GrestNew;
			polyZonotope<Number> res = polyZonotope<Number>(center + zono.center(), G, expMatNew, GrestNew);
			// res.set_center(center + zono.center());
			// res.set_generators(G);
			// res.set_expMat(expMatNew);
			// res.set_Grest(GrestNew);

			return(res);

		}
		return(polyZono);
	}
	template <typename Number>
	Matrix_t<Number> polyshape(Matrix_t<Number>& matrix){
		Matrix_t<Number> res = matrix;
		RemoveColumn(res, res.cols() - 1);
		int i = 1;
		while(i+1 < res.cols()){
			//判断是否是顶点
			Number error = (res(0,i+1) - res(0,i-1)) * (res(1,i) - res(1,i-1)) - (res(0,i) - res(0,i-1)) * (res(1,i+1) - res(1,i-1));
			if(abs(error) <= 1e-12){
				RemoveColumn(res, i);
			}else{
				i++;
			}
		}
		//判断第一个
		Number error = (res(0,1) - res(0,res.cols()-1)) * (res(1,0) - res(1,res.cols()-1)) - (res(0,0) - res(0,res.cols()-1)) * (res(1,1) - res(1,res.cols()-1));
		if(abs(error) <= 1e-12){
			RemoveColumn(res, 0);
		}
		//判断最后一个
		error = (res(0,0) - res(0,res.cols()-2)) * (res(1,res.cols()-1) - res(1,res.cols()-2)) - (res(0,res.cols()-1) - res(0,res.cols()-2)) * (res(1,0) - res(1,res.cols()-2));
		if(abs(error) <= 1e-12){
			RemoveColumn(res, res.cols()-1);
		}
		return res;
	}
	template <typename Number>
	void getPolygon(polyZonotope<Number> const & pZ, Matrix_t<Number>& V, Matrix_t<Number>& poly){
		
		//zonotope over-approximation
		Zonotope<Number> zono = pZ.toZonotope();

		//calculate vertices of zonotope
		V = convert2polygon(zono);

		//transform to 2D polytope
		poly = polyshape(V);

		Matrix_t<Number> Vtemp = V;
		Matrix_t<Number> polytemp = poly;
		V = Vtemp.adjoint();
		poly = polytemp.adjoint();
	}

	template <typename Number>
	void splitLongestGen(polyZonotope<Number> const & pZ, polyZonotope<Number>& pZ1, polyZonotope<Number>& pZ2){

		//determine longest generator
		Matrix_t<Number> len = (pZ.generators().cwiseProduct(pZ.generators())).colwise().sum();
		int maxRow, inds;
		int s = len.maxCoeff(&maxRow, &inds);

		// cout << inds << endl;

		//find factor with the largest exponent
		int index;
		s = pZ.expMat().col(inds).maxCoeff(&maxRow, &index);
		// index = pZ.id()[index];

		// cout << index << endl;
		// cout << endl << "Good" << endl << endl;
		//determine all generators in which the selected dependent factor occurs
		Matrix_t<int> findexpMat(1,pZ.expMat().cols());
		findexpMat.row(0) = pZ.expMat().row(index);
		vector<int> ind = findPositive(findexpMat);

		//parse input arguments
		Matrix_t<int> maxexpMat(1,ind.size());
		for(int i = 0; i < ind.size(); i++){
			maxexpMat(0,i) = pZ.expMat()(index,ind[i]);
		}
		int polyOrd = maxexpMat.maxCoeff();

		// cout << polyOrd << endl;

		//determine coefficients of polynomials that correspond to (0.5 + 0.5x)^p
		vector<Matrix_t<Number>> polyCoeff1(max(2,polyOrd));
		Matrix_t<Number> polyCoeff11(1,2);
		polyCoeff11 << 1, 1;
		Matrix_t<Number> polyCoeff12(1,3);
		polyCoeff12 << 1, 2, 1;

		polyCoeff1[0] = 0.5 * polyCoeff11;
		polyCoeff1[1] = 0.25 * polyCoeff12;

		if(polyOrd >= 3){
			for(int i = 3; i <= polyOrd; i++){
				Matrix_t<Number> coeffs(1,i+1);
				for(int j = 0; j <= i; j++){
					coeffs(0,j) = combo(i,j);
				}
				polyCoeff1[i-1] = pow(0.5,i) * coeffs;
			}
		}

		//determine coefficients of polynomials that correspond to (-0.5 + 0.5x)^p
		vector<Matrix_t<Number>> polyCoeff2(max(2,polyOrd));
		Matrix_t<Number> polyCoeff21(1,2);
		polyCoeff21 << -1, 1;
		Matrix_t<Number> polyCoeff22(1,3);
		polyCoeff22 << 1, -2, 1;

		polyCoeff2[0] = 0.5 * polyCoeff21;
		polyCoeff2[1] = 0.25 * polyCoeff22;

		if(polyOrd >= 3){
			for(int i = 3; i <= polyOrd; i++){
				Matrix_t<Number> coeffs(1,i+1);
				for(int j = 0; j <= i; j++){
					int factor;
					if((i - j) % 2 == 0){
						factor = 1;
					}else{
						factor = -1;
					}
					coeffs(0,j) = factor * combo(i, j);
				}
				polyCoeff2[i-1] = pow(0.5,i) * coeffs;
			}
		}

		// cout << polyCoeff1 << endl;
		// cout << polyCoeff2 << endl;

		//construct the modified generators for the splitted zonotopes
		Matrix_t<int> expMat1 = pZ.expMat();
		Matrix_t<int> expMat2 = pZ.expMat();
		Matrix_t<Number> G1 = pZ.generators();
		Matrix_t<Number> G2 = pZ.generators();
		Vector_t<Number> c1 = pZ.center();
		Vector_t<Number> c2 = pZ.center();

		for(int i = 0; i < polyOrd; i++){
			Matrix_t<int> findExpMat(1,ind.size());
			for(int i = 0; i < ind.size(); i++){
				findExpMat(0,i) = pZ.expMat()(index,ind[i]);
			}
			vector<int> indTemp = findNumber(findExpMat, i+1);

			Matrix_t<Number> g(pZ.generators().rows(), 1);
			Matrix_t<int> e(pZ.expMat().rows(), 1);

			for(int j = 0; j < indTemp.size(); j++){
				g.col(0) = pZ.generators().col(ind[indTemp[j]]);
				e.col(0) = pZ.expMat().col(ind[indTemp[j]]);

				// cout << g << endl;
				// cout << e << endl;

				//first splitted zonotope
				Matrix_t<Number> coef = polyCoeff1[i];

				G1.col(ind[indTemp[j]]) = coef(0,0) * g;
				expMat1(index,ind[indTemp[j]]) = 0;

				Matrix_t<Number> G1temp = G1;
				G1.conservativeResize(G1.rows(), G1.cols() + coef.cols() - 1);				
				G1 << G1temp, g * coef.middleCols(1, coef.cols() - 1);

				// cout << G1 << endl;

				Matrix_t<int> expMatOnes;
				expMatOnes.setOnes(1,coef.cols() - 1);
				Matrix_t<int> E_ = e * expMatOnes;
				// std::vector<int> idx(coef.cols() - 1);
				// std::iota(idx.begin(), idx.end(), 1);
				
				for(int k = 0; k < coef.cols() - 1; k++){
					E_(index,k) = k + 1;
				}

				Matrix_t<int> expMat1temp = expMat1;
				expMat1.conservativeResize(expMat1.rows(), expMat1.cols() + coef.cols() - 1);
				expMat1 << expMat1temp, E_;

				//second splitted zonotope
				coef = polyCoeff2[i];

				G2.col(ind[indTemp[j]]) = coef(0,0) * g;
				expMat2(index,ind[indTemp[j]]) = 0;

				Matrix_t<Number> G2temp = G2;
				G2.conservativeResize(G2.rows(), G2.cols() + coef.cols() - 1);
				G2 << G2temp, g * coef.middleCols(1, coef.cols() - 1);
				expMatOnes.setOnes(1,coef.cols() - 1);
				E_ = e * expMatOnes;
				// std::vector<int> idx(coef.cols() - 1);
				// std::iota(idx.begin(), idx.end(), 1);
				
				for(int k = 0; k < coef.cols() - 1; k++){
					E_(index,k) = k + 1;
				}

				Matrix_t<int> expMat2temp = expMat2;
				expMat2.conservativeResize(expMat2.rows(), expMat2.cols() + coef.cols() - 1);
				expMat2 << expMat2temp, E_;				
			}
			
			// cout << G1 << endl;
			// cout << G2 << endl;
			// cout << expMat1 << endl;
			// cout << expMat2 << endl;

			//remove the finished indizes from the list
			// for(int k = 0; k < indTemp.size(); k++){
			// 	ind.erase(ind.begin() + indTemp[k] - k);
			// }
			for (auto r_it = indTemp.rbegin(); r_it != indTemp.rend(); ++r_it){
				ind.erase(ind.begin() + *r_it);
			}

			// cout << ind << endl;

			if(ind.size() == 0){
				break;
			}
		}

		//over-approximate all selected generators that did not get splitted
		if(ind.size() != 0){
			Matrix_t<int> expMat1Temp = expMat1;
			expMat1.conservativeResize(expMat1.rows() + 1, expMat1.cols());
			Matrix_t<int> expMatZero;
			expMatZero.setZero(1,expMat1.cols());
			expMat1 << expMat1Temp,
					   expMatZero;

			Matrix_t<int> expMat2Temp = expMat2;
			expMat2.conservativeResize(expMat2.rows() + 1, expMat2.cols());
			expMatZero.setZero(1,expMat2.cols());
			expMat2 << expMat2Temp,
					   expMatZero;

			for(int i = 0; i < ind.size(); i++){
				expMat1(expMat1.rows() - 1, ind[i]) = expMat1(index,ind[i]);
				expMat1(index,ind[i]) = 0;

				expMat2(expMat2.rows() - 1, ind[i]) = expMat2(index,ind[i]);
				expMat2(index,ind[i]) = 0;				
			}	   
		}

		//add every generator with all-zero exponent matrix to the zonotope center
		Matrix_t<int> temp = expMat1.colwise().sum();
		ind = findZero(temp);

		Matrix_t<Number> G1temp(G1.rows(),ind.size());
		for(int i = 0; i < ind.size(); i++){
			G1temp.col(i) = G1.col(ind[i]);
		}
		c1 = c1 + G1temp.rowwise().sum();

		// cout << c1 << endl;

		// for(int i = 0; i < ind.size(); i++){
		// 	RemoveColumn(G1, ind[i] - i);
		// 	RemoveColumn(expMat1, ind[i] - i);
		// }
		RemoveColumn_Vec(G1, ind);
		RemoveColumn_Vec(expMat1, ind);

		temp = expMat2.colwise().sum();
		ind = findZero(temp);

		Matrix_t<Number> G2temp(G2.rows(),ind.size());
		for(int i = 0; i < ind.size(); i++){
			G2temp.col(i) = G2.col(ind[i]);
		}
		c2 = c2 + G2temp.rowwise().sum();
		// for(int i = 0; i < ind.size(); i++){
		// 	RemoveColumn(G2, ind[i] - i);
		// 	RemoveColumn(expMat2, ind[i] - i);
		// }
		RemoveColumn_Vec(G2, ind);
		RemoveColumn_Vec(expMat2, ind);

		//remove redundant exponents
		Matrix_t<Number> G1New;
		Matrix_t<int> expMat1New;
		Matrix_t<Number> G2New;
		Matrix_t<int> expMat2New;	

		removeRedundantExponents(expMat1,G1,expMat1New,G1New);
		removeRedundantExponents(expMat2,G2,expMat2New,G2New);

		//construct the resulting polynomial zonotopes
		Matrix_t<Number> Grest = pZ.Grest();
		pZ1 = polyZonotope<Number>(c1,G1New,expMat1New,Grest);
		pZ2 = polyZonotope<Number>(c2,G2New,expMat2New,Grest);		
	}

	template <typename Number>
	bool isFullDim(polyZonotope<Number> const & pZ){

		//get dimension and rank
		int n = pZ.center().rows();

		Matrix_t<Number> Gnew(pZ.generators().rows(), pZ.generators().cols() + pZ.Grest().cols());
		Gnew << pZ.generators(), pZ.Grest();
		Eigen::ColPivHouseholderQR<Matrix_t<Number>> QR_decompGnew(Gnew);
		int rk = QR_decompGnew.rank();

		//compare dimension and rank
		return n == rk;
	}

	template <typename Number>
	polyZonotope<Number> project(polyZonotope<Number> const & pZ, int dim){

		dim--;
		
		polyZonotope<Number> res;

		res.set_center(pZ.center().row(dim));

		if(pZ.generators().cols() != 0){
			res.set_generators(pZ.generators().row(dim));
		}

		if(pZ.Grest().cols() != 0){
			res.set_Grest(pZ.Grest().row(dim));
		}

		res.set_expMat(pZ.expMat());
		res.set_id(pZ.id());

		return res;
	}
} // namespace reachSolver
