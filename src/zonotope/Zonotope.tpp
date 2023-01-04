/**
 * @file   Zonotope.cpp
 * @brief  Zonotope implementation
 * @author Yunjun
 * @date April 2021
 * @version 1.0
 *
 * Reference:
 *   CORA ../contSet/zonotope/all.m
 */

/**
 * @author Changyuan Zhao
 * @date  Nov 2021
 * @version 1.1
 */

#include <zonotope/Zonotope.h>
#include <zonotope/globalfunctions.h>

using namespace std;
namespace reachSolver
{

	// template class Zonotope<double>;
	/*****************************************************************************
 *                                                                           *
 *                           Constructors and Destructors                    *
 *                                                                           *
 *****************************************************************************/
	/**
 * @brief Constructor with no params
 */
	template <typename Number>
	Zonotope<Number>::Zonotope()
		: dimension_(0)
		, center_(0)
		, generators_(0, 0)
	{}

	/**
 * @brief Constructor with a given dimension
 * @param dimension Dimensionality of Zonotope
 */
	template <typename Number>
	Zonotope<Number>::Zonotope(size_t dimension)
		: dimension_(dimension)
		, center_(Vector_t<Number>::Zero(dimension))
		, generators_(Matrix_t<Number>::Zero(dimension, 1))
	{
		assert(dimension != 0);
	}

	/**
 * @brief: constructor with a given interval matrix
 * @param interval_m an interval matrix
 */
	template <typename Number>
	Zonotope<Number>::Zonotope(IntervalMatrix<Number> interval_m)
	{
		center_		= 0.5 * (inf(interval_m) + sup(interval_m));
		generators_ = (0.5 * (sup(interval_m) - inf(interval_m))).asDiagonal();
		dimension_	= center_.rows();
	}

	/**
 * @brief Constructor with center and generators.
 * @param center A  vector
 * @param generators A  matrix
 */
	template <typename Number>
	Zonotope<Number>::Zonotope(const Vector_t<Number>& center, const Matrix_t<Number>& generators)
		: dimension_(center.rows())
		, center_(center)
		, generators_(generators)
	{
		assert(center.rows() == generators.rows());
	}

	template <typename Number>
	Zonotope<Number>::Zonotope(const Vector_t<Number>& center)
	{
		dimension_ = center.rows();
		center_	   = center;
		generators_.resize(dimension_, 0);
	}

	template <typename Number>

	/**
 * @brief: destructor
 */
	Zonotope<Number>::~Zonotope()
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
	size_t Zonotope<Number>::dimension() const
	{
		return dimension_;
	}

	/**
 * @brief Get the current center
 * @return center a nx1 matrix
 */
	template <typename Number>
	const Vector_t<Number>& Zonotope<Number>::center() const
	{
		return center_;
	}

	/**
 * @brief Replaces the current center with the parameter center
 * @param center a nx1 matrix
 */
	template <typename Number>
	void Zonotope<Number>::set_center(const Vector_t<Number>& center)
	{
		if (dimension_ == 0)
		{
			dimension_	= center.rows();
			generators_ = Matrix_t<Number>::Zero(dimension_, 1);
		}
		assert((std::size_t) center.rows() == dimension_ && "Center has to have same dimensionality as zonotope.");
		center_ = center;
		// uniteEqualVectors();
		DeleteZeroGenerators();
	}

	/**
 * @brief Get the current generators
 * @return center a nxm matrix
 */
	template <typename Number>
	const Matrix_t<Number>& Zonotope<Number>::generators() const
	{
		return generators_;
	}

	/**
 * @brief Replaces the current matrix of generators with the parameter
 * generators
 * @param generators a nxm matrix
 */
	template <typename Number>
	void Zonotope<Number>::set_generators(const Matrix_t<Number>& generators)
	{
		if (dimension_ == 0)
		{
			dimension_ = generators.rows();
			center_	   = Vector_t<Number>::Zero(dimension_);
		}
		assert((std::size_t) generators.rows() == dimension_
			   && "Generators have to have same dimensionality as zonotope");
		generators_ = generators;
		// uniteEqualVectors();
		DeleteZeroGenerators();
	}

	/**
 * @brief Add generators to Zonotope. Simply performs setGenerators if
 * generators was previously not initialized.
 * @param generators a nxm matrix
 * @return true if able to add generators
 */
	template <typename Number>
	bool Zonotope<Number>::AddGenerators(const Matrix_t<Number>& generators)
	{
		if (dimension_ == 0)
		{
			dimension_	= generators.rows();
			center_		= Vector_t<Number>::Zero(dimension_);
			generators_ = generators;
			return true;
		}
		assert((std::size_t) generators.rows() == dimension_
			   && "Added generators have to have same dimensionality as zonotope");
		Matrix_t<Number> tmp = generators_;
		generators_.resize(tmp.rows(), generators.cols() + tmp.cols());
		generators_ << tmp, generators;

		// uniteEqualVectors();
		DeleteZeroGenerators();
		return true;
	}

	/**
 * @brief Get the order
 * @return zonotope order
 */
	template <typename Number>
	Number Zonotope<Number>::order() const
	{
		// object empty.
		if (generators_.rows() == 0)
		{
			return Number(0);
		}
		return Number(generators_.cols()) / Number(generators_.rows());
	}

	/**
 * @brief Number of generators
 * @return number of generators
 */
	template <typename Number>
	size_t Zonotope<Number>::numGenerators() const
	{
		return generators_.cols();
	}

	/**
 * @brief: remove specified column generators
 * @param colomn a given column
 * @return the zonotope without the specified column generator
 */
	template <typename Number>
	void Zonotope<Number>::RemoveGenerator(unsigned int colomn)
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

	/**
 * @brief Removes zero generators in generator matrix
 */
	template <typename Number>
	void Zonotope<Number>::DeleteZeroGenerators()
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
			RemoveGenerator(*r_it);
		}

		// RemoveColumn_Vec(generators_, zeroIndex);
	}

	/**
 * @brief Changes the dimension of a Zonotope. if new_dim > old dim, new rows
 * are initialized with null
 * @param new_dim The new dimension of the Zonotope
 * @return True, if change in dimension was successful
 */
	template <typename Number>
	bool Zonotope<Number>::ChangeDimension(size_t new_dim)
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

			// If new dim > old dim, initialize all new rows to zero
			for (unsigned i = dimension_; i < new_dim; i++)
			{
				center_.row(i).setZero();
				generators_.row(i).setZero();
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
	template <typename Number>
	bool compareVectors(const Vector_t<Number>& v1, const Vector_t<Number>& v2)
	{
		Number v1_sum = v1.array().abs().matrix().sum();
		Number v2_sum = v2.array().abs().matrix().sum();
		Number v1_inf = v1.array().abs().maxCoeff();
		Number v2_inf = v2.array().abs().maxCoeff();

		return (v1_sum - v1_inf) < (v2_sum - v2_inf);
	}

	template <typename Number>
	struct comVec{
		Vector_t<Number> v;
		Number v_sum;
		Number v_inf;
		comVec(Vector_t<Number> & vec)
		: v(vec){
			v_sum = v.array().abs().matrix().sum();
			v_inf = v.array().abs().maxCoeff();
		}

		bool operator<(const comVec& other){
			return (v_sum - v_inf) < (other.v_sum - other.v_inf);
		}
		
	};
	
// 	/**
//  * @brief Reduces the order of a zonotope
//  * @param limitOrder order of reduced zonotope
//  */
// 	template <typename Number>
// 	void Zonotope<Number>::Reduce(unsigned limitOrder)
// 	{
// 		// std::cout << __func__ << ": Current order: " << this->order() << std::endl;
// 		while (this->order() > limitOrder)
// 		{
// 			Matrix_t<Number> generators = generators_;
// 			Eigen::Index	 dim		= generators_.rows();

// 			// create duplicates of generators to sort them
// 			// std::vector<reachSolver::Vector_t<Number>> sortedGenerators1;
// 			std::vector<comVec<Number>> sortedGenerators;
// 			for (Eigen::Index i = 0; i < generators.cols(); i++)
// 			{
// 				// sortedGenerators1.push_back(generators.col(i));
// 				Vector_t<Number> vecSort = generators.col(i);
// 				comVec<Number> comV(vecSort); 
// 				sortedGenerators.push_back(comV);

// 			}

// 			// Sort generators according to the difference between 1-norm and infty-norm
// 			// (i.e.: How suitable are the generators to
// 			// be overapproximated by an interval hull)

// 			// clock_t start, end;

// 			// start = clock();
// 			// std::sort(sortedGenerators1.begin(), sortedGenerators1.end(), compareVectors<Number>);
// 			std::sort(sortedGenerators.begin(), sortedGenerators.end());
// 			// end = clock();
// 			// cout << "sort time cost: " << (double) (end - start) / CLOCKS_PER_SEC << endl;

// 			Eigen::Index	 numRemainingGenerators = Eigen::Index(dim * (limitOrder - 1));
// 			Eigen::Index	 numReducedGenerators	= Eigen::Index(sortedGenerators.size() - numRemainingGenerators);
// 			Matrix_t<Number> remainingGenerators(generators.rows(), numRemainingGenerators);
// 			// Row-wise sum of all chosen generators (absolute value)
// 			Vector_t<Number> sumVector = Vector_t<Number>::Zero(dim);
// 			for (unsigned i = 0; i < numReducedGenerators; i++)
// 			{
// 				sumVector += sortedGenerators[i].v.array().abs().matrix();
// 			}

// 			// inserts the original remaining vectors
// 			for (Eigen::Index i = 0; i < numRemainingGenerators; i++)
// 			{
// 				remainingGenerators.col(i) = sortedGenerators[i + numReducedGenerators].v;
// 			}

// 			// calculate interval hull of reduced generators
// 			Matrix_t<Number> intervalHull = sumVector.asDiagonal();

// 			Matrix_t<Number> reducedGenerators = Matrix_t<Number>(dim, remainingGenerators.cols() + dim);

// 			reducedGenerators << intervalHull, remainingGenerators;
// 			generators_ = reducedGenerators;
// 		}
// 		// std::cout << __func__ << ": Reduced order: " << this->order() << std::endl;
// 	}

		/**
 * @brief Reduces the order of a zonotope (use Girard)
 * @param limitOrder order of reduced zonotope
 */
	template <typename Number>
	void Zonotope<Number>::Reduce(unsigned limitOrder){
		//initialize Z_red
		Zonotope<Number> Zred = *this;

		//pick generators to reduce
		Vector_t<Number> c;
		Matrix_t<Number> Gunred;
		Matrix_t<Number> Gred;
		pickedGenerators(Zred, limitOrder, c, Gunred, Gred);
		Matrix_t<Number> G(Zred.generators().rows(), Gunred.cols());

		//box remaining generators
		if(Gred.cols() != 0){
			Vector_t<Number> d = Gred.cwiseAbs().rowwise().sum();
			Matrix_t<Number> Gbox = d.asDiagonal();
			G.resize(Gred.rows(), Gbox.cols() + Gunred.cols());
			G << Gunred, Gbox;
		}else{
			G << Gunred;
		}
		set_center(c);
		set_generators(G);
	}

	/**
 * @brief Reduces the order of a zonotope using PCA
 * @param limitOrder order of reduced zonotope
 */
	template <typename Number>
	void Zonotope<Number>::ReducePCA(unsigned limitOrder){
		
		//initialize Z_red
		Zonotope<Number> Zred = *this;

		// cout << Zred << endl;

		Vector_t<Number> c;
		Matrix_t<Number> Gunred;
		Matrix_t<Number> Gred;
		pickedGenerators(Zred, limitOrder, c, Gunred, Gred);

		// cout << c << endl;
		// cout << Gunred << endl;
		// cout << Gred << endl;

		if(Gred.cols() != 0){
			//obtain matrix of points from generator matrix
			Matrix_t<Number> V(Gred.rows(),2 * Gred.cols());
			V << Gred, -Gred;

			//compute the arithmetic mean of the vertices
			Matrix_t<Number> mean = V.rowwise().sum()/V.cols();

			// cout << Gred << endl;
			// cout << V << endl;
			// cout << mean << endl;

			//obtain sampling matrix
			Matrix_t<Number> Vones;
			Vones.setOnes(1,V.cols());
			Matrix_t<Number> translation = mean * Vones;
			Matrix_t<Number> sampleMatrix = V - translation; 

			//compute the covariance matrix
			Matrix_t<Number> sampleMatrix_ad = sampleMatrix.adjoint();
			Matrix_t<Number> C = cov(sampleMatrix_ad);

			// cout << C << endl;

			//singular value decomposition
			Eigen::JacobiSVD<Matrix_t<Number>> svd(C, Eigen::ComputeThinU | Eigen::ComputeThinV);  
    		Matrix_t<Number> U = svd.matrixU(); 

			// cout << U << endl;

			//map generators
			Matrix_t<Number> Gtrans = U.adjoint() * Gred;

			//box generators
			Matrix_t<Number> Gtranstemp = Gtrans.cwiseAbs().rowwise().sum();

			// cout << Gtranstemp << endl;

			Matrix_t<Number> Gbox = Gtranstemp.asDiagonal();

			// cout << Gbox << endl;

			//transform generators back
			Gred = U * Gbox;
		}
		//build reduced zonotope
		Matrix_t<Number> G(Gred.rows(), Gred.cols() + Gunred.cols());
		G << Gred, Gunred;
		set_center(c);
		set_generators(G);
	}

	/**
 * @brief judge wether the zonotope is empty
 * @return True if it's empty
 */
	template <typename Number>
	int Zonotope<Number>::Empty()
	{
		if (this->dimension_ == 0 && this->center_ == Vector_t<Number>(0, 1)
			&& this->generators_ == Matrix_t<Number>(0, 0))
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
	void Zonotope<Number>::Clear()
	{
		generators_.resize(0, 0);
		center_.resize(0, 1);
		dimension_ = 0;
	}

	/**
 * @brief display the zonotope
 */
	template <typename Number>
	void Zonotope<Number>::Display() const
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
	Zonotope<Number> Zonotope<Number>::operator*(Number num) const
	{
		Zonotope<Number> result;
		result.set_center(num * center_);
		result.set_generators(num * generators_);
		return result;
	}

	/**
 * @brief implement the linear maps of a set, i.e., "*" operator
 *        since zonotope could only be times by right
 * @param matrix a matrix
 * @return a  zonotope = matrix * this zonotope
 */
	template <typename Number>
	Zonotope<Number> Zonotope<Number>::Times(const Matrix_t<Number>& matrix) const
	{
		assert(matrix.cols() == center_.rows());
		assert(matrix.cols() == center_.rows());
		Zonotope<Number> result;
		result.set_center(matrix * center_);
		result.set_generators(matrix * generators_);
		return result;
	}

	template <typename Number>
	Zonotope<Number> Zonotope<Number>::operator*(const Matrix_t<Number>& matrix) const
	{
		return Times(matrix);
	}

	template <typename Number>
	Zonotope<Number> Zonotope<Number>::Times(const IntervalMatrix<Number> int_matrix) const
	{
		Matrix_t<Number> M_min = inf(int_matrix);
		Matrix_t<Number> M_max = sup(int_matrix);
		// get center of interval matrix
		Matrix_t<Number> T = 0.5 * (M_max + M_min);
		// get symmetric interval matrix
		Matrix_t<Number> S = 0.5 * (M_max - M_min);
		Matrix_t<Number> Z = Matrix_t<Number>(center_.rows(), center_.cols() + generators_.cols());
		Z << center_, generators_;
		Matrix_t<Number> Zabssum = Z.cwiseAbs().rowwise().sum();
		// compute new zonotope
		Zonotope<Number> result;
		result.set_center(T * center_);
		Matrix_t<Number> W = Matrix_t<Number>(S.rows(), generators_.cols() + S.rows());
		W << T * generators_, (Matrix_t<Number>) (S * Zabssum).asDiagonal();
		result.set_generators(W);
		return result;
	}

	template <typename Number>
	Zonotope<Number> Zonotope<Number>::operator*(const IntervalMatrix<Number> int_matrix) const
	{
		return Times(int_matrix);
	}

	/**
 * @brief Get the minkowski addition of two zonotope,i.e., "+" operator
 * @param another_zonotope
 * @return a  zonotope = zonotope1 + zonotope2
 */
	template <typename Number>
	Zonotope<Number> Zonotope<Number>::Plus(const Zonotope& another_zonotope) const
	{
		assert(dimension_ == another_zonotope.dimension());
		Zonotope<Number> sum;
		sum.set_center(this->center_ + another_zonotope.center());
		Matrix_t<Number> sum_generators;
		sum_generators.resize(dimension_, generators_.cols() + another_zonotope.generators().cols());
		sum_generators << generators_, another_zonotope.generators();
		sum.set_generators(sum_generators);
		// sum.uniteEqualVectors();
		sum.DeleteZeroGenerators();
		// sum.reduce();
		return sum;
	}

	template <typename Number>
	Zonotope<Number> Zonotope<Number>::operator+(const Zonotope& another_zonotope) const
	{
		return Plus(another_zonotope);
	}

	template <typename Number>
	Zonotope<Number> Zonotope<Number>::operator+(const Vector_t<Number>& vector) const
	{
		Zonotope<Number> newZ = Zonotope(this->center_, this->generators_);
		newZ.set_center(this->center_ + vector);
		return newZ;
	}

	template <typename Number>
	Zonotope<Number> Zonotope<Number>::operator+(const Number num) const
	{
		Zonotope<Number> newZ	= Zonotope(this->center_, this->generators_);
		size_t			 dim	= newZ.dimension();
		Vector_t<Number> vector = Vector_t<Number>(dim);
		for (size_t i = 0; i < dim; i++)
		{
			vector[i] = num;
		}
		newZ.set_center(this->center_ + vector);
		return newZ;
	}

	template <typename Number>
	Zonotope<Number> Zonotope<Number>::operator-(const Vector_t<Number>& vector) const
	{
		Zonotope<Number> newZ = Zonotope(this->center_, this->generators_);
		newZ.set_center(this->center_ - vector);

		return newZ;
	}

	template <typename Number>
	Zonotope<Number> Zonotope<Number>::operator-(const Number num) const
	{
		// compute center
		size_t			 dim	= this->dimension_;
		Vector_t<Number> vector = Vector_t<Number>(dim);
		for (size_t i = 0; i < dim; i++)
		{
			vector[i] = num;
		}
		return (*this) - vector;
	}

	/**
 * @brief Generates a zonotope that encloses a zonotopes and its linear
 * transformation
 * @param another_zonotope second zonotope object, satisfying Z2 = (M * Z1) +
 * Zplus
 * @return zonotope, that encloses Z1 and Z2
 */
	template <typename Number>
	Zonotope<Number> Zonotope<Number>::enclose(const Zonotope& another_zonotope) const
	{
		Matrix_t<Number> Z1 = Matrix_t<Number>(center_.rows(), center_.cols() + generators_.cols());
		Z1 << center_, generators_;
		Matrix_t<Number> Z2 = Matrix_t<Number>(another_zonotope.center().rows(),
											   another_zonotope.center().cols() + another_zonotope.generators().cols());
		Z2 << another_zonotope.center(), another_zonotope.generators();

		size_t			 generators1 = this->generators().cols() + 1;
		size_t			 generators2 = another_zonotope.generators().cols() + 1;
		Matrix_t<Number> Zcut, Zadd, Zequal;

		if (generators2 <= generators1)
		{
			Zcut   = Z1.middleCols(0, generators2);
			Zadd   = Z1.middleCols(generators2, generators1 - generators2);
			Zequal = Z2;
		}
		else
		{
			Zcut   = Z2.middleCols(0, generators1);
			Zadd   = Z2.middleCols(generators1, generators2 - generators1);
			Zequal = Z1;
		}
		Matrix_t<Number> newZ = Matrix_t<Number>(Zadd.rows(), Zadd.cols() + Zcut.cols() + Zequal.cols());
		newZ << (Zcut + Zequal) / 2, (Zcut - Zequal) / 2, Zadd;
		Zonotope<Number> newZonotope = Zonotope(newZ.middleCols(0, 1), newZ.middleCols(1, newZ.cols() - 1));
		return newZonotope;
	}

	template <typename Number>
	IntervalMatrix<Number> Zonotope<Number>::interval() const
	{
		Matrix_t<Number> Z = Matrix_t<Number>(center_.rows(), center_.cols() + generators_.cols());
		Z << center_, generators_;
		// determine left and right limit
		Matrix_t<Number> delta		= Z.cwiseAbs().rowwise().sum() - center_.cwiseAbs();
		Matrix_t<Number> leftLimit	= center_ - delta;
		Matrix_t<Number> rightLimit = center_ + delta;

		// instantiate interval
		IntervalMatrix<Number> I = setIntervalMatrix(leftLimit, rightLimit);

		return I;
	}
	template <typename Number>
	Zonotope<Number> Zonotope<Number>::cartProd(const Zonotope<Number> Z2)
	{
		Vector_t<Number> c(dimension_ + Z2.dimension());
		c << center_, Z2.center();
		Matrix_t<Number> G(dimension_ + Z2.dimension(), generators_.cols() + Z2.generators().cols());
		G << generators_, Eigen::MatrixXd::Zero(dimension_, Z2.generators().cols()),
			Eigen::MatrixXd::Zero(Z2.dimension(), generators_.cols()), Z2.generators();

		Zonotope<Number> result(c, G);
		return result;
	}

	template <typename Number>
	Matrix_t<Number> Zonotope<Number>::Z()
	{
		Matrix_t<Number> result(dimension_, center_.cols() + generators_.cols());
		result << center_, generators_;
		return result;
	}

	template <typename Number>
	Zonotope<Number> Zonotope<Number>::quadMap(std::vector<Matrix_t<Number>> Q)
	{
		// compute an over-approximation of the quadratic map
		// std::cout<<"there is quadMap"<<std::endl;
		// get matrix of zonotope
		Matrix_t<Number> Zmat = Z();
		// std::cout<<"Zmat"<<std::endl;
		// std::cout<<Zmat<<std::endl;
		int dimQ = Q.size();
		int gens = generators_.cols();

		// init solution
		Vector_t<Number> c = Eigen::MatrixXd::Zero(dimQ, 1);
		Matrix_t<Number> G = Eigen::MatrixXd::Zero(dimQ, 0.5 * (gens * gens + gens) + gens);

		//cout << Zmat;
		// count empty matrices
		Vector_t<int> Qnonempty;
		Qnonempty.setZero(dimQ, 1);

		// for each dimension compute generator elements
		for (auto i = 0; i < dimQ; i++)
		{
			Qnonempty(i) = Q[i].any();
			if (Qnonempty(i))
			{
				// pure quadratic evaluation
				// std::cout<<"where?1"<<std::endl;
				Matrix_t<Number> quadMat = Zmat.adjoint() * Q[i] * Zmat;
				// std::cout<<"where?2"<<std::endl;

				// std::cout<<"quadMat"<<std::endl;
				//std::cout<<quadMat<<std::endl;
				// faster method
				// diag elemengs
				////std::cout<<"0.5*quadMat.block(1,1,gens,gens).diagonal();"<<std::endl;
				////std::cout<<0.5*quadMat.block(1,1,gens,gens).diagonal()<<std::endl;
				G.block(i, 0, 1, gens) = 0.5 * quadMat.block(1, 1, gens, gens).diagonal().adjoint();
				//G.block(i, 0, 1, gens) = 0.5 * quadMat.block(1, 1, gens, gens).diagonal();
				
				//cout << G;

				// center
				// std::cout<<"where?3"<<std::endl;
				c(i, 0) = quadMat(0, 0) + G.block(i, 0, 1, gens).sum();
				// std::cout<<"where?"<<std::endl;
				// off-diagonal elements added, pick via logical indexing
				Matrix_t<Number> quadMatoffdiag = quadMat + quadMat.adjoint();
				//cout << quadMatoffdiag << endl;
				//quadMatoffdiag.conservativeResize(quadMatoffdiag.rows() * quadMatoffdiag.cols(), 1);
				quadMatoffdiag.resize(quadMatoffdiag.rows() * quadMatoffdiag.cols(), 1);
				// std::cout<<"where5?"<<std::endl;
				Matrix_t<bool> kInd;
				kInd.setOnes(gens + 1, gens + 1);
				//kInd.setOnes(gens, gens);
				kInd = kInd.triangularView<Eigen::StrictlyLower>();
				//cout << kInd << endl;
				// std::cout<<"where7?"<<std::endl;
				//kInd.conservativeResize(kInd.rows() * kInd.cols(), 1);
				kInd.resize(kInd.rows() * kInd.cols(), 1);
				// std::cout<<"where8?"<<std::endl;
				auto ktemp = gens;
				for (auto mtemp = 0; mtemp < kInd.rows(); mtemp++)
				{
					if (kInd(mtemp, 0) == 1)
					{
						//G.row(i)(ktemp) = quadMatoffdiag(ktemp, 0);
						G.row(i)(ktemp) = quadMatoffdiag(mtemp, 0);
						ktemp++;
					}
				}
				// std::cout<<"where9?"<<std::endl;
			}
		}
		//cout << G;
		// generate new zonotope
		if (Qnonempty.sum() <= 1)
		{
			// std::cout<<"where10?"<<std::endl;
			Zonotope<Number> result(c, G.cwiseAbs().rowwise().sum());
			return result;
		}
		else
		{
			// std::cout<<"where11?"<<std::endl;
			Zonotope<Number> result(c, nonzeroFilter(G));
			return result;
		}
	}

	template <typename Number>
	Zonotope<Number> Zonotope<Number>::quadMapMixed(Zonotope<Number> Z2, std::vector<Matrix_t<Number>> Q)
	{
		// compute an over-approximation of the quadratic map
		// get matrix of zonotope
		Matrix_t<Number> Zmat1 = Z();

		Matrix_t<Number> Zmat2 = Z2.Z();

		int dimQ = Q.size();

		// init soulution (center + generator matrix)
		Matrix_t<Number> Z = Eigen::MatrixXd::Zero(dimQ, Zmat1.cols() * Zmat2.cols());

		// count empty matrices
		Matrix_t<int> Qnonempty;
		Qnonempty.setZero(dimQ, 1);

		// for each dimension compute center+generatpr elements
		for (auto i = 0; i < dimQ; i++)
		{
			Qnonempty(i) = Q[i].any();
			if (Qnonempty(i))
			{
				// pure quadratic evaluation

				// cout << Zmat1 << endl;
				// cout << Q[i] << endl;
				// cout << Zmat2 << endl;

				Matrix_t<Number> quadMat = Zmat1.adjoint() * Q[i] * Zmat2;
				// quadMat.conservativeResize(quadMat.rows() * quadMat.cols(), 1);
				// Z.row(i) = quadMat.adjoint();

				quadMat.resize(1, quadMat.rows() * quadMat.cols());
				Z.row(i) = quadMat;
			}
		}

		// cout << Z << endl;
		// cout << Qnonempty << endl;
		Matrix_t<Number> G = Z.middleCols(1, Z.cols() - 1);

		// generate new zonotope

		// cout << Qnonempty.sum() << endl;

		if (Qnonempty.sum() <= 1)
		{
			Zonotope<Number> result(Z.col(0), G.cwiseAbs().rowwise().sum());
			//cout << result;
			return result;
		}
		else
		{
			Zonotope<Number> result(Z.col(0), nonzeroFilter(G));
			return result;
		}
	}

	template <typename Number>
	std::vector<std::vector<Zonotope<Number>>> Zonotope<Number>::split()
	{
		// initialize
		// input number is one
		// split all dimensions
		Zonotope<Number>						   Z(interval());
		std::vector<std::vector<Zonotope<Number>>> result;
		for (auto i = 0; i < Z.center().rows(); i++)
		{
			// split one dimension
			std::vector<Zonotope<Number>> temp = splitOneDim(Z, i);
			result.push_back(temp);
		}
		return result;
	}

	template <typename Number>
	std::vector<Zonotope<Number>> Zonotope<Number>::split(int dim)
	{
		// initialize
		// no splitting halfspace is passed
		Zonotope<Number>			  Z(interval());
		std::vector<Zonotope<Number>> result = splitOneDim(Z, dim);
		return result;
	}

	template <typename Number>
	std::vector<Zonotope<Number>> Zonotope<Number>::splitOneDim(Zonotope<Number> Z, int dim)
	{
		// center and generator matrix
		Vector_t<Number> c = Z.center();
		Matrix_t<Number> G = Z.generators();

		// compute cnters of splitted parallelpiped
		Vector_t<Number> c1 = c - G.col(dim) / 2;
		Vector_t<Number> c2 = c + G.col(dim) / 2;

		// computer new set of generators
		Matrix_t<Number> Gnew = G;
		Gnew.col(dim)		  = Gnew.col(dim) / 2;

		// generate sp;otted parallepipeds
		Zonotope<Number>			  Zsplit1(c1, Gnew);
		Zonotope<Number>			  Zsplit2(c2, Gnew);
		std::vector<Zonotope<Number>> result;
		result.push_back(Zsplit1);
		result.push_back(Zsplit2);

		return result;
	}

	// template <typename Number>
	// std::ostream& operator<<(std::ostream& os, const Zonotope<Number>& zono)
	// {
	// 	os.precision(5);
	// 	if (zono.dimension() != 0)
	// 	{
	// 		os << "Dimension： " << zono.dimension() << std::endl;
	// 		os << "Center:" << std::endl;
	// 		os << "[";
	// 		for (auto i = 0; i < zono.center().rows()-1; i++)
	// 		{
	// 			os << zono.center()(i, 0) << "; ";
	// 		}
	// 		os << zono.center()(zono.center().rows()-1, 0) << "];";
	// 		os << std::endl;
	// 		os << "Genertors:"
	// 		   << "Number " << zono.generators().cols() << std::endl;

	// 		if( zono.generators().cols() != 0){
	// 			os << "[";
	// 			for (auto i = 0; i < zono.center().rows(); i++)
	// 			{
	// 				for (auto j = 0; j < zono.generators().cols()-1; j++)
	// 				{
	// 					os << zono.generators()(i, j) << ", ";
	// 				}
	// 				os << zono.generators()(i, zono.generators().cols()-1) << ";" << std::endl;
	// 			}
	// 			os << "];";
	// 		}
	// 	}
	// 	return os;
	// }

	template <typename Number>
	std::ostream& operator<<(std::ostream& os, const Zonotope<Number>& zono)
	{
		os.precision(5);

		if (zono.dimension() != 0)
		{
			Matrix_t<Number> Z(zono.dimension(), zono.generators().cols()+1);

			Z << zono.center(), zono.generators();
			
			os << Z;
		}

		return os;
	}

	// template <typename Number>
	// void RemoveColumn(Matrix_t<Number> & matrix, unsigned int colToRemove) {
 	// 	unsigned int numRows = matrix.rows();
  	// 	unsigned int numCols = matrix.cols() - 1;
 
  	// 	if( colToRemove < numCols ) {
    // 		matrix.block(0, colToRemove, numRows, numCols - colToRemove) =
    // 	  	matrix.block(0, colToRemove + 1, numRows, numCols - colToRemove);
  	// 	}
	
  	// 	matrix.conservativeResize(numRows,numCols);
	// }

	// template <typename Number>
	// void RemoveRow(Matrix_t<Number> & matrix, unsigned int rowToRemove) {
	// 	unsigned int numRows = matrix.rows() - 1;
	// 	unsigned int numCols = matrix.cols();
		
	// 	if( rowToRemove < numRows ) {
	// 	matrix.block(rowToRemove,0,numRows-rowToRemove,numCols) =
	// 		matrix.block(rowToRemove+1,0,numRows-rowToRemove,numCols);
	// 	}
		
	// 	matrix.conservativeResize(numRows,numCols);
	// }


	template <typename Number>
	void deleteZeros(Matrix_t<Number> & matrix) {
		Matrix_t<double> col_min = matrix.colwise().minCoeff();
		Matrix_t<double> col_max = matrix.colwise().maxCoeff();
		vector<unsigned int> zeroCol;
		unsigned int numCols = matrix.cols();
		for(auto i = 0; i < numCols; i++){
			if(col_min(0,i) == 0 && col_max(0,i) == 0){
				zeroCol.push_back(i);
			}
		}
		// unsigned int zeroColnum = zeroCol.size();
		// for(auto i = 0; i < zeroColnum; i++){
		// 	RemoveColumn(matrix,zeroCol[i]-i);
		// }
		RemoveColumn_Vec(matrix, zeroCol);
	}

	template <typename Number>
    Matrix_t<Number> convert2polygon(const Zonotope<Number> & zono){
		Vector_t<Number> c = zono.center();
		Matrix_t<Number> matrix = zono.generators();

		//delete zero generators
		deleteZeros(matrix);

		unsigned int colNum = matrix.cols();//obtain number of generators
		
		//obtain size of enclosing intervalhull of first two dimensions
		Number xmax = matrix.cwiseAbs().rowwise().sum()(0,0);
		Number ymax = matrix.cwiseAbs().rowwise().sum()(1,0);

		//Z with normalized direction: All generators pointing "up"
		Matrix_t<Number> matrixNorm = matrix;
		for(auto i = 0; i < colNum; i++){
			if(matrixNorm(1,i) < 0){
				matrixNorm(0,i) = matrixNorm(0,i) * -1;
				matrixNorm(1,i) = matrixNorm(1,i) * -1;
			}
		}
		vector<double> angles;
		for(auto i = 0; i < colNum; i++){
			angles.push_back(atan2(matrixNorm(1,i),matrixNorm(0,i)));
			if(angles[i] < 0){
				angles[i] = angles[i] + 2*M_PI;
			}
			//cout << angles[i] << " ";
		}
		cout << endl;
		//sort all generators by their angle
		vector<int> Index(colNum); //序号
		iota(Index.begin(),Index.end(),0);//递增赋值
		sort(Index.begin(),Index.end(),
			[&angles](double a,double b){return angles[a] < angles[b]; });//此处对数据判断，然后对序号排列

		//cumsum the generators in order of angle
		vector<Number> p_1(colNum+1);
		vector<Number> p_2(colNum+1);

		for(auto i = 0; i < colNum; i++){
			p_1[i+1] = p_1[i] + 2 * matrixNorm(0,Index[i]);
			p_2[i+1] = p_2[i] + 2 * matrixNorm(1,Index[i]);
		}

		Number p_1Max = *max_element(p_1.begin(), p_1.end());
		for(auto i = 0; i < colNum+1; i++){
			p_1[i] = p_1[i] + xmax - p_1Max;
			p_2[i] = p_2[i] - ymax;
		}

		//flip/mirror upper half to get lower half of zonotope (point symmetry)
		//consider center
		Matrix_t<Number> p(2,2*colNum+1);
		for(auto i = 0; i < colNum+1; i++){
			p(0,i) = p_1[i] + c(0);
			p(1,i) = p_2[i] + c(1);
		}
		for(auto i = colNum+1; i < 2*colNum+1; i++){
			p(0,i) = p_1[colNum] + p_1[0] - p_1[i-colNum] + c(0);
			p(1,i) = p_2[colNum] + p_2[0] - p_2[i-colNum] + c(1);
		}
		return p;
	}

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
                          Matrix_t<Number>& Gred)
	{
		c = Z.center();
		Matrix_t<Number> G = Z.generators();

		if(G.cols() != 0){

			//delete zero-length generators
			G = nonzeroFilter(G);

			//number of generators
			Number d = G.rows();
			Number nrOfGens = G.cols();

			//only reduce if zonotope order is greater than the desired order
			if(nrOfGens > d * order){

				//compute metric of generators
				Matrix_t<Number> G_1norm = G.cwiseAbs().colwise().sum();
				Matrix_t<Number> G_infnorm = G.cwiseAbs().colwise().maxCoeff();
				Matrix_t<Number> htemp = G_1norm - G_infnorm;

				std::vector<Number> h;
				for(int i = 0; i < htemp.cols(); i++){
					h.push_back(htemp(0,i));
				}
				//number of generators that are not reduced
				int nUnreduced = d * (order - 1);

				//number of generators that are reduced
				int nReduced = nrOfGens - nUnreduced;

				//pick generators with smallest h values to be reduced
				std::vector<size_t> indtemp = sort_indexes_ascend(h); 
				std::vector<size_t> ind(indtemp.begin(), indtemp.begin() + nReduced);

				Matrix_t<Number> Gredtemp(G.rows(),ind.size());

				for(int i = 0; i < ind.size(); i++){
					Gredtemp.col(i) = G.col(ind[i]);
				}

				// cout << ind << endl;

				//unreduced generators
				std::vector<size_t> idx(nrOfGens);
  				std::iota(idx.begin(), idx.end(), 0);
				
				sort(ind.begin(),ind.end());
				std::vector<size_t> indRemain = setdiff(idx, ind);
				Matrix_t<Number> Gunredtemp(G.rows(),indRemain.size());

				for(int i = 0; i < indRemain.size(); i++){
					Gunredtemp.col(i) = G.col(indRemain[i]);
				}

				Gred = Gredtemp;
				Gunred = Gunredtemp;
				
				// cout << Gred << endl;
				// cout << Gunred << endl;
			}else{
				Gunred = G;
			}
		}
	}

	template <typename Number>
	std::vector<Zonotope<Number>> Zono_border(const Zonotope<Number>& Zono){
		vector<Zonotope<Number>> border;
		int dim = Zono.center().rows();
		Matrix_t<Number> G = Zono.generators();
		Eigen::ColPivHouseholderQR<Matrix_t<Number>> QR_decompG(G);
		Number GRank = QR_decompG.rank();
		if(GRank != dim){
			cout << "No full rank !!!" << endl;
			return border;
		}
		vector<vector<Number>> com;
		vector<Number> temp;
		vector<Number> G_index(G.cols());
		iota(G_index.begin(), G_index.end(), 0);

		comb(G.cols(), temp, dim - 1, com);
		vector<vector<Number>> planes;
		// cout<< com ;
		for(int i = 0; i < com.size(); i++){
			if(isIncludePlanes(com[i], planes)){
				continue;
			}
			vector<Number> plane_indexs = com[i];
			Matrix_t<Number> B(dim, dim-1);
			for(int j = 0; j < com[i].size(); j++){
				B.col(j) = G.col(com[i][j]);
			}

			// cout << B << endl;

			// 选择的 generators 不列满秩
			Eigen::ColPivHouseholderQR<Matrix_t<Number>> QR_decompB(B);
			Number BRank = QR_decompB.rank();

			// cout << BRank << endl;

			if(BRank < dim - 1){
				continue;
			}

			//计算法向量
			Matrix_t<Number> normalVec(1, dim);
			Matrix_t<Number> tempB(dim, dim - 1);
			Vector_t<Number> cNewOne = Zono.center();
			Vector_t<Number> cNewOther = Zono.center();
			for(int j = 0; j < dim; j++){
				tempB = B;
				RemoveRow(tempB, j);
				normalVec(0, j) = pow(-1, j) * tempB.determinant();
			}

			// cout << normalVec << endl;

			vector<Number> remainG_index = setdiff(G_index, com[i]); 
			Matrix_t<Number> direction;
			for(int j = 0; j < remainG_index.size(); j++){
				direction = normalVec * G.col(remainG_index[j]);
				if(direction(0,0) == 0){
					plane_indexs.push_back(remainG_index[j]);
				}else if(direction(0,0) > 0){
					cNewOne += G.col(remainG_index[j]);
					cNewOther -= G.col(remainG_index[j]);
				}else{
					cNewOne -= G.col(remainG_index[j]);
					cNewOther += G.col(remainG_index[j]);
				}
			}
			sort(plane_indexs.begin(), plane_indexs.end());
			Matrix_t<Number> Gnew(G.rows(), plane_indexs.size());
			for(int j = 0; j < plane_indexs.size(); j++){
				Gnew.col(j) = G.col(plane_indexs[j]);
			}

			//两个对称边界
			Zonotope<Number> oneBorder(cNewOne, Gnew);
			Zonotope<Number> otherBorder(cNewOther, Gnew);
			planes.push_back(plane_indexs);
			border.push_back(oneBorder);
			border.push_back(otherBorder);
		}
		return border;
	}

	template <typename Number>
	Matrix_t<Number> Zono_border_matrix(const Zonotope<Number>& Zono){

		int dim = Zono.center().rows();
		Matrix_t<Number> G = Zono.generators();

		Matrix_t<Number> border(0,G.cols());

		Eigen::ColPivHouseholderQR<Matrix_t<Number>> QR_decompG(G);
		Number GRank = QR_decompG.rank();
		if(GRank != dim){
			cout << "No full rank !!!" << endl;
			return border;
		}
		vector<vector<Number>> com;
		vector<Number> temp;
		vector<Number> G_index(G.cols());
		iota(G_index.begin(), G_index.end(), 0);

		comb(G.cols(), temp, dim - 1, com);
		vector<vector<Number>> planes;
		// cout<< com ;
		for(int i = 0; i < com.size(); i++){
			if(isIncludePlanes(com[i], planes)){
				continue;
			}

			border.conservativeResize(border.rows() + 2, border.cols());

			for(int k = 0; k < com[i].size(); k++){
				int idx = com[i][k];
				border(border.rows() - 1, idx) = -2;
				border(border.rows() - 2, idx) = -2;
			}

			vector<Number> plane_indexs = com[i];
			Matrix_t<Number> B(dim, dim-1);
			for(int j = 0; j < com[i].size(); j++){
				B.col(j) = G.col(com[i][j]);
			}

			// cout << B << endl;

			// 选择的 generators 不列满秩
			Eigen::ColPivHouseholderQR<Matrix_t<Number>> QR_decompB(B);
			Number BRank = QR_decompB.rank();

			// cout << BRank << endl;

			if(BRank < dim - 1){
				continue;
			}

			//计算法向量
			Matrix_t<Number> normalVec(1, dim);
			Matrix_t<Number> tempB(dim, dim - 1);
			Vector_t<Number> cNewOne = Zono.center();
			Vector_t<Number> cNewOther = Zono.center();
			for(int j = 0; j < dim; j++){
				tempB = B;
				RemoveRow(tempB, j);
				normalVec(0, j) = pow(-1, j) * tempB.determinant();
			}

			// cout << normalVec << endl;

			vector<Number> remainG_index = setdiff(G_index, com[i]); 
			Matrix_t<Number> direction;
			for(int j = 0; j < remainG_index.size(); j++){
				direction = normalVec * G.col(remainG_index[j]);
				int idx = remainG_index[j];
				if(direction(0,0) == 0){
					plane_indexs.push_back(remainG_index[j]);
					border(border.rows() - 1, idx) = -2;
					border(border.rows() - 2, idx) = -2;
				}else if(direction(0,0) > 0){
					// cNewOne += G.col(remainG_index[j]);
					// cNewOther -= G.col(remainG_index[j]);
					border(border.rows() - 1, idx) = 1;
					border(border.rows() - 2, idx) = -1;
				}else{
					// cNewOne -= G.col(remainG_index[j]);
					// cNewOther += G.col(remainG_index[j]);
					border(border.rows() - 1, idx) = -1;
					border(border.rows() - 2, idx) = 1;
				}
			}

			sort(plane_indexs.begin(), plane_indexs.end());
			planes.push_back(plane_indexs);
		}
		return border;
	}

	template <typename Number>
	Matrix_t<Number> Zono_tiling_matrix(Zonotope<Number>& Zono){

		int dim = Zono.center().rows();

		Matrix_t<Number> G = Zono.generators();

		//确保 border 最右侧的方阵满秩
		Eigen::FullPivLU<Matrix_t<double>> lu(G);
		Matrix_t<Number> Q = lu.permutationQ();
		Matrix_t<Number> tempG = G * Q;

		Matrix_t<Number> Gnew(G.rows(), G.cols());
		Gnew << tempG.rightCols(tempG.cols() - dim), tempG.leftCols(dim);  

		Zono.set_generators(Gnew);

		// cout << G << endl;
		// cout << Gnew << endl;

		Matrix_t<Number> border = Zono_border_matrix(Zono);

		// cout << border << endl;

		Matrix_t<Number> tiling(0, border.cols());

		int col = 0;

		while(col + dim < border.cols()){
			
			// cout << border << endl;

			sortrows_oneCol(border, col);

			// cout << col << endl;
			// cout << border << endl;

			while(border.rows() > 0 && border(0,col) == -2 ){
				RemoveRow(border, 0);
			}

			if(border.rows() == 0 || border(0, col) != -1){
				break;
			}

			// if(border(0, col) != -1){
			// 	col++;
			// 	continue;
			// }

			int row = 0;

			while(border(row, col) == -1 && row < border.rows()){

				tiling.conservativeResize(tiling.rows() + 1, tiling.cols());
				tiling.row(tiling.rows() - 1) = border.row(row);

				tiling(tiling.rows() - 1, col) = -2;
				border(row,col) = 1;

				row++;
			}

			col++;
		}

		tiling.conservativeResize(tiling.rows() + 1, tiling.cols());
		tiling.row(tiling.rows() - 1) = border.row(0);

		for(int i = col; i < border.cols(); i++){
			tiling(tiling.rows() - 1, i) = -2;
		}

		return tiling;
	}

	template <typename Number>
	vector<Zonotope<Number>> Zono_tiling(const Zonotope<Number>& Zono){

		Zonotope<Number> ZonoNew = Zono;

		Matrix_t<Number> Matrix = Zono_tiling_matrix(ZonoNew);

		vector<Zonotope<Number>> res = Matrix2Zono(ZonoNew, Matrix);

		return res;
	}


	template <typename Number>
	vector<Zonotope<Number>> Matrix2Zono(const Zonotope<Number>& Zono, const Matrix_t<Number>& Matrix){

		vector<Zonotope<Number>> res;

		Matrix_t<Number> G = Zono.generators();
		for(int i = 0; i < Matrix.rows(); i++){

			Vector_t<Number> cNew = Zono.center();
			vector<int> Gindex;

			for(int j = 0; j < Matrix.cols(); j++){
				
				if(Matrix(i,j) == -2){

					Gindex.push_back(j);

				}else{

					cNew += Matrix(i,j) * G.col(j);

				}
			}

			Matrix_t<Number> Gnew(G.rows(), Gindex.size());

			for(int k = 0; k < Gindex.size(); k++){
				Gnew.col(k) = G.col(Gindex[k]);
			}

			Zonotope<Number> oneZono(cNew, Gnew);
			res.push_back(oneZono);
		}

		return res;
	}

	template <typename Number>
	Zonotope<Number> project(Zonotope<Number> const & Z, int dim){

		dim--;

		Zonotope<Number> res;

		res.set_center(Z.center().row(dim));
		res.set_generators(Z.generators().row(dim));

		return res;
	}
} // namespace reachSolver
