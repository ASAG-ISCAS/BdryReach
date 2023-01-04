/**
 * @file   LinearSys.cpp
 * @brief  LinearSys class
 * @author Yunjun Bai
 * @date April 2021
 * @version 1.0
 *
 * Reference:
 * CORA ../contDynamics/linearSys/all
 */

/**
 * @author Changyuan Zhao
 * @date  Nov 2021
 * @version 1.1
 *
 */

#include <dynamics/LinearSys.h>

namespace reachSolver
{

	// template class LinearSys<double>;

	/*****************************************************************************
 *                                                                           *
 *                           Constructors and Destructors                    *
 *                                                                           *
 *****************************************************************************/
	template <typename Number>
	LinearSys<Number>::LinearSys(Matrix_t<Number>& A, Matrix_t<Number>& B)
		: ContDynamics<Number>(std::string("linearSys"), A.rows(), B.cols(), 1)
		, A_(A)
		, B_(B)
		, c_(0, 0)
		, C_(Eigen::MatrixXd::Identity(1, 1))
		, D_(0, 0)
		, k_(0, 0)
	{}

	template <typename Number>
	LinearSys<Number>::LinearSys(Matrix_t<Number>& A, Matrix_t<Number>& B, Matrix_t<Number>& c)
		: ContDynamics<Number>(std::string("linearSys"), A.rows(), B.cols(), 1)
		, A_(A)
		, B_(B)
		, c_(c)
		, C_(Eigen::MatrixXd::Identity(1, 1))
		, D_(0, 0)
		, k_(0, 0)
	{}

	template <typename Number>
	LinearSys<Number>::LinearSys(Matrix_t<Number>& A, Matrix_t<Number>& B, Matrix_t<Number>& c, Matrix_t<Number>& C)
		: ContDynamics<Number>(std::string("linearSys"), A.rows(), B.cols(), C.rows())
		, A_(A)
		, B_(B)
		, c_(c)
		, C_(C)
		, D_(0, 0)
		, k_(0, 0)
	{}

	template <typename Number>
	LinearSys<Number>::LinearSys(
		Matrix_t<Number>& A, Matrix_t<Number>& B, Matrix_t<Number>& c, Matrix_t<Number>& C, Matrix_t<Number>& D)
		: ContDynamics<Number>(std::string("linearSys"), A.rows(), B.cols(), C.rows())
		, A_(A)
		, B_(B)
		, c_(c)
		, C_(C)
		, D_(D)
		, k_(0, 0)
	{}

	template <typename Number>
	LinearSys<Number>::LinearSys(Matrix_t<Number>& A,
								 Matrix_t<Number>& B,
								 Matrix_t<Number>& c,
								 Matrix_t<Number>& C,
								 Matrix_t<Number>& D,
								 Matrix_t<Number>& k)
		: ContDynamics<Number>(std::string("linearSys"), A.rows(), B.cols(), C.rows())
		, A_(A)
		, B_(B)
		, c_(c)
		, C_(C)
		, D_(D)
		, k_(k)
	{}

	template <typename Number>
	LinearSys<Number>::LinearSys(std::string name, Matrix_t<Number>& A, Matrix_t<Number>& B)
		: ContDynamics<Number>(name, A.rows(), B.cols(), 1)
		, A_(A)
		, B_(B)
		, c_(0, 0)
		, C_(Eigen::MatrixXd::Identity(1, 1))
		, D_(0, 0)
		, k_(0, 0)
	{}

	template <typename Number>
	LinearSys<Number>::LinearSys(std::string name, Matrix_t<Number>& A, Matrix_t<Number>& B, Matrix_t<Number>& c)
		: ContDynamics<Number>(name, A.rows(), B.cols(), 1)
		, A_(A)
		, B_(B)
		, c_(c)
		, C_(Eigen::MatrixXd::Identity(1, 1))
		, D_(0, 0)
		, k_(0, 0)
	{}

	template <typename Number>
	LinearSys<Number>::LinearSys(
		std::string name, Matrix_t<Number>& A, Matrix_t<Number>& B, Matrix_t<Number>& c, Matrix_t<Number>& C)
		: ContDynamics<Number>(name, A.rows(), B.cols(), C.rows())
		, A_(A)
		, B_(B)
		, c_(c)
		, C_(C)
		, D_(0, 0)
		, k_(0, 0)
	{}

	template <typename Number>
	LinearSys<Number>::LinearSys(std::string	   name,
								 Matrix_t<Number>& A,
								 Matrix_t<Number>& B,
								 Matrix_t<Number>& c,
								 Matrix_t<Number>& C,
								 Matrix_t<Number>& D)
		: ContDynamics<Number>(name, A.rows(), B.cols(), C.rows())
		, A_(A)
		, B_(B)
		, c_(c)
		, C_(C)
		, D_(D)
		, k_(0, 0)
	{}

	template <typename Number>
	LinearSys<Number>::LinearSys(std::string	   name,
								 Matrix_t<Number>& A,
								 Matrix_t<Number>& B,
								 Matrix_t<Number>& c,
								 Matrix_t<Number>& C,
								 Matrix_t<Number>& D,
								 Matrix_t<Number>& k)
		: ContDynamics<Number>(name, A.rows(), B.cols(), C.rows())
		, A_(A)
		, B_(B)
		, c_(c)
		, C_(C)
		, D_(D)
		, k_(k)
	{}

	template <typename Number>
	LinearSys<Number>::~LinearSys()
	{}

	/*****************************************************************************
 *                                                                           *
 *                       Compute Reachable Set                               *
 *                                                                           *
 *****************************************************************************/

	template <typename Number>
	LinearReachableSet<Number> LinearSys<Number>::initReach(Zonotope<Number>& Rinit, ReachOptions<Number>& options)
	{
		if (options.usekrylovError() == 0)
		{
			// compute in Krylov space
			return initReach_Krylov(Rinit, options);
		}
		else
		{
			return initReach_Euclidean(Rinit, options);
		}
	}

	template <typename Number>
	polyLinearReachableSet<Number> LinearSys<Number>::initReach(polyZonotope<Number>&	  Rinit,
																polyReachOptions<Number>& options)
	{
		if (options.usekrylovError() == 0)
		{
			// compute in Krylov space
			return initReach_Krylov(Rinit, options);
		}
		else
		{
			return initReach_Euclidean(Rinit, options);
		}
	}

	template <typename Number>
	LinearReachableSet<Number> LinearSys<Number>::initReach_Krylov(Zonotope<Number>&	 Rinit,
																   ReachOptions<Number>& options)
	{
		return LinearReachableSet<Number>();
	}

	template <typename Number>
	polyLinearReachableSet<Number> LinearSys<Number>::initReach_Krylov(polyZonotope<Number>&	 Rinit,
																	   polyReachOptions<Number>& options)
	{
		return polyLinearReachableSet<Number>();
	}

	template <typename Number>
	void LinearSys<Number>::exponential(ReachOptions<Number>& options)
	{
		// load data from object/options structure
		Matrix_t<Number> A			 = A_;
		Matrix_t<Number> A_abs		 = A_.cwiseAbs();
		size_t			 taylorTerms = options.taylor_terms();
		size_t			 dim_local	 = ContDynamics<Number>::dim();
		Number*			 factors	 = options.factor();

		// initialize
		std::vector<Matrix_t<Number>> Apower	 = std::vector<Matrix_t<Number>>(taylorTerms + 1);
		Apower[0]								 = A;
		std::vector<Matrix_t<Number>> Apower_abs = std::vector<Matrix_t<Number>>(taylorTerms + 1);
		Apower_abs[0]							 = A_abs;
		Matrix_t<Number> M						 = Eigen::MatrixXd::Identity(dim_local, dim_local);
		// compute powers for each term and sum of these
		for (size_t i = 0; i < taylorTerms; i++)
		{
			Apower[i + 1]	  = Apower[i] * A;
			Apower_abs[i + 1] = Apower_abs[i] * A_abs;
			M				  = M + Apower_abs[i] * factors[i];
		}

		// determine error due to finite Taylor series
		Matrix_t<Number> tmp	 = A_abs * options.time_step();
		Matrix_t<Number> tmp_exp = tmp.exp();
		Matrix_t<Number> W		 = tmp_exp - M;

		// compute absolute value of W for numerical stability
		W						  = W.cwiseAbs();
		Matrix_t<Number>	   W2 = -1. * W;
		IntervalMatrix<Number> E  = setIntervalMatrix(W2, W);

		// write to object structure
		taylor_.powers = Apower;
		taylor_.error  = E;
	}

	template <typename Number>
	void LinearSys<Number>::exponential(polyReachOptions<Number>& options)
	{
		// load data from object/options structure
		Matrix_t<Number> A			 = A_;
		Matrix_t<Number> A_abs		 = A_.cwiseAbs();
		size_t			 taylorTerms = options.taylor_terms();
		size_t			 dim_local	 = ContDynamics<Number>::dim();
		Number*			 factors	 = options.factor();

		// initialize
		std::vector<Matrix_t<Number>> Apower	 = std::vector<Matrix_t<Number>>(taylorTerms + 1);
		Apower[0]								 = A;
		std::vector<Matrix_t<Number>> Apower_abs = std::vector<Matrix_t<Number>>(taylorTerms + 1);
		Apower_abs[0]							 = A_abs;
		Matrix_t<Number> M						 = Eigen::MatrixXd::Identity(dim_local, dim_local);
		// compute powers for each term and sum of these
		for (size_t i = 0; i < taylorTerms; i++)
		{
			Apower[i + 1]	  = Apower[i] * A;
			Apower_abs[i + 1] = Apower_abs[i] * A_abs;
			M				  = M + Apower_abs[i] * factors[i];
		}

		// determine error due to finite Taylor series
		Matrix_t<Number> tmp	 = A_abs * options.time_step();
		Matrix_t<Number> tmp_exp = tmp.exp();
		Matrix_t<Number> W		 = tmp_exp - M;

		// compute absolute value of W for numerical stability
		W						  = W.cwiseAbs();
		Matrix_t<Number>	   W2 = -1. * W;
		IntervalMatrix<Number> E  = setIntervalMatrix(W2, W);

		// write to object structure
		taylor_.powers = Apower;
		taylor_.error  = E;
	}

	template <typename Number>
	void LinearSys<Number>::tie(ReachOptions<Number>& options)
	{
		// load data from object/options structure
		std::vector<Matrix_t<Number>> Apower	  = taylor_.powers;
		size_t						  taylorTerms = options.taylor_terms();
		Number*						  rbyfac	  = options.factor();
		size_t						  dim_local	  = ContDynamics<Number>::dim();

		// initialize Asum
		Matrix_t<Number> Asum_pos = Eigen::MatrixXd::Zero(dim_local, dim_local);
		Matrix_t<Number> Asum_neg = Eigen::MatrixXd::Zero(dim_local, dim_local);

		for (size_t i = 2; i <= taylorTerms; i++)
		{
			// compute factor
			double exp1	  = -(double) i / (i - 1);
			double exp2	  = -(double) 1 / (i - 1);
			double factor = (std::pow(i, exp1) - std::pow(i, exp2)) * rbyfac[i - 1];
			// init Apos, Aneg
			Matrix_t<Number> Apos = Eigen::MatrixXd::Zero(dim_local, dim_local);
			Matrix_t<Number> Aneg = Eigen::MatrixXd::Zero(dim_local, dim_local);

			// obtain positive and negative parts
			for (size_t j = 0; j < dim_local; j++)
			{
				for (size_t k = 0; k < dim_local; k++)
				{
					if (Apower[i - 1](j, k) > 0)
					{
						Apos(j, k) = Apower[i - 1](j, k);
					}
					else
					{
						Aneg(j, k) = Apower[i - 1](j, k);
					}
				}
			}
			// compute powers; factor is always negative
			Asum_pos = Asum_pos + factor * Aneg;
			Asum_neg = Asum_neg + factor * Apos;
		}
		// instantiate interval matrix
		IntervalMatrix<Number> Asum = setIntervalMatrix(Asum_neg, Asum_pos);

		// write to object structure
		taylor_.F = Asum + taylor_.error;
		// taylor_.F.inf = Asum.inf+taylor_.error.inf;
		// taylor_.F.sup = Asum.sup+taylor_.error.sup;
	}

	template <typename Number>
	void LinearSys<Number>::tie(polyReachOptions<Number>& options)
	{
		// load data from object/options structure
		std::vector<Matrix_t<Number>> Apower	  = taylor_.powers;
		size_t						  taylorTerms = options.taylor_terms();
		Number*						  rbyfac	  = options.factor();
		size_t						  dim_local	  = ContDynamics<Number>::dim();

		// initialize Asum
		Matrix_t<Number> Asum_pos = Eigen::MatrixXd::Zero(dim_local, dim_local);
		Matrix_t<Number> Asum_neg = Eigen::MatrixXd::Zero(dim_local, dim_local);

		for (size_t i = 2; i <= taylorTerms; i++)
		{
			// compute factor
			double exp1	  = -(double) i / (i - 1);
			double exp2	  = -(double) 1 / (i - 1);
			double factor = (std::pow(i, exp1) - std::pow(i, exp2)) * rbyfac[i - 1];
			// init Apos, Aneg
			Matrix_t<Number> Apos = Eigen::MatrixXd::Zero(dim_local, dim_local);
			Matrix_t<Number> Aneg = Eigen::MatrixXd::Zero(dim_local, dim_local);

			// obtain positive and negative parts
			for (size_t j = 0; j < dim_local; j++)
			{
				for (size_t k = 0; k < dim_local; k++)
				{
					if (Apower[i - 1](j, k) > 0)
					{
						Apos(j, k) = Apower[i - 1](j, k);
					}
					else
					{
						Aneg(j, k) = Apower[i - 1](j, k);
					}
				}
			}
			// compute powers; factor is always negative
			Asum_pos = Asum_pos + factor * Aneg;
			Asum_neg = Asum_neg + factor * Apos;
		}
		// instantiate interval matrix
		IntervalMatrix<Number> Asum = setIntervalMatrix(Asum_neg, Asum_pos);

		// write to object structure
		taylor_.F = Asum + taylor_.error;
	} // tie poly

	template <typename Number>
	void LinearSys<Number>::inputSolution(ReachOptions<Number>& options)
	{
		// std::cout << "this inputSolution" << std::endl;
		// set of possible inputs
		Zonotope<Number> V = options.U() * B_(0, 0);

		options.set_isRV(true);
		if (V.center() == Eigen::MatrixXd::Zero(this->dim(), 1) && V.Z().cols() == 1)
		{
			options.set_isRV(false);
		}

		// compute vTrans
		Zonotope<Number> vTrans = options.uTrans_lin() * B_(0, 0);
		// consider constant input
		if (c_ != Matrix_t<Number>(0, 0))
		{
			vTrans = vTrans + c_;
		}
		Matrix_t<Number>			  A			  = A_;
		std::vector<Matrix_t<Number>> Apower	  = taylor_.powers;
		IntervalMatrix<Number>		  E			  = taylor_.error;
		size_t						  taylorTerms = options.taylor_terms();
		double						  r			  = options.time_step();
		size_t						  dim		  = std::max(A.cols(), A.rows());
		Number*						  factors	  = options.factor();
		Zonotope<Number>			  inputSolV;
		Matrix_t<Number>			  Asum;

		if (options.isRV())
		{
			// init Vsum
			Zonotope<Number> Vsum = V * r;
			Asum				  = r * Eigen::MatrixXd::Identity(dim, dim);
			for (size_t i = 0; i < taylorTerms; i++)
			{
				Matrix_t<Number> temp	 = Apower[i] * factors[i + 1];
				Zonotope<Number> tempzzz = V * temp;
				Vsum					 = Vsum + tempzzz;
				Asum					 = Asum + Apower[i] * factors[i + 1];
			}
			inputSolV = Vsum + V * E * r;
		}
		else
		{
			// only Asum since V == origin(0)
			Asum = r * Eigen::MatrixXd::Identity(dim, dim);
			// computer higher order terms
			for (auto i = 0; i < taylorTerms; i++)
			{
				Asum = Asum + Apower[i] * factors[i + 1];
			}
		}

		// compute solution due to constant input
		IntervalMatrix<Number> eAtInt = setIntervalMatrix(Asum) + E * r;
		/*Matrix_t<Number> eAtInt

  for(auto i=0;i<Asum_neg.rows();i++){
      for(auto j=0;j<Asum_neg.cols();j++){
          Asum(i,j).setLeftBound(Asum_neg(i,j));
          Asum(i,j).setRightBound(Asum_pos(i,j));
      }
  }

  eAtInt.inf = Asum+E.inf*r;
  eAtInt.sup = Asum+E.sup*r;*/
		// std::cout<<"333333"<<std::endl;
		Zonotope<Number> inputSolVtrans = vTrans * eAtInt;

		// compute additional uncertainty if origin is not contained in input set
		Zonotope<Number> inputCorr;
		if (options.originContained() != 0)
		{
			Matrix_t<Number> inputCorr = Eigen::MatrixXd::Zero(dim, 1);
		}
		else
		{
			// compute inputF
			inputTie(options);
			IntervalMatrix<Number> inputF = taylor_.inputF;
			inputCorr					  = vTrans * inputF;
		}

		// write to object structure
		taylor_.V = V;
		// Matrix_t<Number> SolV_Z = Matrix_t<Number>(inputSolV.center().rows(),
		// inputSolV.center().cols()+inputSolV.generators().cols());

		// SolV_Z << inputSolV.center(), inputSolV.generators();

		// cout<<"Solvez"<<endl;
		// cout<<SolV_Z<<endl;

		// cout<<SolV_Z.any()<<endl;

		// cout<<inputSolV.center().any()<<endl;

		if (options.isRV() && (inputSolV.center().any() || inputSolV.generators().any()))
		{
			taylor_.RV = inputSolV;
		}
		else
		{
			taylor_.RV = Zonotope<Number>(dim);
		}
		// std::cout<<"4444444444444"<<std::endl;
		// Matrix_t<Number> SolVtrans_Z =
		// Matrix_t<Number>(inputSolVtrans.center().size(), Eigen::Dynamic);
		// SolVtrans_Z << inputSolVtrans.center(), inputSolV.generators();
		if (inputSolVtrans.center().any() || inputSolVtrans.generators().any())
		{
			taylor_.Rtrans = inputSolVtrans;
		}
		else
		{
			taylor_.Rtrans = Zonotope<Number>(dim);
		}
		taylor_.inputCorr = inputCorr;
		taylor_.eAtInt	  = eAtInt;

		// delete[] factors;
	}

	template <typename Number>
	void LinearSys<Number>::inputSolution(polyReachOptions<Number>& options)
	{
		// std::cout << "this inputSolution" << std::endl;
		// set of possible inputs
		Zonotope<Number> V = options.U() * B_(0, 0);

		options.set_isRV(true);
		if (V.center() == Eigen::MatrixXd::Zero(this->dim(), 1) && V.Z().cols() == 1)
		{
			options.set_isRV(false);
		}

		// compute vTrans
		Zonotope<Number> vTrans = options.uTrans_lin() * B_(0, 0);

		// cout << V << endl;
		// cout << vTrans << endl;

		// consider constant input
		if (c_ != Matrix_t<Number>(0, 0))
		{
			vTrans = vTrans + c_;
		}
		Matrix_t<Number>			  A			  = A_;
		std::vector<Matrix_t<Number>> Apower	  = taylor_.powers;
		IntervalMatrix<Number>		  E			  = taylor_.error;
		size_t						  taylorTerms = options.taylor_terms();
		double						  r			  = options.time_step();
		size_t						  dim		  = std::max(A.cols(), A.rows());
		Number*						  factors	  = options.factor();
		Zonotope<Number>			  inputSolV;
		Matrix_t<Number>			  Asum;

		if (options.isRV())
		{
			// init Vsum
			Zonotope<Number> Vsum = V * r;
			Asum				  = r * Eigen::MatrixXd::Identity(dim, dim);
			for (size_t i = 0; i < taylorTerms; i++)
			{
				Matrix_t<Number> temp	 = Apower[i] * factors[i + 1];
				Zonotope<Number> tempzzz = V * temp;
				Vsum					 = Vsum + tempzzz;
				Asum					 = Asum + Apower[i] * factors[i + 1];
			}
			//compute overall solution
			inputSolV = Vsum + V* E* r;
		}
		else
		{
			// only Asum since V == origin(0)
			Asum = r * Eigen::MatrixXd::Identity(dim, dim);
			// computer higher order terms
			for (auto i = 0; i < taylorTerms; i++)
			{
				Asum = Asum + Apower[i] * factors[i + 1];
			}
		}

		// compute solution due to constant input
		IntervalMatrix<Number> eAtInt = setIntervalMatrix(Asum) + E * r;

		Zonotope<Number> inputSolVtrans = vTrans * eAtInt;

		// cout << eAtInt << endl;
		// cout << inputSolVtrans << endl;

		// compute additional uncertainty if origin is not contained in input set
		Zonotope<Number> inputCorr;
		if (options.originContained() != 0)
		{
			Matrix_t<Number> inputCorr = Eigen::MatrixXd::Zero(dim, 1);
		}
		else
		{
			// compute inputF
			inputTie(options);
			IntervalMatrix<Number> inputF = taylor_.inputF;
			inputCorr					  = vTrans * inputF;
		}

		// write to object structure
		taylor_.V = V;

		if (options.isRV() && (inputSolV.center().any() || inputSolV.generators().any()))
		{
			taylor_.RV = inputSolV;
		}
		else
		{
			taylor_.RV = Zonotope<Number>(dim);
		}
		if (inputSolVtrans.center().any() || inputSolVtrans.generators().any())
		{
			taylor_.Rtrans = inputSolVtrans;
		}
		else
		{
			taylor_.Rtrans = Zonotope<Number>(dim);
		}
		taylor_.inputCorr = inputCorr;
		taylor_.eAtInt	  = eAtInt;

		// delete[] factors;
	}

	template <typename Number>
	void LinearSys<Number>::inputTie(ReachOptions<Number>& options)
	{
		// load data from object structure
		std::vector<Matrix_t<Number>> Apower	  = taylor_.powers;
		IntervalMatrix<Number>		  E			  = taylor_.error;
		size_t						  taylorTerms = options.taylor_terms();
		double						  r			  = options.time_step();
		size_t						  dim		  = std::max(A_.cols(), A_.rows());

		// initialize Asum
		Matrix_t<Number> Asum_pos = Eigen::MatrixXd::Zero(dim, dim);
		Matrix_t<Number> Asum_neg = Eigen::MatrixXd::Zero(dim, dim);
		for (size_t i = 2; i <= taylorTerms + 1; i++)
		{
			// compute factor
			double exp1	  = -(double) i / (i - 1);
			double exp2	  = -(double) 1 / (i - 1);
			double factor = (std::pow(i, exp1) - std::pow(i, exp2)) * options.factor()[i - 1];

			// init Apos, Aneg
			Matrix_t<Number> Apos = Eigen::MatrixXd::Zero(dim, dim);
			Matrix_t<Number> Aneg = Eigen::MatrixXd::Zero(dim, dim);
			// obtain positive and negative parts
			for (size_t j = 0; j < dim; j++)
			{
				for (size_t k = 0; k < dim; k++)
				{
					if (Apower[i - 2](j, k) > 0)
					{
						Apos(j, k) = Apower[i - 2](j, k);
					}
					else
					{
						Aneg(j, k) = Apower[i - 2](j, k);
					}
				}
			}

			// compute powers; factor is always negative
			Asum_pos = Asum_pos + factor * Aneg;
			Asum_neg = Asum_neg + factor * Apos;
		}
		// instantiate interval matrix
		IntervalMatrix<Number> Asum = setIntervalMatrix(Asum_neg, Asum_pos);

		// compute error due to finite Taylor series according to internal document
		// "Input Error Bounds in Reachability Analysis"
		IntervalMatrix<Number> Einput = E * r;

		/*Einput.inf = E.inf*r;
  Einput.sup = E.sup*r;*/

		// write to object structure
		taylor_.inputF = Asum + Einput;
		/*taylor_.inputF.inf = Asum.inf+Einput.inf;
  taylor_.inputF.sup = Asum.sup+Einput.sup;*/
	}

	template <typename Number>
	void LinearSys<Number>::inputTie(polyReachOptions<Number>& options)
	{
		// load data from object structure
		std::vector<Matrix_t<Number>> Apower	  = taylor_.powers;
		IntervalMatrix<Number>		  E			  = taylor_.error;
		size_t						  taylorTerms = options.taylor_terms();
		double						  r			  = options.time_step();
		size_t						  dim		  = std::max(A_.cols(), A_.rows());

		// initialize Asum
		Matrix_t<Number> Asum_pos = Eigen::MatrixXd::Zero(dim, dim);
		Matrix_t<Number> Asum_neg = Eigen::MatrixXd::Zero(dim, dim);
		for (size_t i = 2; i <= taylorTerms + 1; i++)
		{
			// compute factor
			double exp1	  = -(double) i / (i - 1);
			double exp2	  = -(double) 1 / (i - 1);
			double factor = (std::pow(i, exp1) - std::pow(i, exp2)) * options.factor()[i - 1];

			// init Apos, Aneg
			Matrix_t<Number> Apos = Eigen::MatrixXd::Zero(dim, dim);
			Matrix_t<Number> Aneg = Eigen::MatrixXd::Zero(dim, dim);
			// obtain positive and negative parts
			for (size_t j = 0; j < dim; j++)
			{
				for (size_t k = 0; k < dim; k++)
				{
					if (Apower[i - 2](j, k) > 0)
					{
						Apos(j, k) = Apower[i - 2](j, k);
					}
					else
					{
						Aneg(j, k) = Apower[i - 2](j, k);
					}
				}
			}

			// compute powers; factor is always negative
			Asum_pos = Asum_pos + factor * Aneg;
			Asum_neg = Asum_neg + factor * Apos;
		}
		// instantiate interval matrix
		IntervalMatrix<Number> Asum = setIntervalMatrix(Asum_neg, Asum_pos);

		// compute error due to finite Taylor series according to internal document
		// "Input Error Bounds in Reachability Analysis"
		IntervalMatrix<Number> Einput = E * r;

		/*Einput.inf = E.inf*r;
  Einput.sup = E.sup*r;*/

		// write to object structure
		taylor_.inputF = Asum + Einput;
	}

	template <typename Number>
	LinearReachableSet<Number> LinearSys<Number>::initReach_Euclidean(Zonotope<Number>&		Rinit,
																	  ReachOptions<Number>& options)
	{
		// std::cout << "this initReach_Euclidean" << std::endl;
		// compute exponential matrix
		exponential(options);
		// std::cout << "end exponential" << std::endl;
		// compute time interval error (tie)
		tie(options);
		// std::cout << "end tie" << std::endl;
		// compute reachable set due to input
		inputSolution(options);
		// std::cout << "end inputsolution" << std::endl;
		// change the time step
		taylor_.timeStep = options.time_step();
		// compute reachable set of first time interval
		Matrix_t<Number> eAt = (A_ * options.time_step()).exp();
		// save data to object structure
		taylor_.eAt						 = eAt;
		IntervalMatrix<Number> F		 = taylor_.F;
		Zonotope<Number>	   RV		 = taylor_.RV;
		Zonotope<Number>	   inputCorr = taylor_.inputCorr;
		Zonotope<Number>	   Rtrans	 = taylor_.Rtrans;

		// first time step homogeneous solution
		Zonotope<Number> Rhom_tp = Rinit * eAt + Rtrans;

		Zonotope<Number> Rhom = Rinit.enclose(Rhom_tp) + Rinit * F + inputCorr;
		// reduce zonotope
		Rhom.Reduce(options.zonotope_order());
		Rhom_tp.Reduce(options.zonotope_order());
		RV.Reduce(options.zonotope_order());

		// save homogeneous and particulate solution
		options.set_Rhom(Rhom);
		options.set_Rhom_tp(Rhom_tp);
		options.set_Raux(RV);
		// if (strcmp(options.linAlg,'wrapping-free') == 0){
		//     options.Rpar=interval(RV);
		// }else{
		options.set_Rpar(RV);
		// }
		options.set_Rtrans(taylor_.Rtrans);

		Zonotope<Number> Rtotal;
		Zonotope<Number> Rtotal_tp;
		// total solution
		// if isa(Rinit,'mptPolytope'){
		//     // convert zonotopes to polytopes
		//     Radd=mptPolytope(RV);
		//     Rtotal=Rhom+Radd;
		//     Rtotal_tp=Rhom_tp+Radd;
		// }else{
		// original computation
		Rtotal	  = Rhom + RV;
		Rtotal_tp = Rhom_tp + RV;
		// }

		// write results to reachable set struct Rfirst
		LinearReachableSet<Number> Rfirst = LinearReachableSet<Number>();
		Rfirst.set_time_point(Rtotal_tp);
		Rfirst.set_time_interval(Rtotal);
		// std::cout << "end initReach_Euclidean" << std::endl;
		return Rfirst;
	}

	template <typename Number>
	polyLinearReachableSet<Number> LinearSys<Number>::initReach_Euclidean(polyZonotope<Number>&		Rinit,
																		  polyReachOptions<Number>& options)
	{
		std::cout << "this initReach_Euclidean" << std::endl;
		// compute exponential matrix
		exponential(options);
		// std::cout << "end exponential" << std::endl;
		// compute time interval error (tie)
		tie(options);
		// std::cout << "end tie" << std::endl;
		// compute reachable set due to input
		inputSolution(options);
		// std::cout << "end inputsolution" << std::endl;
		// change the time step
		taylor_.timeStep = options.time_step();
		// compute reachable set of first time interval
		Matrix_t<Number> eAt = (A_ * options.time_step()).exp();
		// save data to object structure
		taylor_.eAt						 = eAt;
		IntervalMatrix<Number> F		 = taylor_.F;
		Zonotope<Number>	   RV		 = taylor_.RV;
		Zonotope<Number>	   inputCorr = taylor_.inputCorr;
		Zonotope<Number>	   Rtrans	 = taylor_.Rtrans;
		// first time step homogeneous solution

		// cout << Rinit << endl;
		// cout << eAt << endl;
		// cout << Rtrans << endl;

		// cout << Rinit << endl;
		// cout << F << endl;
		// cout << inputCorr << endl;
		
		polyZonotope<Number> Rhom_tp = Rinit * eAt + Rtrans;
		polyZonotope<Number> Rhom	 = Rinit.enclose(Rhom_tp) + zonotope(Rinit) * F + inputCorr;

		// reduce zonotope

		// cout << Rinit.enclose(Rhom_tp) << endl;
		// cout << Rhom << endl;
		// cout << Rhom_tp << endl;
		// cout << RV << endl;

		Rhom.Reduce(options.zonotope_order());
		// std::cout << endl << "HERE GOOD" << endl << std::endl;
		Rhom_tp.Reduce(options.zonotope_order());
		RV.Reduce(options.zonotope_order());

		// cout << Rhom << endl;
		// cout << Rhom_tp << endl;
		// cout << RV << endl;

		// std::cout << endl << "HERE GOOD" << endl << std::endl;

		// save homogeneous and particulate solution
		options.set_Rhom(Rhom);
		options.set_Rhom_tp(Rhom_tp);
		options.set_Raux(RV);
		// if (strcmp(options.linAlg,'wrapping-free') == 0){
		//     options.Rpar=interval(RV);
		// }else{
		options.set_Rpar(RV);
		// }
		options.set_Rtrans(taylor_.Rtrans);

		polyZonotope<Number> Rtotal;
		polyZonotope<Number> Rtotal_tp;
		// total solution
		// if isa(Rinit,'mptPolytope'){
		//     // convert zonotopes to polytopes
		//     Radd=mptPolytope(RV);
		//     Rtotal=Rhom+Radd;
		//     Rtotal_tp=Rhom_tp+Radd;
		// }else{
		// original computation
		Rtotal	  = Rhom + RV;
		Rtotal_tp = Rhom_tp + RV;
		// }

		// write results to reachable set struct Rfirst*/
		polyLinearReachableSet<Number> Rfirst = polyLinearReachableSet<Number>();
		Rfirst.set_time_point(Rtotal_tp);
		Rfirst.set_time_interval(Rtotal);
		std::cout << "end initReach_Euclidean" << std::endl;
		return Rfirst;
	}

	template <typename Number>
	Zonotope<Number> LinearSys<Number>::error_solution(ReachOptions<Number>& options,
													   Zonotope<Number>		 Vdyn,
													   Zonotope<Number>		 Vstat)
	{
		Zonotope<Number> errorStat;
		if (Vstat.Empty())
		{
			errorStat = Zonotope<Number>(this->dim());
		}
		else
		{
			errorStat = Vstat * taylor_.eAtInt;
		}
		// load data from object/options structure
		std::vector<Matrix_t<Number>> Apower	  = taylor_.powers;
		IntervalMatrix<Number>		  E			  = taylor_.error;
		size_t						  taylorTerms = options.taylor_terms();
		double						  r			  = options.time_step();
		Number*						  factors	  = options.factor();
		// initialize Asum

		Zonotope<Number> Asum = Vdyn * r;
		for (size_t i = 0; i < taylorTerms; i++)
		{
			// compute powers
			Matrix_t<Number> temp	 = factors[i + 1] * Apower[i];
			Zonotope<Number> ApowerV = Vdyn * temp;
			// compute sums
			Asum = Asum + ApowerV;
		}
		// get error due to finite Taylor series
		Zonotope<Number> F;
		// if(isa(Vdyn,'zonoBundle')){
		//     F=E*Vdyn.Z{1}*r;
		// }else{
		F = Vdyn * E * r;
		// }
		// Compute error solution (dyn. + stat.)
		Zonotope<Number> Rerror = Asum + F + errorStat;

		return Rerror;
	}

	template <typename Number>
	Zonotope<Number> LinearSys<Number>::error_solution(polyReachOptions<Number>& options,
													   Zonotope<Number>			 Vdyn,
													   Zonotope<Number>			 Vstat)
	{
		Zonotope<Number> errorStat;
		if (Vstat.Empty())
		{
			errorStat = Zonotope<Number>(this->dim());
		}
		else
		{
			errorStat = Vstat * taylor_.eAtInt;
		}
		// load data from object/options structure
		std::vector<Matrix_t<Number>> Apower	  = taylor_.powers;
		IntervalMatrix<Number>		  E			  = taylor_.error;
		size_t						  taylorTerms = options.taylor_terms();
		double						  r			  = options.time_step();
		Number*						  factors	  = options.factor();
		// initialize Asum

		Zonotope<Number> Asum = Vdyn * r;
		for (size_t i = 0; i < taylorTerms; i++)
		{
			// compute powers
			Matrix_t<Number> temp	 = factors[i + 1] * Apower[i];
			Zonotope<Number> ApowerV = Vdyn * temp;
			// compute sums
			Asum = Asum + ApowerV;
		}
		// get error due to finite Taylor series
		Zonotope<Number> F;
		// if(isa(Vdyn,'zonoBundle')){
		//     F=E*Vdyn.Z{1}*r;
		// }else{
		F = Vdyn * E * r;
		// }
		// Compute error solution (dyn. + stat.)
		Zonotope<Number> Rerror = Asum + F + errorStat;

		return Rerror;
	}

	template <typename Number>
	polyZonotope<Number> LinearSys<Number>::error_solution(polyReachOptions<Number>& options,
														   Zonotope<Number>			 Vdyn,
														   polyZonotope<Number>		 Vstat)
	{
		polyZonotope<Number> errorStat;

		// cout << taylor_.eAtInt << endl;

		if (Vstat.Empty())
		{
			errorStat = polyZonotope<Number>(this->dim());
		}
		else
		{
			errorStat = Vstat * taylor_.eAtInt;
		}
		// load data from object/options structure
		std::vector<Matrix_t<Number>> Apower	  = taylor_.powers;
		IntervalMatrix<Number>		  E			  = taylor_.error;
		size_t						  taylorTerms = options.taylor_terms();
		double						  r			  = options.time_step();
		Number*						  factors	  = options.factor();
		// initialize Asum

		Zonotope<Number> Asum = Vdyn * r;
		for (size_t i = 0; i < taylorTerms; i++)
		{
			// compute powers
			Matrix_t<Number> temp	 = factors[i + 1] * Apower[i];
			Zonotope<Number> ApowerV = Vdyn * temp;
			// compute sums
			Asum = Asum + ApowerV;
		}
		// get error due to finite Taylor series
		Zonotope<Number> F;
		// if(isa(Vdyn,'zonoBundle')){
		//     F=E*Vdyn.Z{1}*r;
		// }else{
		F = Vdyn * E * r;
		// }
		// Compute error solution (dyn. + stat.)

		// cout << errorStat << endl;
		// cout << Asum << endl;
		// cout << F << endl;

		polyZonotope<Number> Rerror = errorStat + Asum + F;
		// polyZonotope<Number> Rerror;
		return Rerror;
	}

	template <typename Number>
	polyZonotope<Number> LinearSys<Number>::deltaReach(polyReachOptions<Number>& options, polyZonotope<Number> Rinit)
	{
		// compute delta reachable set
		Matrix_t<Number>	   eAt		 = taylor_.eAt;
		IntervalMatrix<Number> F		 = taylor_.F;
		Zonotope<Number>	   RV		 = taylor_.RV;
		Zonotope<Number>	   inputCorr = taylor_.inputCorr;
		Zonotope<Number>	   Rtrans	 = taylor_.Rtrans;

		// first time step homogensous solution
		size_t				 dim		   = F.rows();
		Matrix_t<Number>	 tempeat	   = (eAt - Eigen::MatrixXd::Identity(dim, dim));
		polyZonotope<Number> Rhom_tp_delta = Rinit * tempeat + Rtrans;
		polyZonotope<Number> Rhom;
		// if(isa(Rint "zonotope")){
		// polyZonotope<Number> O(Eigen::MatrixXd::Zero(dim,1));
		// Rhom = enclose(O, Rhom_tp_delta) + F*Rinit + inputCorr;
		//}
		// else if isa(Rinit, "polyZonotope")
		// else{
		Matrix_t<Number> temp1 = Eigen::MatrixXd::Zero(dim, dim);

		polyZonotope<Number> O = Rhom_tp_delta * temp1;
		Rhom				   = O.enclose(Rhom_tp_delta) + Rinit.toZonotope() * F + inputCorr;
		//}

		// cout << Rhom << endl;

		// reduce solution
		Rhom.Reduce(options.intermediate_order());
		RV.Reduce(options.intermediate_order());

		// cout << Rhom << endl;
		// cout << RV << endl;

		// total solution

		// original compitation
		polyZonotope<Number> Rdelta = Rhom + RV;
		return Rdelta;
	}

	template <typename Number>
	void LinearSys<Number>::lin_print()
	{
		std::cout << "linearsys" << std::endl;
		std::cout << "taylor_type" << std::endl;
		std::cout << "timestep: " << taylor_.timeStep << std::endl;
		std::cout << "eAt" << std::endl;
		std::cout << taylor_.eAt << std::endl;
		std::cout << "F" << std::endl;
		std::cout << "inf" << std::endl;
		std::cout << taylor_.F.inf << std::endl;
		std::cout << "sup" << std::endl;
		std::cout << taylor_.F.sup << std::endl;
		std::cout << "RV" << std::endl;
		std::cout << taylor_.RV << std::endl;
		std::cout << "inputCorr" << std::endl;
		std::cout << taylor_.inputCorr << std::endl;
		std::cout << "Rtrans" << std::endl;
		std::cout << taylor_.Rtrans << std::endl;
		std::cout << "v" << std::endl;
		std::cout << taylor_.V << std::endl;
		std::cout << "error" << std::endl; /*taylor_.inputF.inf = Asum.inf+Einput.inf;
   taylor_.inputF.sup = Asum.sup+Einput.sup;*/
		std::cout << "inf" << std::endl;
		std::cout << taylor_.error.inf << std::endl;
		std::cout << "sup" << std::endl;
		std::cout << taylor_.error.sup << std::endl;
		std::cout << "inputF" << std::endl;
		std::cout << "inf" << std::endl;
		std::cout << taylor_.inputF.inf << std::endl;
		std::cout << "sup" << std::endl;
		std::cout << taylor_.inputF.sup << std::endl;
		std::cout << "eAtInt" << std::endl;
		std::cout << "inf" << std::endl;
		std::cout << taylor_.eAtInt.inf << std::endl;
		std::cout << "sup" << std::endl;
		std::cout << taylor_.eAtInt.sup << std::endl;
		std::cout << "Powers" << std::endl;
		for (auto i : taylor_.powers)
		{
			std::cout << i << std::endl;
		}
		std::cout << "A" << std::endl;
		std::cout << A_ << std::endl;
		std::cout << "B" << std::endl;
		std::cout << B_ << std::endl;
		std::cout << "c" << std::endl;
		std::cout << c_ << std::endl;
		std::cout << "C" << std::endl;
		std::cout << C_ << std::endl;
		std::cout << "D_" << std::endl;
		std::cout << D_ << std::endl;
		std::cout << "k" << std::endl;
		std::cout << k_ << std::endl;
		std::cout << "krylov_" << std::endl;
		std::cout << krylov_ << std::endl;
	}

} // namespace reachSolver