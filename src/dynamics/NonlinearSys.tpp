/*
 * @Author: your name
 * @Date: 2021-11-02 19:20:19
 * @LastEditTime: 2021-12-30 22:39:33
 * @LastEditors: Please set LastEditors
 * @Description: In User Settings Edit
 * @FilePath: /Solver/src/dynamics/NonlinearSys.cpp
 */
/**
 * @file   NonlinearSys.cpp
 * @brief  NonlinearSys class
 * @author Yunjun Bai
 * @date April 2021
 * @version 1.0
 *
 * Reference:
 * CORA ../contDynamics/nonlinearSys/all
 */

/**
 * @author Changyuan Zhao
 * @date  Nov 2021
 * @version 1.1
 *
 */

#include <dynamics/NonlinearSys.h>

using namespace std;

namespace reachSolver
{

	// template class NonlinearSys<double>;

	/*****************************************************************************
 *                                                                           *
 *                           Constructors and Destructors                    *
 *                                                                           *
 *****************************************************************************/
	// template <typename Number>
	// NonlinearSys<Number>::NonlinearSys(function_type fun_handle)
	//     :ContDynamics(new String("nonlinearSys"), 10, 1, 1)
	//     :mFile_(fun_handle){
	//     std::vector<size_t> tmp_obj = numberOfInputs(fun_handle, 2);
	//     num_states_ = tmp_obj[0];
	//     num_inputs_ = tmp_obj[1];
	//     assert(fun_handle!=NULL);
	// }

	// template <typename Number>
	// NonlinearSys<Number>::NonlinearSys(std::string name, function_type
	// fun_handle):
	//     mFile_(fun_handle),
	//     name_(name){
	//     std::vector<size_t> tmp_obj = numberOfInputs(fun_handle, 2);
	//     num_states_ = tmp_obj[0];
	//     num_inputs_ = tmp_obj[1];
	//     assert(fun_handle!=NULL);
	// }
	/*
template <typename Number>
NonlinearSys<Number>::NonlinearSys()
{

}

template <typename Number>
NonlinearSys<Number>::NonlinearSys(size_t dimension)
{}
*/
	template <typename Number>
	NonlinearSys<Number>::NonlinearSys(capd::IMap fun_handle, size_t num_states, size_t num_inputs, size_t dimension)
		: ContDynamics<Number>(std::string("nonlinearSys"), num_states, num_inputs, dimension)
		, out_dimesion(dimension)
		, myFun(fun_handle)
		, state_num(num_states)
		, input_num(num_inputs)
	{}

	template <typename Number>
	NonlinearSys<Number>::NonlinearSys(
		std::string name, capd::IMap fun_handle, size_t num_states, size_t num_inputs, size_t dimension)
		: ContDynamics<Number>(name, num_states, num_inputs, dimension)
		, out_dimesion(dimension)
		, myFun(fun_handle)
		, state_num(num_states)
		, input_num(num_inputs)
	{}

	template <typename Number>
	NonlinearSys<Number>::~NonlinearSys()
	{}

	/*****************************************************************************
 *                                                                           *
 *                       Public Functions on Properties                      *
 *                                                                           *
 *****************************************************************************/

	// template <typename Number>
	// const typename NonlinearSys<Number>::function_type
	// NonlinearSys<Number>::mFile() const{
	//     return mFile_;
	// }

	// template <typename Number>
	// const typename NonlinearSys<Number>::function_type
	// NonlinearSys<Number>::jacobian() const{
	//     return jacobian_;
	// }

	// template <typename Number>
	// const typename NonlinearSys<Number>::function_type
	// NonlinearSys<Number>::hessian() const{
	//     return hessian_;
	// }

	// template <typename Number>
	// const typename NonlinearSys<Number>::function_type
	// NonlinearSys<Number>::thirdOrderTensor() const{
	//     return thirdOrderTensor_;
	// }

	// template <typename Number>
	// const typename NonlinearSys<Number>::function_type
	// NonlinearSys<Number>::tensors() const{
	//     return tensors_;
	// }

	/*****************************************************************************
 *                                                                           *
 *                       Compute Reachable Set                               *
 *                                                                           *
 *****************************************************************************/

	/*template <typename Number>
std::vector<ReachableSetElement<Number>>
NonlinearSys<Number>::UniformOutput(TimeInt<Number>& time_int){
    std::vector<ReachableSetElement<Number>>;

}

template <typename Number>
std::vector<ReachableSetElement<Number>>
NonlinearSys<Number>::UniformOutput(TimePoint<Number>& time_point){


}*/

	template <typename Number>
	ReachableSet<Number> NonlinearSys<Number>::createReachSetObject(TimeInt<Number>&   time_int,
																	TimePoint<Number>& time_point)
	{
		// ReachableSet<Number> R(time_point,time_int);
		ReachableSet<Number> R;
		R.set_Rtime_point(time_point);
		R.set_Rtime_interval(time_int);

		return R;
	}

	template <typename Number>
	std::vector<ReachableSet<Number>> NonlinearSys<Number>::createReachSetObjectSplit(TimeInt<Number>&   time_int, 
																   TimePoint<Number>& time_point)
	{
		vector<ReachableSet<Number>> res;
		
		//splitted sets -> bring to correct format
		vector<int> ind(time_point.rs(1).size());
		iota(ind.begin(), ind.end(), 0);
		vector<int> parent(time_point.rs(1).size(), 0);

		vector<TimeInt<Number>> timeInt_;
		vector<TimePoint<Number>> timePoint_;

		//loop over all time steps
		for(int i = 0; i < time_point.set().size(); i++){

			//modify index vector
			vector<int> ind_(time_point.rs(i).size(), 0);
			int maxInd = *max_element(ind.begin(), ind.end());

			for(int j = 0; j < time_point.rs(i).size(); j++){
				if(time_point.set()[i][j].parent() > -1 && i > 0){
					ind_[j] = maxInd + 1;
					parent.push_back(ind[time_point.set()[i][j].parent()]);
					maxInd = maxInd + 1;
				}else{
					ind_[j] = ind[time_point.set()[i][j].prev()];
				}
			}
			ind = ind_;
			
			// cout << i << endl;
			// cout << ind << endl;

			//copy entries
			for(int j = 0; j < time_point.rs(i).size(); j++){

				vector<ReachableSetElement<Number>> IntSet;
				ReachableSetElement<Number> ElementInt;
				
				vector<ReachableSetElement<Number>> PointSet;
				ReachableSetElement<Number> ElementPoint;

				if(ind[j] < timeInt_.size()){
					//timeInt_
					Zonotope<Number> Zono = time_int.set()[i][j].rs();
					ElementInt.set_rs(Zono);
					IntSet.push_back(ElementInt);
					timeInt_[ind[j]].addSet(IntSet);

					timeInt_[ind[j]].addTime(time_int.time(i));

					//timePoint_
					Zono = time_point.set()[i][j].rs();
					ElementPoint.set_rs(Zono);
					PointSet.push_back(ElementPoint);
					timePoint_[ind[j]].addSet(PointSet);

					timePoint_[ind[j]].addTime(time_point.time(i));
				}else{
					
					TimeInt<Number> timeIntTemp;
					TimePoint<Number> timePointTemp;

					//timeInt_
					Zonotope<Number> Zono = time_int.set()[i][j].rs();
					ElementInt.set_rs(Zono);
					IntSet.push_back(ElementInt);
					timeIntTemp.addSet(IntSet);

					timeIntTemp.addTime(time_int.time(i));

					//timePoint_
					Zono = time_point.set()[i][j].rs();
					ElementPoint.set_rs(Zono);
					PointSet.push_back(ElementPoint);
					timePointTemp.addSet(PointSet);

					timePointTemp.addTime(time_point.time(i));

					timeInt_.push_back(timeIntTemp);
					timePoint_.push_back(timePointTemp);
				}
			}
		}

		//generate reachSet object
		for(int i = 0; i < timeInt_.size(); i++){

			ReachableSet<Number> R_ = createReachSetObject(timeInt_[i], timePoint_[i]);
			if(parent[i] > -1){
				R_.set_parent_rs(parent[i]);
			}
			res.push_back(R_);
		}

		return res;
	}

	template <typename Number>
	polyReachableSet<Number> NonlinearSys<Number>::createReachSetObject(polyTimeInt<Number>&   time_int,
																		polyTimePoint<Number>& time_point)
	{
		// ReachableSet<Number> R(time_point,time_int);
		polyReachableSet<Number> R;
		R.set_Rtime_point(time_point);
		R.set_Rtime_interval(time_int);

		return R;
	}

	template <typename Number>
	void NonlinearSys<Number>::checkOptionsReach(ReachOptions<Number>& options, int hyb)
	{
		noinput = false;
		// if there is no input
		if (this->num_inputs() == 0)
		{
			std::cout << "no input" << std::endl;
			noinput = true;
			this->set_num_inputs(1);
			input_num = 1;
			options.set_U(Zonotope<Number>(1));
			options.set_uTrans_nonlin(Eigen::MatrixXd::Zero(1, 1));
			// options.set_uTrans_nonlin(Eigen::MatrixXd::Zero(1,1));
		}
		/*capd::IMap f(_f, this->dim() + this->num_inputs() , this->num_outputs(), 0,
  5); myFun = f;*/
	}

	template <typename Number>
	void NonlinearSys<Number>::checkOptionsReach(polyReachOptions<Number>& options, int hyb)
	{
		noinput = false;
		// if there is no input
		if (this->num_inputs() == 0)
		{
			std::cout << "no input" << std::endl;
			noinput = true;
			this->set_num_inputs(1);
			input_num = 1;
			options.set_U(Zonotope<Number>(1));
			options.set_uTrans_nonlin(Eigen::MatrixXd::Zero(1, 1));
			// options.set_uTrans_nonlin(Eigen::MatrixXd::Zero(1,1));
		}
		/*capd::IMap f(_f, this->dim() + this->num_inputs() , this->num_outputs(), 0,
  5); myFun = f;*/
	}

	// template <typename Number>
	// void NonlinearSys<Number>::derivatives(ReachOptions<Number>& options){

	// }

	// typedef Vector_t<double> (*mFile_type)(Vector_t<double> vector1,
	// Vector_t<double> vector2); typedef double (**mFile_f_type)(Vector_t<double>
	// vector1, Vector_t<double> vector2);

	/*
template <typename Number>
void NonlinearSys<Number>::jacobian_(Vector_t<Number> vector1, Vector_t<Number>
vector2, Matrix_t<Number>& matrix1, Matrix_t<Number>& matrix2){
// void jacobian(mFile_type mFile_, Vector_t<double> vector1, Vector_t<double>
vector2, Matrix_t<double> matrix1, Matrix_t<double> matrix2){

    autodiff::VectorXreal F;       // the output vector F = f(x, p, q) evaluated
together with Jacobian below

    matrix1  = autodiff::jacobian(this->mFile_, autodiff::wrt(vector1),
autodiff::at(vector1, vector2), F);       // evaluate the function and the
Jacobian matrix Jx = dF/dx matrix2  = autodiff::jacobian(this->mFile_,
autodiff::wrt(vector2), autodiff::at(vector1, vector2), F);       // evaluate
the function and the Jacobian matrix Jx = dF/du

    // std::cout << "F = \n" << F << std::endl;     // print the evaluated
output vector F
    // std::cout << "Jx = \n" << Jx << std::endl;   // print the evaluated
Jacobian matrix dF/dx
    // std::cout << "Jp = \n" << Jp << std::endl;   // print the evaluated
Jacobian matrix dF/dp
    // std::cout << "Jq = \n" << Jq << std::endl;   // print the evaluated
Jacobian matrix dF/dq
    // std::cout << "Jqpx = \n" << Jqpx << std::endl; // print the evaluated
Jacobian matrix [dF/dq, dF/dp, dF/dx]
}
*/
	template <typename Number>
	Vector_t<Number> NonlinearSys<Number>::mFile_(Vector_t<Number> vector1, Vector_t<Number> vector2)
	{
		int dimOut = this->num_outputs();
		int dimIn  = this->dim() + this->num_inputs();

		Vector_t<Number> result(dimOut);

		capd::IVector x(dimIn);
		// capd::Interval v(dimIn);
		for (auto i = 0; i < this->dim(); i++)
		{
			x[i] = vector1(i, 0);
		}
		for (auto i = 0; i < this->num_inputs(); i++)
		{
			x[i + this->dim()] = vector2(i, 0);
		}

		// capd::IVector x(dimIn,v);
		capd::IVector y = myFun(x);

		for (auto i = 0; i < dimOut; i++)
		{
			result(i, 0) = y[i].leftBound();
		}
		// delete []v;
		return result;
	}

	template <typename Number>
	void NonlinearSys<Number>::jacobian_(IntervalMatrix<Number>	 im1,
										 IntervalMatrix<Number>	 im2,
										 IntervalMatrix<Number>& matrix1,
										 IntervalMatrix<Number>& matrix2)
	{
		// void jacobian(mFile_type mFile_, Vector_t<double> vector1, Vector_t<double>
		// vector2, Matrix_t<double> matrix1, Matrix_t<double> matrix2){
		int dimOut = this->num_outputs();
		int dimIn  = this->dim() + this->num_inputs();

		// capd::Interval *v = new capd::Interval(dimIn);
		capd::IVector x(dimIn);
		for (auto i = 0; i < this->dim(); i++)
		{
			x[i] = im1(i, 0);
		}
		for (auto i = 0; i < this->num_inputs(); i++)
		{
			x[i + this->dim()] = im2(i, 0);
		}

		capd::IMatrix Df(dimOut, dimIn);
		// capd::IVector x(dimIn);
		capd::IVector y = myFun(x, Df);

		matrix1.resize(dimOut, this->dim());
		matrix2.resize(dimOut, this->num_inputs());
		for (auto i = 0; i < dimOut; i++)
		{
			for (auto j = 0; j < this->dim(); j++)
			{
				matrix1(i, j) = Df[i][j];
			}
			for (auto j = this->dim(); j < dimIn; j++)
			{
				matrix2(i, j - this->dim()) = Df[i][j];
			}
		}
	}

	template <typename Number>
	void NonlinearSys<Number>::jacobian_(Vector_t<Number>  vector1,
										 Vector_t<Number>  vector2,
										 Matrix_t<Number>& A,
										 Matrix_t<Number>& B)
	{
		int dimOut = this->num_outputs();
		int dimIn  = this->dim() + this->num_inputs();

		// capd::Interval *v = new capd::Interval(dimIn);
		capd::IVector x(dimIn);
		for (auto i = 0; i < this->dim(); i++)
		{
			x[i] = vector1(i, 0);
		}
		for (auto i = 0; i < this->num_inputs(); i++)
		{
			x[i + this->dim()] = vector2(i, 0);
		}

		capd::IMatrix Df(dimOut, dimIn);
		// capd::IVector x(dimIn, v);
		capd::IVector y = myFun(x, Df);

		A.resize(dimOut, this->dim());
		B.resize(dimOut, this->num_inputs());
		for (auto i = 0; i < dimOut; i++)
		{
			for (auto j = 0; j < this->dim(); j++)
			{
				A(i, j) = Df[i][j].leftBound();
			}
			for (auto j = this->dim(); j < dimIn; j++)
			{
				B(i, j - this->dim()) = Df[i][j].leftBound();
			}
		}
	}

	/*
template <typename Number>
std::vector<IntervalMatrix> NonlinearSys<Number>::hessian_(IntervalMatrix im1,
IntervalMatrix im2){
// std::vector<IntervalMatrix> hessian(mFile_f_type mFile_f_, IntervalMatrix
im1, IntervalMatrix im2){ typedef autodiff::dual2nd
(*func_p)(autodiff::ArrayXdual2nd& vector1, autodiff::ArrayXdual2nd& vector2);
    func_p mFile_f[6];
    for(int i=0; i<6; i++){
        mFile_f[i] = (func_p)this->mFile_f_[i];
    }
    // using Eigen::MatrixXd;
    // using Eigen::VectorXd;

    // Seprate inf and sup
    Matrix_t<double>im1_inf = im1.inf;
    Matrix_t<double>im1_sup = im1.sup;
    Matrix_t<double>im2_inf = im2.inf;
    Matrix_t<double>im2_sup = im2.sup;

    autodiff::ArrayXdual2nd tmp_im1_inf, tmp_im1_sup, tmp_im2_inf, tmp_im2_sup;
    tmp_im1_inf = autodiff::ArrayXdual2nd(im1_inf.rows());
    tmp_im1_sup = autodiff::ArrayXdual2nd(im1_sup.rows());
    tmp_im2_inf = autodiff::ArrayXdual2nd(im2_inf.rows());
    tmp_im2_sup = autodiff::ArrayXdual2nd(im2_sup.rows());
    for(int i=0; i<im1_inf.rows(); i++){
        tmp_im1_inf[i] = im1_inf(i, 0);
        tmp_im1_sup[i] = im1_sup(i, 0);
    }
    for(int i=0; i<im2_inf.rows(); i++){
        tmp_im2_inf[i] = im2_inf(i, 0);
        tmp_im2_sup[i] = im2_sup(i, 0);
    }

    // MatrixXd Jdyn_comb_inf = jacobian(this->mfile_, wrt(im1_inf, im2_inf),
at(im1_inf, im2_inf), F);
    // MatrixXd Jdyn_comb_sup = jacobian(this->mfile_, wrt(im1_sup, im2_sup),
at(im1_sup, im2_sup), F);
    // assert(Jdyn_comb_inf.rows() == Jdyn_comb_sup.rows());
    std::vector<IntervalMatrix> result =
std::vector<IntervalMatrix>(im1_inf.rows()); IntervalMatrix tmp_result;
    autodiff::dual2nd u; // the output scalar u = f(x) evaluated together with
Hessian below autodiff::VectorXdual g; // gradient of f(x) evaluated together
with Hessian below for(int k = 0; k<im1_inf.rows(); k++){ tmp_result.inf =
autodiff::hessian(mFile_f[k], autodiff::wrt(tmp_im1_inf, tmp_im2_inf),
autodiff::at(tmp_im1_inf, tmp_im2_inf), u, g); tmp_result.sup =
autodiff::hessian(mFile_f[k], autodiff::wrt(tmp_im1_sup, tmp_im2_sup),
autodiff::at(tmp_im1_sup, tmp_im2_sup), u, g); result[k] = tmp_result;
    }
    return result;

}
*/

	template <typename Number>
	std::vector<IntervalMatrix<Number>> NonlinearSys<Number>::hessian_(IntervalMatrix<Number> im1,
																	   IntervalMatrix<Number> im2)
	{
		int dimOut = this->num_outputs();
		int dimIn  = this->dim() + this->num_inputs();

		// capd::Interval *v = new capd::Interval(dimIn);
		capd::IVector x(dimIn);
		for (auto i = 0; i < this->dim(); i++)
		{
			x[i].setLeftBound(im1(i, 0).leftBound());
			x[i].setRightBound(im1(i, 0).rightBound());
		}
		for (auto i = 0; i < this->num_inputs(); i++)
		{
			x[i + this->dim()].setLeftBound(im2(i, 0).leftBound());
			x[i + this->dim()].setRightBound(im2(i, 0).rightBound());
		}

		capd::IMatrix  Df(dimOut, dimIn);
		capd::IHessian Hf(dimOut, dimIn);
		// capd::IVector x(dimIn, v);
		capd::IVector y = myFun(x, Df, Hf);

		std::vector<IntervalMatrix<Number>> result(dimOut);

		for (auto i = 0; i < dimOut; i++)
		{
			IntervalMatrix<Number> temp(dimIn, dimIn);
			for (auto j = 0; j < dimIn; j++)
			{
				for (auto k = 0; k < dimIn; k++)
				{
					temp(j, k).setLeftBound(2*Hf(i, j, k).leftBound());
					temp(j, k).setRightBound(2*Hf(i, j, k).rightBound());
				}
			}
			result[i] = temp;
		}

		return result;
	}

	template <typename Number>
	std::vector<Matrix_t<Number>> NonlinearSys<Number>::hessian_(Vector_t<Number> im1, Vector_t<Number> im2)
	{	
		// clock_t start, end;
		// start = clock();
		int dimOut = this->num_outputs();
		int dimIn  = this->dim() + this->num_inputs();

		// capd::Interval *v = new capd::Interval(dimIn);
		capd::IVector x(dimIn);
		for (auto i = 0; i < this->dim(); i++)
		{
			x[i] = im1(i);
		}
		for (auto i = 0; i < this->num_inputs(); i++)
		{
			x[i + this->dim()] = im2(i);
		}

		capd::IMatrix  Df(dimOut, dimIn);
		capd::IHessian Hf(dimOut, dimIn);
		// capd::IVector x(dimIn, v);
		capd::IVector y = myFun(x, Df, Hf);

		std::vector<Matrix_t<Number>> result(dimOut);

		for (auto i = 0; i < dimOut; i++)
		{
			Matrix_t<Number> temp(dimIn, dimIn);
			for (auto j = 0; j < dimIn; j++)
			{
				for (auto k = 0; k < dimIn; k++)
				{
					if(j == k){
						temp(j, k) = 2*Hf(i, j, k).leftBound();
					}else{
						temp(j, k) = Hf(i, j, k).leftBound();
					}

				}
			}
			result[i] = temp;
		}
	
	// end = clock();
	// cout << "hessian_ time cost: " << (double) (end - start) / CLOCKS_PER_SEC << endl;
		return result;
	}

	template <typename Number>
	Matrix_t<IntervalMatrix<Number>> NonlinearSys<Number>::thirdOrderTensor(IntervalMatrix<Number>		   im1,
																			IntervalMatrix<Number>		   im2,
																			Vector_t<std::vector<size_t>>& ind)
	{
		// clock_t start, end;
		// start = clock();

		int dimIn  = this->dim() + this->num_inputs();
		int dimOut = this->num_outputs();
		// Interval_t<Number> *v = new Interval_t<Number>(dimIn);
		capd::IVector x(dimIn);

		for (auto i = 0; i < this->dim(); i++)
		{
			x[i].setLeftBound(im1(i, 0).leftBound());
			x[i].setRightBound(im1(i, 0).rightBound());
		}
		for (auto i = 0; i < this->num_inputs(); i++)
		{
			x[i + this->dim()].setLeftBound(im2(i, 0).leftBound());
			x[i + this->dim()].setRightBound(im2(i, 0).rightBound());
		}
		// capd::IMap f(_f, dimIn, dimOut, 0, 3);
		Matrix_t<IntervalMatrix<Number>> result(dimOut, dimIn);
		// capd::IVector x(dimIn, v);
		capd::IJet jet(dimOut, dimIn, 3);
		// std::cout<<x<<std::endl;

		// cout << x << endl;

		myFun(x, jet);
		// std::cout<<jet.toString()<<std::endl;
		// std::cout<<jet(2,2,2,2)<<std::endl;
		// std::cout<<"before loop"<<std::endl;
		for (auto i = 0; i < dimOut; i++)
		{
			for (auto m = 0; m < dimIn; m++)
			{
				IntervalMatrix<Number> temp(dimIn, dimIn);
				result(i, m) = temp;
			}
		}

		int factor = 1;
		for (auto i = 0; i < dimOut; i++)
		{	
			for (auto j = 0; j < dimIn; j++)
			{
				for (auto k = j; k < dimIn; k++)
				{
					for (auto m = k; m < dimIn; m++)
					{
						if( m != k && k != j){
							factor = 1;
						}else if(m == k && k == j){
							factor = 6;
						}else{
							factor = 2;
						}
						result(i, j)(k, m).setLeftBound(factor*jet(i, j, k, m).leftBound());
						result(i, j)(k, m).setRightBound(factor*jet(i, j, k, m).rightBound());
						// if(abs(result(i,j)(k,m).rightBound())<1e-300){
						//     result(i,j)(k,m).setRightBound(0.);
						// }
						// if(abs(result(i,j)(k,m).leftBound())<1e-300){
						//     result(i,j)(k,m).setLeftBound(0.);
						// }
						result(i, k)(j, m) = result(i, j)(k, m);
						result(i, j)(m, k) = result(i, j)(k, m);
						result(i, k)(m, j) = result(i, j)(k, m);
						result(i, m)(j, k) = result(i, j)(k, m);
						result(i, m)(k, j) = result(i, j)(k, m);
					}
				}
			}
		}
		// std::cout<<"end loop"<<std::endl;
		// Vector_t<std::vector<Number>> &ind;
		ind.resize(dimOut);
		for (auto i = 0; i < dimOut; i++)
		{
			std::vector<size_t> temp;
			for (auto j = 0; j < dimIn; j++)
			{
				int flag = 0;
				for (auto m = 0; m < result(i, j).rows() && flag == 0; m++)
				{
					for (auto n = 0; n < result(i, j).cols(); n++)
					{
						if (result(i, j)(m, n).leftBound() != 0. || result(i, j)(m, n).rightBound() != 0.)
						{
							flag = 1;
							break;
						}
					}
				}
				if (flag == 1)
				{
					temp.push_back(j);
				}
			}
			ind(i) = temp;
		}

		// end = clock();
		// cout << "thirdOrderTensor time cost: " << (double) (end - start) / CLOCKS_PER_SEC << endl;
		return result;
	}

	// template <typename Number>
	// int NonlinearSys<Number>::reach(ReachOptions<Number>& options,
	// ReachSpecification& spec, ReachableSet<Number> & R){

	template <typename Number>
	int NonlinearSys<Number>::reach(ReachOptions<Number>& options, ReachableSet<Number>& R)
	{
		int res = 1;
		// options preprocessing
		checkOptionsReach(options, 0);

		// compute symbolic derivatives
		// derivatives(options);

		// obtain factors for initial state and input solution time step
		double	r	   = options.time_step();
		Number* factor = new double[options.taylor_terms() + 1];
		for (size_t i = 1; i <= options.taylor_terms() + 1; i++)
		{
			factor[i - 1] = std::pow(r, i) / std::tgamma(i + 1);
		}
		options.set_factor(factor);
		options.set_t(options.tStart());

		// if a trajectory should be tracked
		if (options.uTrans_Vec().size() != 0)
		{
			options.set_uTrans_nonlin(options.uTrans_Vec().col(0));
		}

		// time period
		int					tmp_t_count = (options.tFinal() - options.tStart()) / options.time_step();
		std::vector<Number> tVec		= std::vector<Number>(tmp_t_count + 1);
		for (size_t i = 0; i <= tmp_t_count; i++)
		{
			tVec[i] = options.tStart() + options.time_step() * i;
		}

		// initialize cell-arrays that store the reachable set
		TimeInt<Number>	  time_int	 = TimeInt<Number>(tmp_t_count);
		TimePoint<Number> time_point = TimePoint<Number>(tmp_t_count);

		// initialize reachable set computations
		ReachableSetElement<Number> Rinit_element = ReachableSetElement<Number>(options.R0(), options.max_error());

		std::vector<ReachableSetElement<Number>> Rinit = {Rinit_element};
		ReachableSet<Number>					 Rnext = ReachableSet<Number>();

		Rinit[0].set_error(0 * Vector_t<Number>(this->dim()));

		verboseLog(1, options.t(), options);
		try
		{
			Rnext = initReach(Rinit, options);
		}
		catch (SetExplosionException e1)
		{
			std::cout << e1.what() << std::endl;
			return 0;
		}
		std::cout << "after first initReach" << std::endl;

		//cout << "Rnext.time_point().size()" << Rnext.time_point().size() << endl;
		// loop over all reachability steps
		Number t;
		for (size_t i = 1; i < tVec.size() - 1; i++)
		{
			time_int.set_rs(i - 1, Rnext.time_interval());
			Interval_t<Number> tmp_interval;
			tmp_interval.setLeftBound(tVec[i - 1]);
			tmp_interval.setRightBound(tVec[i]);
			time_int.set_time(i - 1, tmp_interval);
			time_point.set_rs(i - 1, Rnext.time_point());
			time_point.set_time(i - 1, tVec[i]);

			// check specificaiton
			// if(!spec.empty()){
			//     if(!spec.check(Rnext.ti)){
			//         R = createReachSetObject(time_int, time_point);
			//         return 0;
			//     }
			// }

			// increment time and set counter
			t = tVec[i];
			options.set_t(t);
			if (options.verbose() == 1)
			{
				std::cout << t << std::endl;
			}

			// if a trajectory should be tracked
			if (options.uTrans_Vec().size() != 0)
			{
				options.set_uTrans_nonlin(options.uTrans_Vec().col(i));
			}

			// compute next reachable set
			try
			{
				Rnext = post(Rnext, options);
			}
			catch (SetExplosionException e1)
			{
				R = createReachSetObject(time_int, time_point);
				std::cout << e1.what() << std::endl;
				return 0;
			}
		}
		// check specificaiton
		// if(!spec.empty()){
		//     if(!spec.check(Rnext.ti)){
		//         res = 0;
		//     }
		// }
		time_int.set_rs(tmp_t_count - 1, Rnext.time_interval());
		Interval_t<Number> tmp_interval;
		tmp_interval.setLeftBound(tVec[tmp_t_count - 1]);
		tmp_interval.setRightBound(tVec[tmp_t_count]);
		time_int.set_time(tmp_t_count - 1, tmp_interval);
		time_point.set_rs(tmp_t_count - 1, Rnext.time_point());
		time_point.set_time(tmp_t_count - 1, tVec[tmp_t_count]);

		// cout << time_point.set().size();
		// cout << time_point.rs(2).size();

		R = createReachSetObject(time_int, time_point);
		delete []factor;
		return res;

	} // NonlinearSys<Number>::reach

	template <typename Number>
	int NonlinearSys<Number>::reachSplit(ReachOptions<Number>& options, vector<ReachableSet<Number>>& R)
	{
		int res = 1;
		// options preprocessing
		checkOptionsReach(options, 0);

		// compute symbolic derivatives
		// derivatives(options);

		// obtain factors for initial state and input solution time step
		double	r	   = options.time_step();
		Number* factor = new double[options.taylor_terms() + 1];
		for (size_t i = 1; i <= options.taylor_terms() + 1; i++)
		{
			factor[i - 1] = std::pow(r, i) / std::tgamma(i + 1);
		}
		options.set_factor(factor);
		options.set_t(options.tStart());

		// if a trajectory should be tracked
		if (options.uTrans_Vec().size() != 0)
		{
			options.set_uTrans_nonlin(options.uTrans_Vec().col(0));
		}

		// time period
		int					tmp_t_count = (options.tFinal() - options.tStart()) / options.time_step();
		std::vector<Number> tVec		= std::vector<Number>(tmp_t_count + 1);
		for (size_t i = 0; i <= tmp_t_count; i++)
		{
			tVec[i] = options.tStart() + options.time_step() * i;
		}

		// initialize cell-arrays that store the reachable set
		TimeInt<Number>	  time_int	 = TimeInt<Number>(tmp_t_count);
		TimePoint<Number> time_point = TimePoint<Number>(tmp_t_count);

		// initialize reachable set computations
		ReachableSetElement<Number> Rinit_element = ReachableSetElement<Number>(options.R0(), options.max_error());

		std::vector<ReachableSetElement<Number>> Rinit = {Rinit_element};
		ReachableSet<Number>					 Rnext = ReachableSet<Number>();

		Rinit[0].set_error(0 * Vector_t<Number>(this->dim()));

		verboseLog(1, options.t(), options);
		try
		{
			Rnext = initReach(Rinit, options);
		}
		catch (SetExplosionException e1)
		{
			std::cout << e1.what() << std::endl;
			return 0;
		}
		std::cout << "after first initReach" << std::endl;

		//cout << "Rnext.time_point().size()" << Rnext.time_point().size() << endl;
		// loop over all reachability steps
		Number t;
		for (size_t i = 1; i < tVec.size() - 1; i++)
		{
			time_int.set_rs(i - 1, Rnext.time_interval());
			Interval_t<Number> tmp_interval;
			tmp_interval.setLeftBound(tVec[i - 1]);
			tmp_interval.setRightBound(tVec[i]);
			time_int.set_time(i - 1, tmp_interval);
			time_point.set_rs(i - 1, Rnext.time_point());
			time_point.set_time(i - 1, tVec[i]);

			// check specificaiton
			// if(!spec.empty()){
			//     if(!spec.check(Rnext.ti)){
			//         R = createReachSetObject(time_int, time_point);
			//         return 0;
			//     }
			// }

			// increment time and set counter
			t = tVec[i];
			options.set_t(t);
			if (options.verbose() == 1)
			{
				std::cout << t << std::endl;
			}

			// if a trajectory should be tracked
			if (options.uTrans_Vec().size() != 0)
			{
				options.set_uTrans_nonlin(options.uTrans_Vec().col(i));
			}

			// compute next reachable set
			try
			{
				Rnext = post(Rnext, options);
			}
			catch (SetExplosionException e1)
			{
				R = createReachSetObjectSplit(time_int, time_point);
				std::cout << e1.what() << std::endl;
				return 0;
			}
		}
		// check specificaiton
		// if(!spec.empty()){
		//     if(!spec.check(Rnext.ti)){
		//         res = 0;
		//     }
		// }
		time_int.set_rs(tmp_t_count - 1, Rnext.time_interval());
		Interval_t<Number> tmp_interval;
		tmp_interval.setLeftBound(tVec[tmp_t_count - 1]);
		tmp_interval.setRightBound(tVec[tmp_t_count]);
		time_int.set_time(tmp_t_count - 1, tmp_interval);
		time_point.set_rs(tmp_t_count - 1, Rnext.time_point());
		time_point.set_time(tmp_t_count - 1, tVec[tmp_t_count]);

		// cout << time_point.set().size();
		// cout << time_point.rs(139).size();

		R = createReachSetObjectSplit(time_int, time_point);
		delete []factor;
		return res;

	} // NonlinearSys<Number>::reach

	template <typename Number>
	int NonlinearSys<Number>::reach(polyReachOptions<Number>& options, polyReachableSet<Number>& R)
	{
		int res = 1;
		// options preprocessing
		if(isFullDim(options.R0())){
			options.set_volApproxMethod("interval");
		}else{
			options.set_volApproxMethod("pca");
		}

		//边界算法用 pca 效果不好
		options.set_volApproxMethod("interval");

		checkOptionsReach(options, 0);

		// compute symbolic derivatives
		// derivatives(options);

		// obtain factors for initial state and input solution time step
		double	r	   = options.time_step();
		Number* factor = new double[options.taylor_terms() + 1];
		for (size_t i = 1; i <= options.taylor_terms() + 1; i++)
		{
			factor[i - 1] = std::pow(r, i) / std::tgamma(i + 1);
		}
		options.set_factor(factor);
		options.set_t(options.tStart());

		// if a trajectory should be tracked
		if (options.uTrans_Vec().size() != 0)
		{
			options.set_uTrans_nonlin(options.uTrans_Vec().col(0));
		}

		// time period
		int					tmp_t_count = (options.tFinal() - options.tStart()) / options.time_step();
		std::vector<Number> tVec		= std::vector<Number>(tmp_t_count + 1);
		for (size_t i = 0; i <= tmp_t_count; i++)
		{
			tVec[i] = options.tStart() + options.time_step() * i;
		}

		// initialize cell-arrays that store the reachable set
		polyTimeInt<Number>	  time_int	 = polyTimeInt<Number>(tmp_t_count);
		polyTimePoint<Number> time_point = polyTimePoint<Number>(tmp_t_count);

		// initialize reachable set computations
		polyReachableSetElement<Number> Rinit_element
			= polyReachableSetElement<Number>(options.R0(), options.max_error());

		std::vector<polyReachableSetElement<Number>> Rinit = {Rinit_element};
		polyReachableSet<Number>					 Rnext = polyReachableSet<Number>();

		Rinit[0].set_error(0 * Vector_t<Number>(this->dim()));

		verboseLog(1, options.t(), options);
		try
		{
			Rnext = initReach(Rinit, options);
		}
		catch (SetExplosionException e1)
		{
			std::cout << e1.what() << std::endl;
			return 0;
		}
		// std::cout << "after first initReach" << std::endl;
		// loop over all reachability steps
		Number t;
		for (size_t i = 1; i < tVec.size() - 1; i++)
		{
			time_int.set_rs(i - 1, Rnext.time_interval());
			Interval_t<Number> tmp_interval;
			tmp_interval.setLeftBound(tVec[i - 1]);
			tmp_interval.setRightBound(tVec[i]);
			time_int.set_time(i - 1, tmp_interval);
			time_point.set_rs(i - 1, Rnext.time_point());
			time_point.set_time(i - 1, tVec[i]);

			// check specificaiton
			// if(!spec.empty()){
			//     if(!spec.check(Rnext.ti)){
			//         R = createReachSetObject(time_int, time_point);
			//         return 0;
			//     }
			// }

			// increment time and set counter
			t = tVec[i];
			options.set_t(t);
			if (options.verbose() == 1)
			{
				std::cout << t << std::endl;
			}

			// if a trajectory should be tracked
			if (options.uTrans_Vec().size() != 0)
			{
				options.set_uTrans_nonlin(options.uTrans_Vec().col(i));
			}

			// cout << Rnext.time_point().size() << endl;
			// cout << Rnext.time_point()[0].rs() << endl;
			// cout << Rnext.time_point()[0].error() << endl;
			// cout << Rnext.time_point()[0].prev() << endl;
			// cout << Rnext.time_point()[0].parent()<< endl;

			// polyReachableSet<Number> Rnexttemp;
			// Rnexttemp = initReach(Rnext.time_point(), options);

			// compute next reachable set
			try
			{
				Rnext = post(Rnext, options);
				// std::cout << "after " << i << "th"
						//   << " initReach" << std::endl;
			}
			catch (SetExplosionException e1)
			{
				R = createReachSetObject(time_int, time_point);
				std::cout << e1.what() << std::endl;
				return 0;
			}

			// cout << Rnext.time_point().size() << endl;
			// cout << Rnext.time_interval().size() << endl;

		}
		// check specificaiton
		// if(!spec.empty()){
		//     if(!spec.check(Rnext.ti)){
		//         res = 0;
		//     }
		// }

		time_int.set_rs(tmp_t_count - 1, Rnext.time_interval());
		Interval_t<Number> tmp_interval;
		tmp_interval.setLeftBound(tVec[tmp_t_count - 1]);
		tmp_interval.setRightBound(tVec[tmp_t_count]);
		time_int.set_time(tmp_t_count - 1, tmp_interval);
		time_point.set_rs(tmp_t_count - 1, Rnext.time_point());
		time_point.set_time(tmp_t_count - 1, tVec[tmp_t_count]);

		R = createReachSetObject(time_int, time_point);
		delete[] factor;
		return res;

	} // NonlinearSys<Number>::reach

	template <typename Number>
	ReachableSet<Number> NonlinearSys<Number>::initReach_linRem(std::vector<ReachableSetElement<Number>>& Rinit,
																ReachOptions<Number>&					  options)
	{
		return ReachableSet<Number>();
	}

	// template <typename Number>
	// std::vector<Zonotope<Number>> splitOneDim(Zonotope<Number> Zbundle,
	// Matrix_t<double> leftLimit, Matrix_t<double> rightLimit, int dim){
	//     // split limits for a given dimension
	//     Matrix_t<double> leftLimitMod = leftLimit;
	//     leftLimitMod[dim] = 0.5*(leftLimit[dim]+rightLimit[dim]);
	//     Matrix_t<double>  rightLimitMod = rightLimit;
	//     rightLimitMod[dim] = 0.5*(leftLimit[dim]+rightLimit[dim]);

	//     // construct zonotopes which are the left and right boxes
	//     IntervalMatrix IM_Zleft, IM_Zright;
	//     IM_Zleft.inf = leftLimit; IM_Zleft.sup = rightLimitMod;
	//     IM_Zright.inf = leftLimitMod; IM_Zright.sup = rightLimit;
	//     Zonotope<Number> Zleft = Zonotope<Number>(IM_Zleft);
	//     Zonotope<Number> Zright = Zonotope<Number>(IM_Zright);

	//     // generate splitted zonotope bundles
	//     std::vector<Zonotope<Number>> Zsplit = new
	//     std::vector<Zonotope<Number>>(2); Zsplit[0] = Zbundle & Zleft; Zsplit[1]
	//     = Zbundle & Zright; return Zsplit;
	// }

	// template <typename Number>
	// std::vector<Zonotope<Number>> NonlinearSys<Number>::split(Zonotope<Number>
	// Zbundle, int number){
	//     std::vector<Zonotope<Number>> Zsplit;
	//     // split given dimension
	//     if(number > 0){
	//         int N = number;
	//         // obtain enclosing interval hull
	//         IntervalMatrix IH = Zbundle.interval();
	//         // obtain limits
	//         Matrix_t<double> leftLimit = IH.inf;
	//         Matrix_t<double> rightLimit = IH.sup;
	//         // split one dimension
	//         Zsplit = splitOneDim(Zbundle,leftLimit,rightLimit,N);
	//     }else{
	//         Zsplit = new std::vector<Zonotope<Number>>();
	//     }
	//     return Zsplit;
	// }

	template <typename Number>
	ReachableSet<Number> NonlinearSys<Number>::initReach(std::vector<ReachableSetElement<Number>> Rinit,
														 ReachOptions<Number>&					  options)
	{
		/*if(Rinit.size()==1){
      cout<<"adadad"<<endl;
      Rinit[0].set_error(Vector_t<Number>(Rinit[0].error().rows()));
  }*/
		// std::cout << "this is initReach" << std::endl;
		ReachableSet<Number> Rnext = ReachableSet<Number>();
		// compute reachable set using the options.alg = 'linRem' algorithm
		if (options.alg() == std::string("linRem"))
		{
			Rnext = initReach_linRem(Rinit, options);
		}

		// loop over all parallel sets
		int										 setCounter = 0;
		std::vector<ReachableSetElement<Number>> Rtp		= std::vector<ReachableSetElement<Number>>();
		std::vector<ReachableSetElement<Number>> Rti		= std::vector<ReachableSetElement<Number>>();
		std::vector<ReachableSetElement<Number>> R0			= std::vector<ReachableSetElement<Number>>();
		ReachableSetElement<Number>				 Rtemp_tp	= ReachableSetElement<Number>();
		ReachableSetElement<Number>				 Rtemp_ti	= ReachableSetElement<Number>();

		for (size_t i = 0; i < Rinit.size(); i++)
		{

			// cout << Rinit[i].rs() << endl;
			// cout << Rinit[i].error() << endl;

			// compute reachable set of abstraction
			int dimForSplit = linReach(options, Rinit[i], Rtemp_ti, Rtemp_tp);
			// std::cout<<"there"<<std::endl;
			// std::cout << "dimForSplit" << std::endl;
			// std::cout << dimForSplit << std::endl;
			// check if initial set has to be split
			if (dimForSplit == -1)
			{
				Rtp.push_back(Rtemp_tp);
				Rtp[setCounter].set_prev(i);
				Rti.push_back(Rtemp_ti);
				R0.push_back(Rinit[i]);
				setCounter = setCounter + 1;
			}
			else
			{
				std::cout << "split! ...number of parallel sets: " << Rinit.size() + 1 << std::endl;

				// split the initial set
				std::vector<Zonotope<Number>>			 Rtmp = Rinit[i].rs().split(dimForSplit);
				std::vector<ReachableSetElement<Number>> Rsplit(2);
				Rsplit[0].set_rs(Rtmp[0]);
				Rsplit[1].set_rs(Rtmp[1]);
				
				// cout << Rtmp << endl;

				// reset the linearization error
				Rsplit[0].set_error(Eigen::MatrixXd::Zero(options.max_error().rows(), 1));
				Rsplit[1].set_error(Eigen::MatrixXd::Zero(options.max_error().rows(), 1));

				// compute the reachable set for the splitted sets
				ReachableSet<Number> Rres;
				Rres = initReach(Rsplit, options);
				// std::cout << "setCounter" << std::endl;
				// std::cout << setCounter << std::endl;
				// std::cout << "Rtp.size()" << std::endl;
				// std::cout << Rtp.size() << std::endl;
				// std::cout << "Rres.time_point().size()" << std::endl;
				// std::cout << Rres.time_point().size() << std::endl;
				for (size_t j = 0; j < Rres.time_point().size(); j++)
				{
					// std::cout<<"1"<<std::endl;
					Rtp.push_back(Rres.time_point()[j]);
					// std::cout<<"2"<<std::endl;

					cout << "before " << Rtp[setCounter].parent() << endl;

					Rtp[setCounter].set_parent(i);

					cout << "after " << Rtp[setCounter].parent() << endl;

					// std::cout<<"3"<<std::endl;
					Rti.push_back(Rres.time_interval()[j]);
					R0.push_back(Rres.R0()[j]);
					setCounter = setCounter + 1;
				} // for j
			}	  // else
		}		  // for i

		Rnext.set_time_point(Rtp);
		Rnext.set_time_interval(Rti);
		Rnext.set_R0(R0);
		// std::cout << "end initReach" << std::endl;
		return Rnext;
	} // initReach

	template <typename Number>
	polyReachableSet<Number> NonlinearSys<Number>::initReach(std::vector<polyReachableSetElement<Number>> Rinit,
															 polyReachOptions<Number>&					  options)
	{
		/*if(Rinit.size()==1){
      cout<<"adadad"<<endl;
      Rinit[0].set_error(Vector_t<Number>(Rinit[0].error().rows()));
  }*/
		// std::cout << "this is initReach" << std::endl;
		polyReachableSet<Number> Rnext = polyReachableSet<Number>();
		// compute reachable set using the options.alg = 'linRem' algorithm

		// loop over all parallel sets
		int											 setCounter = 0;
		std::vector<polyReachableSetElement<Number>> Rtp		= std::vector<polyReachableSetElement<Number>>();
		std::vector<polyReachableSetElement<Number>> Rti		= std::vector<polyReachableSetElement<Number>>();
		std::vector<polyReachableSetElement<Number>> R0			= std::vector<polyReachableSetElement<Number>>();
		polyReachableSetElement<Number>				 Rtemp_tp	= polyReachableSetElement<Number>();
		polyReachableSetElement<Number>				 Rtemp_ti	= polyReachableSetElement<Number>();

		for (size_t i = 0; i < Rinit.size(); i++)
		{
			// compute reachable set of abstraction

			// cout << Rinit[i].rs() << endl;
			// cout << Rinit[i].error() << endl;
			// cout << Rtemp_ti.rs() << endl;			
			// cout << Rtemp_tp.rs() << endl;
			// cout << Rtemp_tp.error() << endl;

			int dimForSplit = linReach(options, Rinit[i], Rtemp_ti, Rtemp_tp);
			// std::cout<<"there"<<std::endl;
			// std::cout << "dimForSplit" << std::endl;
			// std::cout << dimForSplit << std::endl;

			// cout << Rtemp_ti.rs() << endl;			
			// cout << Rtemp_tp.rs() << endl;
			// cout << Rtemp_tp.error() << endl;

			// check if initial set has to be split
			if (dimForSplit == -1)
			{
				Rtp.push_back(Rtemp_tp);
				Rtp[setCounter].set_prev(i);
				Rti.push_back(Rtemp_ti);
				R0.push_back(Rinit[i]);
				setCounter = setCounter + 1;
			}
			else
			{
				std::cout << "split! ...number of parallel sets: " << Rinit.size() + 1 << std::endl;

			} // else
		}	  // for i

		Rnext.set_time_point(Rtp);
		Rnext.set_time_interval(Rti);
		Rnext.set_R0(R0);
		// std::cout << "end initReach" << std::endl;
		return Rnext;
	} // polyinitReach

	template <typename Number>
	ReachableSet<Number> NonlinearSys<Number>::post(ReachableSet<Number>& R, ReachOptions<Number>& options)
	{
		// In contrast to the linear system: the nonlinear system has to be constantly
		// initialized due to the linearization procedure
		ReachableSet<Number> Rnext;
		Rnext = initReach(R.time_point(), options);

		// reduce zonotopes
		for (size_t i = 0; i < Rnext.time_point().size(); i++)
		{
			if (!Rnext.time_point()[i].rs().Empty())
			{
				// clock_t start, end;
				// start = clock();
				std::vector<ReachableSetElement<Number>> timepoint1 = Rnext.time_point();
				Zonotope<Number>						 timepoint2 = timepoint1[i].rs();
				timepoint2.Reduce(options.zonotope_order());
				timepoint1[i].set_rs(timepoint2);
				Rnext.set_time_point(timepoint1);
				std::vector<ReachableSetElement<Number>> timeinterval1 = Rnext.time_interval();
				Zonotope<Number>						 timeinterval2 = timeinterval1[i].rs();
				timeinterval2.Reduce(options.zonotope_order());
				timeinterval1[i].set_rs(timeinterval2);
				Rnext.set_time_interval(timeinterval1);
				// end = clock();
				// cout<<"time cost: "<<(double)(end-start)/CLOCKS_PER_SEC<<endl;
				// Rnext.time_point()[i].rs().Reduce(options.zonotope_order());
				// Rnext.time_interval()[i].rs().Reduce(options.zonotope_order());
			}
		}

		// delete redundant reachable sets
		Rnext.deleteRedundantSets(R, options);
		return Rnext;
	}

	template <typename Number>
	polyReachableSet<Number> NonlinearSys<Number>::post(polyReachableSet<Number>& R, polyReachOptions<Number>& options)
	{	
		//potentially restructure the polynomial zonotope
		for(int i = 0; i < R.time_point().size(); i++){
			Number ratio = approxVolumeRatio(R.time_point()[i].rs(), options);
			cout << R.internalCount() << endl;
			cout << ratio << endl;

			if(ratio > options.maxPolyZonoRatio()){
				cout << "restructure !!!" << endl;
				polyZonotope<Number> temp = restructure(R.time_point()[i].rs(), options.maxDepGenOrder());
				R.set_time_point_rs(i,temp);
				// cout << R.time_point()[i].rs() << endl;
			}
			// cout << R.time_point()[i].rs() << endl;
			// cout << R.time_point()[i].error() << endl;
		}


		// In contrast to the linear system: the nonlinear system has to be constantly
		// initialized due to the linearization procedure
		polyReachableSet<Number> Rnext;

		// cout << R.time_point()[0].rs() <<endl;
		// cout << R.time_point()[0].error() <<endl;

		Rnext = initReach(R.time_point(), options);
		// std::cout << "before post" << std::endl;

		// cout << Rnext.time_point()[0].rs() <<endl;
		// cout << Rnext.time_point()[0].error() <<endl;
		// cout << Rnext.time_interval()[0].rs()<<endl;
		// cout << Rnext.R0()[0].rs() <<endl;
		// cout << Rnext.R0()[0].error() <<endl;

		// reduce zonotopes
		for (size_t i = 0; i < Rnext.time_point().size(); i++)
		{
			if (!Rnext.time_point()[i].rs().Empty())
			{
				std::vector<polyReachableSetElement<Number>> timepoint1 = Rnext.time_point();
				polyZonotope<Number>						 timepoint2 = timepoint1[i].rs();
				timepoint2.Reduce(options.zonotope_order());
				timepoint1[i].set_rs(timepoint2);
				Rnext.set_time_point(timepoint1);
				std::vector<polyReachableSetElement<Number>> timeinterval1 = Rnext.time_interval();
				polyZonotope<Number>						 timeinterval2 = timeinterval1[i].rs();
				timeinterval2.Reduce(options.zonotope_order());
				timeinterval1[i].set_rs(timeinterval2);
				Rnext.set_time_interval(timeinterval1);
			}
		}

		// delete redundant reachable sets
		Rnext.deleteRedundantSets(R, options);
		// std::cout << "end post" << std::endl;
		return Rnext;
	}

	template <typename Number>
	LinearSys<Number> NonlinearSys<Number>::linearize(ReachOptions<Number>& options,
													  Zonotope<Number>&		R,
													  ReachOptions<Number>& linOptions)
	{
		// linearization point p.u of the input is the center of the input set
		linerror_p_type p;
		p.u = options.uTrans_nonlin();
		// std::cout<<"(options.uTrans_nonlin())"<<std::endl;
		// std::cout<<(options.uTrans_nonlin())<<std::endl;
		// obtain linearization point
		//  if(options.haslinearizationPoint() == true){
		//      p.x = options.linearizationPoint();
		//  }else{
		//  linearization point p.x of the state is the center of the last reachable
		//  set R translated by 0.5*delta_t*f0

		// std::cout << "this is linearize" << std::endl;
		Vector_t<Number> tmp_vector1 = Vector_t<Number>(this->num_inputs());
		tmp_vector1 << p.u;

		// std::cout<<"(options.uTrans_nonlin())"<<std::endl;
		// std::cout<<(options.uTrans_nonlin())<<std::endl;

		Vector_t<Number> f0prev = mFile_(R.center(), tmp_vector1);

		// std::cout<<"(options.uTrans_nonlin())"<<std::endl;
		// std::cout<<(options.uTrans_nonlin())<<std::endl;

		try
		{
			p.x = R.center() + f0prev * 0.5 * options.time_step();
		}
		catch (const std::exception& e1)
		{
			p.x = R.center();
		}

		// }
		// substitute p into the system equation to obtain the constant input
		Vector_t<Number> f0 = mFile_(p.x, tmp_vector1);

		// substitute p into the Jacobian with respect to x and u to obtain the system
		// matrix A and the input matrix B
		Matrix_t<Number> A = Matrix_t<Number>(this->num_outputs(), this->dim());
		Matrix_t<Number> B = Matrix_t<Number>(this->num_outputs(), this->num_inputs());
		// std::cout << "before jacobian" << std::endl;

		jacobian_(p.x, tmp_vector1, A, B);

		// std::cout << "after jacobian" << std::endl;
		Matrix_t<Number> A_lin = A;
		Matrix_t<Number> B_lin = B;
		linOptions			   = options;

		// if (strcmp(options.alg(),'linRem') == 0){
		//     // in order to compute dA,dB, we use the reachability set computed for
		//     one step in initReach [dA,dB] =
		//     lin_error2dAB(options.Ronestep,options.U,obj.hessian,p); A =
		//     matZonotope(A,{dA}); B = matZonotope(B,{dB}); linSys =
		//     linParamSys(A,1,'constParam'); linOptions.compTimePoint = 1;
		// }else{
		// set up linearized system

		// tmp_B = 1 as input matrix encountered in uncertain inputs
		Matrix_t<Number> tmp_B	 = Matrix_t<Number>(1, 1);
		tmp_B(0, 0)				 = 1;
		LinearSys<Number> linSys = LinearSys<Number>(std::string("linSys"), A, tmp_B);
		// std::cout<<"after1 "<<std::endl;
		// B=1 as input matrix encountered in uncertain inputs
		// }

		// set up options for linearized system
		// std::cout<<"after3zzz"<<std::endl;
		linOptions.set_U((options.U() + options.uTrans_nonlin() - p.u) * B);
		// std::cout<<"after3"<<std::endl;
		Vector_t<Number> Ucenter = linOptions.U().center();
		linOptions.set_U(linOptions.U() - Ucenter);
		// std::cout<<"after4"<<std::endl;
		linOptions.set_uTrans_lin(Zonotope<Number>((f0 + Ucenter), Eigen::MatrixXd::Zero(f0.rows(), 1)));
		// linOptions.uTrans = zonotope(f0 + Ucenter);

		linOptions.set_originContained(0);
		// std::cout<<"after2"<<std::endl;
		// save constant input
		linerror_.f0 = f0;

		// save linearization point
		linerror_.p = p;

		// std::cout << "end linearize" << std::endl;

		return linSys;
	}

	template <typename Number>
	LinearSys<Number> NonlinearSys<Number>::linearize(polyReachOptions<Number>& options,
													  polyZonotope<Number>&		R,
													  polyReachOptions<Number>& linOptions)
	{
		// linearization point p.u of the input is the center of the input set
		linerror_p_type p;
		p.u = options.uTrans_nonlin();

		// linearization point p.x of the state is the center of the last reachable
		// set R translated by 0.5*delta_t*f0
		// std::cout << "this is linearize" << std::endl;
		Vector_t<Number> tmp_vector1 = Vector_t<Number>(this->num_inputs());
		tmp_vector1 << p.u;
		Vector_t<Number> f0prev = mFile_(R.center(), tmp_vector1);

		try
		{
			p.x = R.center() + f0prev * 0.5 * options.time_step();
		}
		catch (const std::exception& e1)
		{
			p.x = R.center();
		}
		// substitute p into the system equation to obtain the constant input
		Vector_t<Number> f0 = mFile_(p.x, tmp_vector1);

		// substitute p into the Jacobian with respect to x and u to obtain the system
		// matrix A and the input matrix B
		Matrix_t<Number> A = Matrix_t<Number>(this->num_outputs(), this->dim());
		Matrix_t<Number> B = Matrix_t<Number>(this->num_outputs(), this->num_inputs());

		jacobian_(p.x, tmp_vector1, A, B);

		Matrix_t<Number> A_lin = A;
		Matrix_t<Number> B_lin = B;
		linOptions			   = options;

		// if (strcmp(options.alg(),'linRem') == 0){
		//     // in order to compute dA,dB, we use the reachability set computed for
		//     one step in initReach [dA,dB] =
		//     lin_error2dAB(options.Ronestep,options.U,obj.hessian,p); A =
		//     matZonotope(A,{dA}); B = matZonotope(B,{dB}); linSys =
		//     linParamSys(A,1,'constParam'); linOptions.compTimePoint = 1;
		// }else{
		// set up linearized system

		// tmp_B = 1 as input matrix encountered in uncertain inputs
		Matrix_t<Number> tmp_B	 = Matrix_t<Number>(1, 1);
		tmp_B(0, 0)				 = 1;
		LinearSys<Number> linSys = LinearSys<Number>(std::string("linSys"), A, tmp_B);
		// std::cout<<"after1 "<<std::endl;
		// B=1 as input matrix encountered in uncertain inputs
		// }

		// set up options for linearized system
		linOptions.set_U((options.U() + options.uTrans_nonlin() - p.u) * B);
		Vector_t<Number> Ucenter = linOptions.U().center();
		linOptions.set_U(linOptions.U() - Ucenter);
		linOptions.set_uTrans_lin(Zonotope<Number>((f0 + Ucenter), Eigen::MatrixXd::Zero(f0.rows(), 1)));
		// linOptions.uTrans = zonotope(f0 + Ucenter);
		linOptions.set_originContained(0);

		// save constant input
		linerror_.f0 = f0;

		// save linearization point
		linerror_.p = p;

		// std::cout << "end linearize" << std::endl;
		return linSys;
	} // linearize poly

	template <typename Number>
	double NonlinearSys<Number>::linReach_linRem(LinearReachableSet<Number>& R,
												 Zonotope<Number>&			 Rinit,
												 Zonotope<Number>&			 Rdelta,
												 ReachOptions<Number>&		 options,
												 Zonotope<Number>&			 Rti,
												 Zonotope<Number>&			 Rtp)
	{
		return 0;
	}

	template <typename Number>
	int NonlinearSys<Number>::linReach(ReachOptions<Number>&		options,
									   ReachableSetElement<Number>& Rstart,
									   ReachableSetElement<Number>& Rti,
									   ReachableSetElement<Number>& Rtp)
	{
		// extract initial set and abstraction error
		// std::cout << "this is linReach" << std::endl;
		Zonotope<Number> Rinit = Rstart.rs();

		Vector_t<Number> abstrerr	  = Rstart.error();
		Zonotope<Number> Rti_internal = Zonotope<Number>();
		Zonotope<Number> Rtp_internal = Zonotope<Number>();

		// linearize the nonlinear system
		ReachOptions<Number> linOptions;
		LinearSys<Number>	 linSys = linearize(options, Rinit, linOptions);

		// translate Rinit by linearization point
		Zonotope<Number> Rdelta = Rinit - linerror_.p.x;
		// std::cout<<"her1231231e"<<std::endl;

		// compute reachable set of the linearized system
		LinearReachableSet<Number> R = linSys.initReach(Rdelta, linOptions);
		// std::cout<<"here1"<<std::endl;
		if (options.alg() == std::string("poly"))
		{
			// polyZonotope<Number> Rdiff = linSys.deltaReach(options, Rdelta);
			/*if(options.tensor_order > 2){
        //precompStatError;
    }*/
		}

		// compute reachable set of the abstracted system including the abstraction
		// error using the selected algorithm
		double perfInd;
		if (options.alg() == std::string("linRem"))
		{
			perfInd = linReach_linRem(R, Rinit, Rdelta, options, Rti_internal, Rtp_internal);
		}
		else
		{
			// loop until the actual abstraction error is smaller than the estimated
			// linearization error
			Rtp_internal	   = R.time_point();
			Rti_internal	   = R.time_interval();
			double perfIndCurr = DBL_MAX;

			Zonotope<Number> VerrorDyn;
			Zonotope<Number> VerrorStat;
			perfInd = 0;
			// std::cout<<"here2"<<std::endl;
			int thistime = 0;
			// std::cout<<"here3"<<std::endl;
			while (perfIndCurr > 1 && perfInd <= 1)
			{
				//  estimate the abstraction error
				Vector_t<Number> appliedError = abstrerr * 1.1;
				/*cout<<"appliedError"<<endl;
      cout<<appliedError<<endl;*/

				Zonotope<Number> Verror	   = Zonotope<Number>(appliedError * 0, appliedError.asDiagonal());
				Zonotope<Number> RallError = linSys.error_solution(options, Verror, Zonotope<Number>());
				// compute the abstraction error using the conservative linearization
				// approach described in [1]
				Vector_t<Number> trueError;
				if (options.alg() == std::string("lin"))
				{
					// compute overall reachable set including linearization error
					Zonotope<Number> Rmax = Rti_internal + RallError;
					/*std::cout<<"Rmax"<<std::endl;
        std::cout<<Rmax<<std::endl;
        std::cout<<"Rti_internal"<<std::endl;
        std::cout<<Rti_internal<<std::endl;
        std::cout<<"RallError"<<std::endl;
        std::cout<<RallError<<std::endl;*/
					// compute linearization error
					VerrorDyn				   = Zonotope<Number>();
					IntervalMatrix<Number> IHx = Rmax.interval();
					// compute intervals of total reachable set
					IntervalMatrix<Number> totalInt_x;
					totalInt_x = IHx + setIntervalMatrix(linerror_.p.x);

					/*totalInt_x.inf = IHx.inf+linerror_.p.x;
        totalInt_x.sup = IHx.sup+linerror_.p.x;*/
					trueError  = abstrerr_lin(options, Rmax, VerrorDyn);
					VerrorStat = Zonotope<Number>();
				}
				else
				{
					// compute overall reachable set including linearization error
					// Zonotope<Number> Rmax = Rdelta + RallError + Rdiff.toZonotope();
					// compute abstraction error
					// abstrerr_poly
				}
				// compare linearization error with the maximum allowed error
				// perfIndCurr = (trueError.array() / appliedError.array()).maxCoeff();
				// // std::cout<<"perfInd"<<std::endl;
				// perfInd = (trueError.array() / options.max_error().array()).maxCoeff();
				// // std::cout<<"abstrerr"<<std::endl;
				// abstrerr = trueError;
				// // if any(abstrerr > 1e+100)
				// //     throw(SetExplosionException());

				Vector_t<Number>  perfIndCurrArray = trueError.array() / appliedError.array();
				for(int i = 0; i < perfIndCurrArray.rows(); i++){
					if(isnan(perfIndCurrArray(i))){
						perfIndCurrArray(i) = -DBL_MAX;
					}
				}
				Vector_t<Number>  perfIndArray = trueError.array() / options.max_error().array();
				for(int i = 0; i < perfIndArray.rows(); i++){
					if(isnan(perfIndArray(i))){
						perfIndArray(i) = -DBL_MAX;
					}
				}
				perfIndCurr = perfIndCurrArray.maxCoeff();
				perfInd		= perfIndArray.maxCoeff();
				abstrerr	= trueError;

				// cout << options.max_error() << endl;
				// cout << perfIndCurr << endl;
				// cout << perfInd << endl;
				// cout << trueError << endl;
			} // while
			// translate reachable sets by linearization point
			Rti_internal = Rti_internal + linerror_.p.x;
			Rtp_internal = Rtp_internal + linerror_.p.x;
			// compute the reachable set due to the linearization error
			// std::cout<<"error_solution"<<std::endl;
			Zonotope<Number> Rerror = linSys.error_solution(options, VerrorDyn, VerrorStat);
			// add the abstraction error to the reachable sets
			Rti_internal = Rti_internal + Rerror;
			Rtp_internal = Rtp_internal + Rerror;
		} // else

		// determine the best dimension to split the set in order to reduce the
		// linearization error

		// std::cout<<"dimForsplit"<<std::endl;
		int dimForSplit = -1;

		// cout << perfInd << endl;

		if (perfInd > 1)
		{
			// std::cout<<"select"<<std::endl;
			// cout << "perfInd " << perfInd << endl;
			dimForSplit = select(options, Rstart);

			//cout << dimForSplit << endl;
		}
		// store the linearization error
		Rtp.set_rs(Rtp_internal);
		Rtp.set_error(abstrerr);
		Rti.set_rs(Rti_internal);
		// std::cout << "end linReach" << std::endl;
		return dimForSplit;
	} // linReach

	template <typename Number>
	int NonlinearSys<Number>::linReach(polyReachOptions<Number>&		options,
									   polyReachableSetElement<Number>& Rstart,
									   polyReachableSetElement<Number>& Rti,
									   polyReachableSetElement<Number>& Rtp)
	{
		// extract initial set and abstraction error
		std::cout << "this is linReach" << std::endl;
		clock_t start, end;

		polyZonotope<Number> Rinit = Rstart.rs();

		Vector_t<Number>	 abstrerr	  = Rstart.error();
		polyZonotope<Number> Rti_internal = polyZonotope<Number>();
		polyZonotope<Number> Rtp_internal = polyZonotope<Number>();

		// linearize the nonlinear system
		polyReachOptions<Number> linOptions;

		// cout << Rinit <<endl;

		LinearSys<Number>		 linSys = linearize(options, Rinit, linOptions);

		// translate Rinit by linearization point
		polyZonotope<Number> Rdelta = Rinit - linerror_.p.x;
		// std::cout<<"her1231231e"<<std::endl;

		// cout << Rdelta <<endl;

		// compute reachable set of the linearized system
		polyLinearReachableSet<Number>	 R = linSys.initReach(Rdelta, linOptions);
	
		// cout << Rinit << endl;
		// cout << R.time_interval() << endl;
		// cout << R.time_point() << endl;
		// cout << endl << "HERE GOOD" << endl << endl;

		std::vector<Matrix_t<Number>>	 H;
		Zonotope<Number>				 Zdelta;
		polyZonotope<Number>			 errorStat;
		Matrix_t<IntervalMatrix<Number>> T;
		std::vector<size_t>				 ind3;
		Zonotope<Number>				 Zdelta3;
		polyZonotope<Number>			 Rdiff;
		// std::cout<<"here1"<<std::endl;

		if (options.alg() == std::string("poly"))
		{
			// cout << Rdelta << endl;
			// cout << endl << "HERE GOOD" << endl << endl;
			Rdiff = linSys.deltaReach(options, Rdelta);

			// cout << Rdiff << endl;
			// cout << endl << "HERE GOOD" << endl << endl;
			if (options.tensor_order() > 2)
			{
				start = clock();
				precompStatError(Rdelta, options, H, Zdelta, errorStat, T, ind3, Zdelta3);
				end = clock();
			}
		}

		cout << "pos1 time cost: " << (double) (end - start) / CLOCKS_PER_SEC << endl;

		// cout << endl << "HERE GOOD" << endl << endl;
		// compute reachable set of the abstracted system including the abstraction
		// error using the selected algorithm
		double perfInd;
		// if(options.alg() == std::string("linRem")){
		//     perfInd = linReach_linRem(R, Rinit, Rdelta, options, Rti_internal,
		//     Rtp_internal);
		// }else{
		// loop until the actual abstraction error is smaller than the estimated
		// linearization error
		Rtp_internal	   = R.time_point();
		Rti_internal	   = R.time_interval();
		double perfIndCurr = DBL_MAX;

		Zonotope<Number>	 VerrorDyn;
		polyZonotope<Number> VerrorStat;
		perfInd = 0;

		// cout << endl << "HERE GOOD" << endl << endl;

		start = clock();

		while (perfIndCurr > 1 && perfInd <= 1)
		{
			//  estimate the abstraction error
			Vector_t<Number> appliedError = abstrerr * 1.1;

			Zonotope<Number> Verror	   = Zonotope<Number>(appliedError * 0, appliedError.asDiagonal());
			Zonotope<Number> RallError = linSys.error_solution(options, Verror, Zonotope<Number>());
			// compute the abstraction error using the conservative linearization
			// approach described in [1]
			Vector_t<Number> trueError;

			// compute overall reachable set including linearization error
			polyZonotope<Number> Rmax = Rdelta + RallError + Rdiff.toZonotope();
			// compute abstraction error

			// cout << Rmax << Rdiff << RallError << H << Zdelta<< errorStat << T << ind3 << Zdelta3;
			// cout << Rmax << H << Zdelta<< errorStat << T << ind3 << Zdelta3;

			// cout << errorStat << endl;

			trueError = abstrerr_poly(
				options, Rmax, Rdiff + RallError, H, Zdelta, errorStat, T, ind3, Zdelta3, VerrorDyn, VerrorStat);
			//}

			// cout << VerrorDyn << endl;
			// cout << VerrorStat << endl;
			// cout << trueError << endl;
			// compare linearization error with the maximum allowed error
			Vector_t<Number>  perfIndCurrArray = trueError.array() / appliedError.array();
			for(int i = 0; i < perfIndCurrArray.rows(); i++){
				if(isnan(perfIndCurrArray(i))){
					perfIndCurrArray(i) = -DBL_MAX;
				}
			}
			Vector_t<Number>  perfIndArray = trueError.array() / options.max_error().array();
			for(int i = 0; i < perfIndArray.rows(); i++){
				if(isnan(perfIndArray(i))){
					perfIndArray(i) = -DBL_MAX;
				}
			}
			perfIndCurr = perfIndCurrArray.maxCoeff();
			perfInd		= perfIndArray.maxCoeff();
			// perfIndCurr = (trueError.array() / appliedError.array()).maxCoeff();
			// perfInd		= (trueError.array() / options.max_error().array()).maxCoeff();
			abstrerr	= trueError;

			//cout << Rmax << Rdiff + RallError << H << Zdelta<< errorStat << T << ind3 << Zdelta3;
			// cout << trueError << endl;
			// cout << appliedError << endl;
			// cout << perfIndCurr << ' ' << perfInd << endl;
			// cout << trueError.array() / appliedError.array() << endl <<trueError.array() / options.max_error().array()<< endl;
			// if any(abstrerr > 1e+100)
			//     throw(SetExplosionException());

		} // while
		// translate reachable sets by linearization point

		end = clock();
		cout << "pos2 time cost: " << (double) (end - start) / CLOCKS_PER_SEC << endl;

		Rti_internal = Rti_internal + linerror_.p.x;
		Rtp_internal = Rtp_internal + linerror_.p.x;

		// cout << Rtp_internal << endl;
		// cout << Rti_internal << endl;
		// cout << VerrorDyn << endl;
		// cout << VerrorStat << endl;

		// compute the reachable set due to the linearization error
		// std::cout<<"error_solution"<<std::endl;
		polyZonotope<Number> Rerror = linSys.error_solution(options, VerrorDyn, VerrorStat);

		// cout << Rerror;

		// add the abstraction error to the reachable sets
		if (options.alg() == std::string("poly"))
		{
			Rti_internal = Rti_internal.exactPlus(Rerror);
			Rtp_internal = Rtp_internal.exactPlus(Rerror);
		}
		else
		{
			Rti_internal = Rti_internal + Rerror;
			Rtp_internal = Rtp_internal + Rerror;
		}
		//} // else

		// determine the best dimension to split the set in order to reduce the
		// linearization error

		// std::cout<<"dimForsplit"<<std::endl;
		int dimForSplit = -1;

		// store the linearization error
		Rtp.set_rs(Rtp_internal);
		Rtp.set_error(abstrerr);
		Rti.set_rs(Rti_internal);

		// cout << Rtp_internal << endl;
		// cout << abstrerr << endl;
		// cout << Rti_internal << endl;

		std::cout << "end linReach" << std::endl;
		return dimForSplit;
	} // linReach-poly

	template <typename Number>
	Vector_t<Number> NonlinearSys<Number>::abstrerr_lin(ReachOptions<Number>& options,
														Zonotope<Number>	  R,
														Zonotope<Number>&	  VerrorDyn)
	{
		// std::cout << "this is abstrerr_lin" << std::endl;
		Vector_t<Number> trueError;
		// compute interval of reachable set
		IntervalMatrix<Number> IHx = R.interval();
		// compute intervals of total reachable set
		IntervalMatrix<Number> totalInt_x;
		// std::cout<<"R"<<std::endl;
		// std::cout<<R<<std::endl;

		// std::cout<<"IHx"<<std::endl;
		// std::cout<<IHx<<std::endl;

		// std::cout<<"p.x"<<std::endl;
		// std::cout<<setIntervalMatrix(linerror_.p.x)<<std::endl;

		totalInt_x = IHx + setIntervalMatrix(linerror_.p.x);

		/*totalInt_x.inf = IHx.inf+linerror_.p.x;
  totalInt_x.sup = IHx.sup+linerror_.p.x;*/

		// compute intervals of input
		IntervalMatrix<Number> IHu = options.U().interval();
		// translate intervals by linearization point
		IntervalMatrix<Number> totalInt_u;

		// Matrix_t<Number> tmp_matrix = Matrix_t<Number>(this->num_inputs(),1);
		// tmp_matrix << linerror_.p.u;

		totalInt_u = IHu + setIntervalMatrix(linerror_.p.u);

		/*totalInt_u.inf = IHu.inf+linerror_.p.u;
  totalInt_u.sup = IHu.sup+linerror_.p.u;*/

		/*cout<<"IHx.inf"<<endl;
  cout<<IHx.inf<<endl;
  cout<<"IHx.sup"<<endl;
  cout<<IHx.sup<<endl;*/
		if (options.tensor_order() == 2)
		{
			// obtain maximum absolute values within IHx, IHu
			Matrix_t<Number> dx = inf(IHx).cwiseAbs().cwiseMax(sup(IHx).cwiseAbs());
			Matrix_t<Number> du = inf(IHu).cwiseAbs().cwiseMax(sup(IHu).cwiseAbs());
			/*cout<<"dx"<<endl;
    cout<<dx<<endl;
    cout<<"du"<<endl;
    cout<<du<<endl;*/
			// evaluate the hessian matrix with the selected range-bounding technique
			std::vector<IntervalMatrix<Number>> H;
			// if(options.lagrangeRem().method()!='interval'){

			// }else{
			//     if(name_ == "nonlinParamSys"){

			//     }else{
			// std::cout << "before Hessian" << std::endl;
			// std::cout<<"totalInt_x"<<std::endl;
			// std::cout<<totalInt_x<<std::endl;
			// std::cout<<"totalInt_u"<<std::endl;
			// std::cout<<totalInt_u<<std::endl;
			H = hessian_(totalInt_x, totalInt_u);
			// std::cout << "after Hessian" << std::endl;
			//     }
			// }

			// calculate the Lagrange remainder (second-order error)
			Vector_t<Number> errorLagr = Eigen::MatrixXd::Zero(H.size(), 1);
			Vector_t<Number> dz		   = Matrix_t<Number>(dx.rows() + du.rows(), dx.cols());
			dz << dx, du;
			/*cout<<"dz"<<endl;
    cout<<dz<<endl;*/
			// std::cout<<"before for"<<std::endl;
			for (size_t i = 0; i < H.size(); i++)
			{
				/*cout<<"H"<<endl;
      cout<<H[i].inf<<endl;
      cout<<H[i].sup<<endl;*/
				Matrix_t<Number> H__inf, H__sup, H_max;
				H__inf = inf(H[i]).cwiseAbs();
				H__sup = sup(H[i]).cwiseAbs();
				H_max  = H__inf.cwiseMax(H__sup);

				errorLagr[i] = 0.5 * dz.adjoint() * H_max * dz;
				/*cout<<"H_max"<<endl;
      cout<<H_max<<endl;
      cout<<"errorlagr,i"<<endl;
      cout<<errorLagr[i]<<endl;*/
			}
			/* cout<<"errorLagr"<<endl;
     cout<<errorLagr<<endl;*/
			// std::cout<<"isthere?"<<std::endl;
			VerrorDyn = Zonotope<Number>(0 * errorLagr, errorLagr.asDiagonal());
			trueError = errorLagr;
			// std::cout << "end absrror" << std::endl;
			return trueError;
		}
		else if (options.tensor_order() == 3)
		{
			// TODO not complete

			// // obtain intervals and combined interval z
			// std::cout << "tensor_order = 3" << std::endl;
			IntervalMatrix<Number> dz(IHx.rows() + IHu.rows(), IHx.cols());
			dz << IHx, IHu;
			// reduce zonotope
			Zonotope<Number> Rred = R;
			// std::cout<<"Reduce"<<std::endl;
			Rred.Reduce(options.error_order());

			// combined zonotope (states + input)
			// std::cout<<"cartProd"<<std::endl;
			Zonotope<Number> Z = Rred.cartProd(options.U());

			// calculate hessian matrix
			std::vector<Matrix_t<Number>> H;
			// std::cout<<"Hessian"<<std::endl;
			H = hessian_(linerror_.p.x, linerror_.p.u);

			// evaluate third-order tensor
			Vector_t<std::vector<size_t>> index;
			// std::cout<<"totalInt_x"<<std::endl;
			// std::cout<<totalInt_x<<std::endl;
			// std::cout<<"thirdOrderTensor"<<std::endl;
			Matrix_t<IntervalMatrix<Number>> T = thirdOrderTensor(totalInt_x, totalInt_u, index);

			// second-order error
			// std::cout<<"quadMap"<<std::endl;
			Zonotope<Number> errorSec = Z.quadMap(H) * 0.5;
			// std::cout<<"endquadMap"<<std::endl;
			Zonotope<Number> errorLagr;

			// calculate the lagrange remainder term (third-order error)
			int flag = 0;
			for (auto i = 0; i < index.rows(); i++)
			{
				if (index(i).size() != 0)
				{
					flag = 1;
					break;
				}
			}
			// std::cout << "flag" << std::endl;
			// std::cout << flag << std::endl;
			if (flag == 0)
			{
				// empty zonotope if all entried in T are empty
				errorLagr = Zonotope<Number>(this->dim());
			}
			else
			{
				// skip tensors with all-zero entries using ind from tensor creation
				IntervalMatrix<Number> errorLagrtemp;
				Matrix_t<Number>	   tempzero = Eigen::MatrixXd::Zero(this->dim(), 1);
				errorLagrtemp					= setIntervalMatrix(tempzero);
				for (auto i = 0; i < index.size(); i++)
				{
					Interval_t<Number> error_sum(0., 0.);
					// std::cout<<"i "<<i<<std::endl;
					// std::cout<<"i.size "<<index[i].size()<<std::endl;
					for (auto j = 0; j < index[i].size(); j++)
					{
						IntervalMatrix<Number> error_tmp = dz.adjoint() * T(i, j) * dz;
						// std::cout<<"error_tmp"<<std::endl;
						// std::cout<<error_tmp<<std::endl;
						error_tmp = error_tmp * dz(j);
						// std::cout<<"error_tmp0.2"<<std::endl;
						// std::cout<<error_tmp<<std::endl;
						Interval_t<Number> error_tmp2 = error_tmp(0, 0);
						error_sum					  = error_sum + error_tmp2;
						// std::cout<<"error_tmp2"<<std::endl;
						// std::cout<<error_tmp2<<std::endl;
					}
					errorLagrtemp(i, 0) = 1. / 6 * error_sum;
				}
				// std::cout<<"errorLagrtemp"<<std::endl;
				// std::cout<<errorLagrtemp<<std::endl;
				errorLagr = Zonotope<Number>(errorLagrtemp);
				// std::cout<<"errorLagrqqqq"<<std::endl;
				// std::cout<<errorLagr<<std::endl;
			}

			// overall linearization error
			// std::cout<<"errorSec"<<std::endl;
			// std::cout<<errorSec<<std::endl;

			// std::cout<<"errorLagr"<<std::endl;
			// std::cout<<errorLagr<<std::endl;

			// std::cout<<"errorLagzzzzzr"<<std::endl;

			VerrorDyn = errorSec + errorLagr;
			VerrorDyn.Reduce(options.intermediate_order());

			// std::cout<<"VerrorDyn"<<std::endl;
			// std::cout<<VerrorDyn<<std::endl;

			IntervalMatrix<Number> VerrorDyntemp;
			Matrix_t<Number>	   H__inf, H__sup, H_max;
			VerrorDyntemp = VerrorDyn.interval();
			H__inf		  = inf(VerrorDyntemp).cwiseAbs();
			H__sup		  = sup(VerrorDyntemp).cwiseAbs();
			H_max		  = H__inf.cwiseMax(H__sup);

			// std::cout<<"trueError"<<std::endl;
			// std::cout<<trueError<<std::endl;
			trueError = H_max;

			// std::cout << "end absrt lin" << std::endl;

			return trueError;
		}
		else
		{
			std::cout << "No abstraction error computation for chosen tensor order!" << std::endl;
			return trueError;
		}
	}

	template <typename Number>
	void NonlinearSys<Number>::verboseLog(int step, double t, ReachOptions<Number>& options)
	{
		if (options.verbose())
		{
			if (abs(t - options.tStart()) < 1e-12)
			{
				std::cout << "Start analysis..." << std::endl;
				std::cout << "-step 1: " << t << std::endl;
				return;
			}
			if (abs(t - options.tFinal()) < 1e-12)
			{
				std::cout << "...time horizon reached, analysis finished." << std::endl;
				return;
			}
			int cycle = 10;
			if (((step - 1) % cycle) == 0)
			{
				std::cout << "- step " << step << ": " << t << std::endl;
			}
		}
	}

	template <typename Number>
	void NonlinearSys<Number>::verboseLog(int step, double t, polyReachOptions<Number>& options)
	{
		if (options.verbose())
		{
			if (abs(t - options.tStart()) < 1e-12)
			{
				std::cout << "Start analysis..." << std::endl;
				std::cout << "-step 1: " << t << std::endl;
				return;
			}
			if (abs(t - options.tFinal()) < 1e-12)
			{
				std::cout << "...time horizon reached, analysis finished." << std::endl;
				return;
			}
			int cycle = 10;
			if (((step - 1) % cycle) == 0)
			{
				std::cout << "- step " << step << ": " << t << std::endl;
			}
		}
	}

	template <typename Number>
	void NonlinearSys<Number>::precompStatError(polyZonotope<Number>			  Rdelta,
												polyReachOptions<Number>&		  options,
												std::vector<Matrix_t<Number>>&	  H,
												Zonotope<Number>&				  Zdelta,
												polyZonotope<Number>&			  errorStat,
												Matrix_t<IntervalMatrix<Number>>& T,
												std::vector<size_t>&			  ind3,
												Zonotope<Number>&				  Zdelta3)
	{

		clock_t start, end;
		start = clock();
		// reduce the reachable set for the initial time point
		polyZonotope<Number> Rred = Rdelta;
		Rred.Reduce(options.error_order());

		// compute a zonotope over-approximation of the reachable set at the initial
		// time point if( isa(Rdelta, "zonotope"){ Rdelta = Rred;
		//}
		// else{
		Zonotope<Number> Rdelta2 = Rdelta.toZonotope();
		Rdelta2.Reduce(options.error_order());

		// cout << endl << "HERE GOOD" << endl << endl;

		//}
		// extend the sets by the input sets
		Vector_t<Number> tempc = Eigen::MatrixXd::Zero(options.U().dimension(), 1);

		// cout << tempc <<endl;

		Zonotope<Number> Ustat(tempc);

		// cout << endl << "HERE GOOD" << endl << endl;

		// cout << Rred << endl;
		// cout << Ustat <<endl;

		polyZonotope<Number> Z = Rred.cartProd(Ustat);

		// cout << endl << "HERE GOOD" << endl << endl;

		Zdelta				   = Rdelta2.cartProd(Ustat);

		// cout << endl << "HERE GOOD" << endl << endl;

		// std::vector<IntervalMatrix<Number>> H;
		H = hessian_(linerror_.p.x, linerror_.p.u);

		// calculate the quadrastic map == static second order error

		// cout << Z << endl;
		// cout << H << endl;

		// clock_t start, end;
		// start = clock();

		polyZonotope<Number> errorSecOrdStat = Z.quadMap(H) * 0.5;

		end = clock();
		cout << "precompStatError time cost: " << (double) (end - start) / CLOCKS_PER_SEC << endl;
		// third-order error

		// cout << endl << "HERE GOOD" << endl << endl;

		if (options.tensor_order() >= 4)
		{
			// reduce the order of the reachable set to speed-up the computations for
			// cubic multiplication
			if (options.error_order() == 3)
			{ // maybe wrong
				Rred.Reduce(options.error_order());
				Rdelta.Reduce(options.error_order());

				Z = Rred.cartProd(options.U());

				polyZonotope<Number> Zdelta3 = Rdelta.cartProd(options.U());
			}
			else
			{
				polyZonotope<Number> Zdelta3 = Zdelta;
			}
			// if isa nonlinParamSys
			// else{
			// T,ind3 = thirdOrderTensor()
			//}

			// errorThirdOrdStat = 1/6 * cubMap(Z,T,ind3);

			// calculate the overall static linearization error
			// if is polyZonotope
			//  errrorStat = exactPlus(errorSecOrdStat,errorThirdOrdStat);
			// else{
			// errorStat = errorSecOrdStat + errorThirdOrdStat
			//}
		}
		else
		{
			errorStat = errorSecOrdStat;
		}
		// cout << endl << "HERE GOOD" << endl << endl;
		// reduce the complexity of the set of static errors
		errorStat.Reduce(options.intermediate_order());
		// cout << endl << "HERE GOOD" << endl << endl;
	}

	template <typename Number>
	Vector_t<Number> NonlinearSys<Number>::abstrerr_poly(polyReachOptions<Number>&		  options,
														 polyZonotope<Number>			  Rall,
														 polyZonotope<Number>			  Rdiff,
														 std::vector<Matrix_t<Number>>	  H,
														 Zonotope<Number>				  Zdelta,
														 polyZonotope<Number>			  VerrorStat,
														 Matrix_t<IntervalMatrix<Number>> T,
														 std::vector<size_t>			  ind3,
														 Zonotope<Number>				  Zdelta3,
														 Zonotope<Number>&				  VerrorDyn,
														 polyZonotope<Number>&			  VerrorStat2)
	{
		// compute interval of reachable set
		IntervalMatrix<Number> dx		  = Rall.interval();
		IntervalMatrix<Number> totalInt_x = dx + setIntervalMatrix(linerror_.p.x);

		// compute intervals of input
		IntervalMatrix<Number> du		  = options.U().interval();
		IntervalMatrix<Number> totalInt_u = du + setIntervalMatrix(linerror_.p.u);

		// obtain intervals and combined interval z
		IntervalMatrix<Number> dz(dx.rows() + du.rows(), dx.cols());
		dz << dx, du;

		// compute zonotope of state and input
		Zonotope<Number> Rred_diff = Rdiff.toZonotope();
		Rred_diff.Reduce(options.error_order());
		Zonotope<Number> Z_diff = Rred_diff.cartProd(options.U());

		// second-order error

		// cout << Z_diff << endl;
		// cout << Zdelta << endl;
		// cout << H << endl;

		Zonotope<Number> error_secondeOrder_dyn
			= (Zdelta.quadMapMixed(Z_diff, H) +  Z_diff.quadMapMixed(Zdelta, H) + Z_diff.quadMap(H)) * 0.5;
		Zonotope<Number> error_thirdOrder_dyn_temp;
		Zonotope<Number> remainder;

		// cout << error_secondeOrder_dyn.generators().cols() << endl;
		// cout << Zdelta.quadMapMixed(Z_diff, H) << endl;
		// cout << Z_diff.quadMap(H) << endl;

		// cout << error_secondeOrder_dyn << endl;

		// third-order error
		if (options.tensor_order() == 3)
		{
			// evaluate the third-oerder tensor
			// judege whether is interval method
			Vector_t<std::vector<size_t>> ind;

			// cout << totalInt_x << endl;
			// cout << totalInt_u << endl;
			
			T = thirdOrderTensor(totalInt_x, totalInt_u, ind);
			// calculate the lagrange remainder term
			Matrix_t<Number>	   setM					= Eigen::MatrixXd::Zero(this->dim(), 1);
			IntervalMatrix<Number> error_thirdOrder_dyn = setIntervalMatrix(setM);

			// cout << dz << endl;
			// cout << T(0,0) << endl;

			for (auto i = 0; i < ind.size(); i++)
			{
				Interval_t<Number> error_sum;
				error_sum.setLeftBound(0.);
				error_sum.setRightBound(0.);
				for (auto j = 0; j < ind(i).size(); j++)
				{
					IntervalMatrix<Number> add = (dz.adjoint() * T(i,ind(i)[j]) * dz) * dz(ind(i)[j]);
					error_sum = error_sum + add(0,0);
					//cout << error_sum;
				}
				error_thirdOrder_dyn(i, 0) = 1.0 / 6.0 * error_sum;
			}

			// cout << error_thirdOrder_dyn;

			error_thirdOrder_dyn_temp = Zonotope<Number>(error_thirdOrder_dyn);

			// no terms of order >= 4
			Vector_t<Number> temptemainder = Eigen::MatrixXd::Zero(this->dim(), 1);
			remainder					   = Zonotope<Number>(temptemainder);
		}
		else
		{
			// tensorOrder >= 4

			// reduce set Zdiff to the desired zonotope order to speed up the
			// computation of cubic multiplication Zonotope<Number> Z_diff3 = Z_diff;

			// third-order error
			// Zonotope<Number> error_thirdOrder_dyn_temp;

			// init hight-order error
			// Matrix_t<Number> setM = Eigen::MatrixXd::Zero(this->dim(),1);
			// IntervalMatrix<Number> remaindertemp = setIntervalMatrix(setM);

			// exact evaluation of intermediate taylor terms
			// for(auto i=4;i<options.tensor_order();i++){
			//    IntervalMatrix<Number> handle;
			// }

			// lagrange remainder over-approxiamting the last taylor term

			// Zonotope<Number> remainder(remaindertemp);
		}

		// combine results

		// cout << error_secondeOrder_dyn << endl;
		// cout << error_thirdOrder_dyn_temp << endl;
		// cout << remainder << endl;

		VerrorDyn = error_secondeOrder_dyn + error_thirdOrder_dyn_temp + remainder;

		// cout << VerrorDyn << endl;

		VerrorDyn.Reduce(options.intermediate_order());

		// cout << VerrorDyn << endl;
		// cout << VerrorStat << endl;

		IntervalMatrix<Number> errorIHabs = VerrorDyn.interval() + VerrorStat.interval();
		for (auto i = 0; i < errorIHabs.rows(); i++)
		{
			for (auto j = 0; j < errorIHabs.cols(); j++)
			{
				errorIHabs(i, j) = capd::intervals::iabs(errorIHabs(i, j));
			}
		}

		Vector_t<Number> trueError = sup(errorIHabs);
		VerrorStat2				   = VerrorStat;

		return trueError;
	}

	template <typename Number>
	int NonlinearSys<Number>::select(ReachOptions<Number> options, ReachableSetElement<Number>& Rinit)
	{
		// compute all possible splits of the maximum reachable set
		std::vector<std::vector<Zonotope<Number>>> Rtmp = Rinit.rs().split();

		// cout << Rtmp << endl;

		std::vector<ReachableSetElement<Number>>   R(Rtmp.size());
		for (auto i = 0; i < Rtmp.size(); i++)
		{
			R[i].set_rs(Rtmp[i][0]); // only rest one of the two split sets
			R[i].set_error(Eigen::MatrixXd::Zero(options.max_error().rows(), 1));
		}

		// adapt the options for reachability analysis
		Vector_t<Number> maxError = options.max_error();

		// cout  << "options.max_error()  " << options.max_error()  << endl;

		options.set_max_error(DBL_MAX * options.max_error());

		if (options.alg() == "linRem")
		{
			options.set_alg("lin");
		}

		// loop over all split sets
		Vector_t<Number> perfInd = Eigen::MatrixXd::Zero(R.size(), 1);
		for (auto i = 0; i < R.size(); i++)
		{
			// compute the reachable set for the splitted set
			ReachableSetElement<Number> Rtp, Rti;
			linReach(options, R[i], Rti, Rtp);

			// compute perfirmance index(mmax lin .error) for the split
			perfInd(i) = Rtp.error().cwiseQuotient(maxError).maxCoeff();
		}

		// find best performance index
		int	   dimForSplit = 0;
		Number tempper	   = perfInd(0);
		for (auto i = 1; i < perfInd.rows(); i++)
		{
			if (perfInd(i) < tempper)
			{
				tempper		= perfInd(i);
				dimForSplit = i;
			}
		}
		return dimForSplit;
	}

	template <typename Number>
	Number approxVolumeRatio(const polyZonotope<Number>& pZ, const polyReachOptions<Number>& options){

		polyZonotope<Number> polyZono = pZ;
		//special cases
		if(polyZono.Grest().cols() == 0){
			return 0;
		}

		if(polyZono.generators().cols() == 0){
			return DBL_MAX;
		}

		//判断是否用 'pca'
		if(options.volApproxMethod() == "pca"){

			//calculate state-space-transformation with pca
			Matrix_t<Number> G(polyZono.generators().rows(), 2 * polyZono.generators().cols() + 2 * polyZono.Grest().cols());
			G << polyZono.generators(), -polyZono.generators(), polyZono.Grest(), -polyZono.Grest();
			Eigen::JacobiSVD<Matrix_t<Number>> svd(G, Eigen::ComputeThinU | Eigen::ComputeThinV);  
    		Matrix_t<Number> U = svd.matrixU(); 
			Matrix_t<Number> Unew = U.adjoint();

			//transform the polynomial zonotope to the new state space
			polyZonotope<Number> polyZonoTemp = polyZono * Unew;
			polyZono = polyZonoTemp;
		}

		//over-approximate the independent generators part with an interval
		int n = polyZono.center().size();
		Vector_t<Number> zeroCenter;
		zeroCenter.setZero(n);
		Zonotope<Number> zono = Zonotope<Number>(zeroCenter, polyZono.Grest());
		IntervalMatrix<Number> Iind = zono.interval();

		//over-approximate the dependent generators part with an interval
		Matrix_t<Number> zeroGrest;
		polyZono.set_Grest(zeroGrest);

		IntervalMatrix<Number> Idep = polyZono.interval();

		// cout << Iind << endl;
		// cout << Idep << endl;

		//remove dimensions that are all-zero
		vector<Number> ind1 = findZero(rad(Idep));
		vector<Number> ind2 = findZero(rad(Iind));

		vector<Number> ind = Unique(ind1, ind2);
		vector<Number> idx(n);
  		iota(idx.begin(), idx.end(), 0);
		ind = setdiff(idx, ind);

		IntervalMatrix<Number> tempIind(ind.size(), Iind.cols());
		IntervalMatrix<Number> tempIdep(ind.size(), Idep.cols());

		for(int i = 0; i < ind.size(); i++){
			tempIind.row(i) = Iind.row(ind[i]);
			tempIdep.row(i) = Idep.row(ind[i]);
		}

		Number Vind = Volume(tempIind);
		Number Vdep = Volume(tempIdep);

		// cout << Vind << endl;
		// cout << Vdep << endl;

		Number ratio = pow(Vind/Vdep,1/(double)n);
		//Number ratio = 0;
		return ratio;
	}

} // namespace reachSolver
