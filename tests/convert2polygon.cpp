/*
 * @Author: your name
 * @Date: 2021-11-11 17:24:27
 * @LastEditTime: 2021-12-28 18:18:30
 * @LastEditors: Please set LastEditors
 * @Description: 打开koroFileHeader查看配置 进行设置:
 * https://github.com/OBKoro1/koro1FileHeader/wiki/%E9%85%8D%E7%BD%AE
 * @FilePath: /Solver/test.cpp
 */


#include <zonotope/Zonotope.h>
#include <eigen3/Eigen/Core>
//#include <functional>
#include <iostream>
#include<cmath>


using namespace reachSolver;

template <typename Number>
void RemoveColumn(Matrix_t<Number> & matrix, unsigned int colToRemove) {
 	unsigned int numRows = matrix.rows();
  	unsigned int numCols = matrix.cols() - 1;
 
  	if( colToRemove < numCols ) {
    	matrix.block(0, colToRemove, numRows, numCols - colToRemove) =
      	matrix.block(0, colToRemove + 1, numRows, numCols - colToRemove);
  	}
 
  	matrix.conservativeResize(numRows,numCols);
}

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

int main()
{
	
	Vector_t<double> c(2);
	Matrix_t<double> matrix(2,7);
	c << -4.5, 8.44;
	matrix << 2.523, 0, 0, 4.35, 2, 2, 1,
		 0, 0, 0, -1, 2, 0, 2;
	Zonotope<double> Z(c,matrix);
	Matrix_t<double> p = convert2polygon(Z);
	cout << p;
	/*
	//delete zero generators
	deleteZeros(matrix);

	unsigned int colNum = matrix.cols();//obtain number of generators
	
	//obtain size of enclosing intervalhull of first two dimensions
	double xmax = matrix.cwiseAbs().rowwise().sum()(0,0);
	double ymax = matrix.cwiseAbs().rowwise().sum()(1,0);
	//cout << xmax << " " << ymax;

	//Z with normalized direction: All generators pointing "up"
	Matrix_t<double> matrixNorm = matrix;
	for(auto i = 0; i < colNum; i++){
		if(matrixNorm(1,i) < 0){
			matrixNorm(0,i) = matrixNorm(0,i) * -1;
			matrixNorm(1,i) = matrixNorm(1,i) * -1;
		}
	}
	//cout << matrixNorm;
	//compute angles
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
	/*
	for(auto i = 0; i < colNum; i++){
		cout << Index[i] << " ";
	}
	*/

	//cumsum the generators in order of angle
	/*
	vector<double> p_1(colNum+1);
	vector<double> p_2(colNum+1);

	for(auto i = 0; i < colNum; i++){
		p_1[i+1] = p_1[i] + 2 * matrixNorm(0,Index[i]);
		p_2[i+1] = p_2[i] + 2 * matrixNorm(1,Index[i]);
	}

	double p_1Max = *max_element(p_1.begin(), p_1.end());
	for(auto i = 0; i < colNum+1; i++){
		p_1[i] = p_1[i] + xmax - p_1Max;
		p_2[i] = p_2[i] - ymax;
	}

	/*
	for(auto i = 0; i < colNum+1; i++){
		cout << p_1[i]<< " ";
	}
	cout << endl;
	for(auto i = 0; i < colNum+1; i++){
		cout << p_2[i]<< " ";
	}
	*/
	//flip/mirror upper half to get lower half of zonotope (point symmetry)
	//consider center
	/*
	Matrix_t<double> p(2,2*colNum+1);
	for(auto i = 0; i < colNum+1; i++){
		p(0,i) = p_1[i] + c(0);
		p(1,i) = p_2[i] + c(1);
	}
	for(auto i = colNum+1; i < 2*colNum+1; i++){
		p(0,i) = p_1[colNum] + p_1[0] - p_1[i-colNum] + c(0);
		p(1,i) = p_2[colNum] + p_2[0] - p_2[i-colNum] + c(1);
	}
	cout << p;
	*/
}
