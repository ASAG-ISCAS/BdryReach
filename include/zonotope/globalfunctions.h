/**
 * @file   globalfunctions.h
 * @brief  global functions
 * @author Changyuan Zhao
 * @date  Dec 2021
 * @version 1.0
 *
 * Reference:
 *   CORA ../global/functions
 */

#ifndef GLOBALFUNCTIONS_H
#define GLOBALFUNCTIONS_H

#include <iostream>
#include <numeric>
#include <zonotope/BasicObject.h>
#include <set>

using namespace std;

namespace reachSolver {

/**
 * @brief: sort rows of a matrix in ascending order as a group
 * @param A a matrix must be a 2-D matrix
 * @return I the index
 */
// template <typename Number> std::vector<size_t> sortrows(Matrix_t<Number> A) {

//   // clock_t start, end;
//   // start = clock();
//   // cout << "rows: " << A.rows() << "cols: " << A.cols() << endl;
//   std::vector<size_t> idx(A.rows());
//   std::iota(idx.begin(), idx.end(), 0);
//   std::sort(idx.begin(), idx.end(),
//             [&A](size_t i1, size_t i2) { 
//               int col = 0;
//               //while(A(i1, col) - A(i2, col) <= 1e-15 && col < A.cols() - 1){
//               while(abs(A(i1, col) - A(i2, col)) <= 1e-12 && col < A.cols() - 1){
//                 col ++;
//               }
//               return A(i1, col) < A(i2, col); });

//   // end = clock();
// 	// cout << "sortrows time cost: " << (double) (end - start) / CLOCKS_PER_SEC << endl;
//   return idx;

// }

template <typename Number> void sortrows_oneCol(Matrix_t<Number>& A, int col) {

  Matrix_t<Number> res(A.rows(), A.cols());
  std::vector<size_t> idx(A.rows());
  std::iota(idx.begin(), idx.end(), 0);
  std::sort(idx.begin(), idx.end(),
            [&A, &col](size_t i1, size_t i2) { 
              return A(i1, col) < A(i2, col); });

  for(int i = 0; i < A.rows(); i++){
    res.row(i) = A.row(idx[i]);
  }

  A = res;

}


template <typename Number> std::vector<size_t> sortrows(const Matrix_t<Number>& A) {
  // clock_t start, end;
  // start = clock();
  //   cout << "rows: " << A.rows() << "cols: " << A.cols() << endl;
  std::vector<size_t> idx(A.rows());
  std::iota(idx.begin(), idx.end(), 0);
  int col = 0;
  std::vector<size_t> firstPos;
  std::vector<size_t> lastPos;
  std::vector<size_t> firstPosNew;
  std::vector<size_t> lastPosNew;
  firstPos.push_back(0);
  lastPos.push_back(A.rows() - 1);

  while(firstPos.size() != 0 && col < A.cols()){
    firstPosNew.clear();
    lastPosNew.clear();
    for(int i = 0; i < firstPos.size(); i++){
        std::sort(idx.begin() + firstPos[i], idx.begin() + lastPos[i] + 1,
            [&A, &col](size_t i1, size_t i2) {return A(i1, col) < A(i2, col); });

        int ref = A(idx[firstPos[i]], col);
        int j = firstPos[i] + 1;
        while(j <= lastPos[i]){
            if(abs(A(idx[j],col) - ref) <= 1e-12){
              firstPosNew.push_back(j - 1);
              while(j <= lastPos[i]){
                if(abs(A(idx[j],col) - ref) <= 1e-12){
                  j++;
                }else{
                  ref = A(idx[j],col);
                  break;
                }
              }
              lastPosNew.push_back(j - 1);
            }else{
              ref = A(idx[j],col);
            }
            j++;
        }
    }
    // cout << firstPosNew << endl;
    // cout << lastPosNew << endl;
    // cout << idx << endl;
    firstPos = firstPosNew;
    lastPos = lastPosNew;
    col++;
  } 
  // end = clock();
	// cout << "sortrows time cost: " << (double) (end - start) / CLOCKS_PER_SEC << endl;
  return idx;
}
// template <typename Number>
// struct A_col{
//   Matrix_t<Number>& A;
//   int col;
//   int row;

//   bool operator < (const A_col& other) const {
//       return A()
//   }

//   }
// };

/**
 * @brief: add up all generators that belong to terms with identical exponents
 * @param ExpMat matrix containing the exponent vectors
 * @param G generator matrix
 * @param ExpMatNew modified exponent matrix
 * @param Gnew modified generator matrix
 */
template <typename Number>
void removeRedundantExponents(Matrix_t<int> ExpMat, Matrix_t<Number> G,
                              Matrix_t<int> &ExpMatNew,
                              Matrix_t<Number> &Gnew) {
  // clock_t start, end;
  // start = clock();

  // remove zero-length generators
  Eigen::MatrixXd idxD = G.colwise().any();

  // std::cout << G << std::endl;
  // std::cout << idxD << std::endl;
  // std::cout << idxD.all() << std::endl;
  // std::cout << !(idxD.all()) << std::endl;
  // std::cout << (idxD.prod() == 0) << std::endl;

  // skip if all non-zero
  if (!(idxD.all())) {
    // if((~idxD).all()){
    if (idxD.maxCoeff() == 0) {
      ExpMatNew.setZero(ExpMat.rows(), 1);
      Gnew = Eigen::MatrixXd::Zero(G.rows(), 1);
      return;
    } else {
      auto j = 0;
      Matrix_t<Number> tempG(G.rows(), G.cols());
      for (auto i = 0; i < G.cols(); i++) {
        if (idxD(i) == true) {
          tempG.col(j) = G.col(i);
          j++;
        }
      }
      tempG.conservativeResize(G.rows(), j);
      //tempG.resize(G.rows(), j);
      j = 0;
      Matrix_t<int> tempExpMat(ExpMat.rows(), ExpMat.cols());
      for (auto i = 0; i < ExpMat.cols(); i++) {
        if (idxD(i) == true) {
          tempExpMat.col(j) = ExpMat.col(i);
          j++;
        }
      }
      tempExpMat.conservativeResize(ExpMat.rows(), j);
      //tempExpMat.resize(ExpMat.rows(), j);
      // std::cout << tempExpMat << std::endl;
      // std::cout << tempG << std::endl;
      ExpMat = tempExpMat;
      G = tempG;
    }
  }

  // std::cout << ExpMat << std::endl;
  // std::cout << G << std::endl;

  // add hash value of the exponent vectors to the exponent matrix
  Matrix_t<int> temp(1, ExpMat.rows());
  for (auto i = 0; i < ExpMat.rows(); i++) {
    temp(0, i) = (int)i + 1;
  }
  Matrix_t<int> temp2 = (temp * ExpMat).adjoint();
  Matrix_t<int> rankMat(temp2.rows(), temp2.cols() + ExpMat.rows());
  rankMat << temp2, ExpMat.adjoint();

  // sort the exponent vectors according to the hash value

  // std::cout << rankMat << std::endl;

  std::vector<size_t> index = sortrows(rankMat);

  // std::cout << index << std::endl;
  
  Matrix_t<int> ExpMatTemp(ExpMat.rows(), ExpMat.cols());
  for (auto i = 0; i < ExpMat.cols(); i++) {
    ExpMatTemp.col(i) = ExpMat.col(index[i]);
  }
  // Matrix_t<Number> ExpMatTemp = ExpMat(index);
  Matrix_t<Number> Gtemp(G.rows(), G.cols());
  for (auto i = 0; i < G.cols(); i++) {
    Gtemp.col(i) = G.col(index[i]);
  }

  // std::cout << ExpMatTemp << std::endl;
  // std::cout << Gtemp << std::endl;

  // initialization
  int counterNew = 0;
  ExpMatNew.setZero(ExpMat.rows(), ExpMat.cols());
  Gnew = Eigen::MatrixXd::Zero(G.rows(), G.cols());

  // first entry
  ExpMatNew.col(counterNew) = ExpMatTemp.col(0);
  Gnew.col(counterNew) = Gtemp.col(0);

  // loop over all exponent vectors
  for (auto counter = 1; counter < ExpMatTemp.cols(); counter++) {
    // flag = 0;
    // for(auto i=0;i<ExpMatNew.rows();i++){
    //     if(ExpMatNew.col(counterNew) == ExpMatTemp(counter))
    // }
    if ((ExpMatNew.col(counterNew) == ExpMatTemp.col(counter))) {
      Gnew.col(counterNew) = Gnew.col(counterNew) + Gtemp.col(counter);
    } else {
      counterNew++;
      Gnew.col(counterNew) = Gtemp.col(counter);
      ExpMatNew.col(counterNew) = ExpMatTemp.col(counter);
    }
  }

  // std::cout << ExpMatNew << std::endl;
  // std::cout << Gnew << std::endl;

  // truncate exponent and generator matrix
  Matrix_t<int> ExpMatNewtemp = ExpMatNew;
  Matrix_t<Number> Gnewtemp = Gnew;
  ExpMatNew = ExpMatNewtemp.middleCols(0, counterNew + 1);
  Gnew = Gnewtemp.middleCols(0, counterNew + 1);

  // end = clock();
	// cout << "time cost: " << (double) (end - start) / CLOCKS_PER_SEC << endl;
}

/**
 * @brief: filters out generators of length0
 * @param G a matrix
 * @return reduced matrix of generators
 */
template <typename Number> Matrix_t<Number> nonzeroFilter(Matrix_t<Number> G) {
  Matrix_t<Number> result(G.rows(), G.cols());
  Vector_t<bool> index;
  index.setZero(G.cols(), 1);
  for (auto i = 0; i < G.cols(); i++) {
    if (G.col(i).any()) {
      index(i) = true;
    }
  }
  auto k = 0;
  for (auto i = 0; i < G.cols(); i++) {
    if (index(i) == true) {
      result.col(k) = G.col(i);
      k++;
    }
  }
  result.conservativeResize(result.rows(), k);
  return result;
}

/**
 * @brief: calculate the inf of interval matrix M
 * @param M intervalMatrix
 * @return a Matrix
 */
template <typename Number> Matrix_t<Number> inf(IntervalMatrix<Number> M) {
  Matrix_t<Number> result(M.rows(), M.cols());
  for (auto i = 0; i < M.rows(); i++) {
    for (auto j = 0; j < M.cols(); j++) {
      result(i, j) = M(i, j).leftBound();
    }
  }
  return result;
}

/**
 * @brief: calculate the sup of interval matrix M
 * @param M intervalMatrix
 * @return a Matrix
 */
template <typename Number> Matrix_t<Number> sup(IntervalMatrix<Number> M) {
  Matrix_t<Number> result(M.rows(), M.cols());
  for (auto i = 0; i < M.rows(); i++) {
    for (auto j = 0; j < M.cols(); j++) {
      result(i, j) = M(i, j).rightBound();
    }
  }
  return result;
}

/**
 * @brief: set interval matrix with inf and sup
 * @param  inf a matrix
 * @param  sup a matrix
 * @return intervalMatrix
 */
template <typename Number>
IntervalMatrix<Number> setIntervalMatrix(Matrix_t<Number> inf,
                                         Matrix_t<Number> sup) {
  IntervalMatrix<Number> result(inf.rows(), inf.cols());
  for (auto i = 0; i < inf.rows(); i++) {
    for (auto j = 0; j < inf.cols(); j++) {
      result(i, j).setLeftBound(inf(i, j));
      result(i, j).setRightBound(sup(i, j));
    }
  }
  return result;
}

/**
 * @brief: set interval matrix with a matrix
 * @param  M a matrix
 * @return intervalMatrix
 */
template <typename Number>
IntervalMatrix<Number> setIntervalMatrix(Matrix_t<Number> M) {
  IntervalMatrix<Number> result(M.rows(), M.cols());
  for (auto i = 0; i < M.rows(); i++) {
    for (auto j = 0; j < M.cols(); j++) {
      result(i, j).setLeftBound(M(i, j));
      result(i, j).setRightBound(M(i, j));
    }
  }
  return result;
}

/**
 * @brief: set interval matrix with inf and sup
 * @param  inf a Vector
 * @param  sup a Vector
 * @return intervalMatrix
 */
template <typename Number>
IntervalMatrix<Number> setIntervalMatrix(Vector_t<Number> inf,
                                         Vector_t<Number> sup) {
  IntervalMatrix<Number> result(inf.rows(), inf.cols());
  for (auto i = 0; i < inf.rows(); i++) {
    for (auto j = 0; j < inf.cols(); j++) {
      result(i, j).setLeftBound(inf(i, j));
      result(i, j).setRightBound(sup(i, j));
    }
  }
  return result;
}

/**
 * @brief: set interval matrix with a Vector
 * @param  M a Vector
 * @return intervalMatrix
 */
template <typename Number>
IntervalMatrix<Number> setIntervalMatrix(Vector_t<Number> M) {
  IntervalMatrix<Number> result(M.rows(), M.cols());
  for (auto i = 0; i < M.rows(); i++) {
    for (auto j = 0; j < M.cols(); j++) {
      result(i, j).setLeftBound(M(i, j));
      result(i, j).setRightBound(M(i, j));
    }
  }
  return result;
}

template <typename Number>
void setinf(IntervalMatrix<Number> &M, Matrix_t<Number> inf) {
  if (M.size() == 0) {
    M.resize(inf.rows(), inf.cols());
  }
  for (auto i = 0; i < M.rows(); i++) {
    for (auto j = 0; j < M.cols(); j++) {
      M(i, j).setLeftBound(inf(i, j));
    }
  }
}

template <typename Number>
void setsup(IntervalMatrix<Number> &M, Matrix_t<Number> sup) {
  if (M.size() == 0) {
    M.resize(sup.rows(), sup.cols());
  }
  for (auto i = 0; i < M.rows(); i++) {
    for (auto j = 0; j < M.cols(); j++) {
      M(i, j).setLeftBound(sup(i, j));
    }
  }
}

/**
 * @brief: Merge the ID vectors of two polyZonotope objects and adpapte the
 * exponent matrices accordingly
 * @param id1 IDvector of first polynomial zonotope
 * @param id2 IDvector of second polynomial zonotope
 * @param expMat1 exponent matrix of the first polynomial zonotope
 * @param expMat2 exponent matrix of the second polynomial zonotope
 * @param Mat1 adapted xponent matrix of the first polynomial zonotope
 * @param Mat2 adapted exponent matrix of the second polynomial zonotope
 * @return id merged of first polynomial zonotope
 */
template <typename Number>
std::vector<size_t> mergeExpMatrix(std::vector<size_t> id1,
                                   std::vector<size_t> id2,
                                   Matrix_t<int> expMat1, Matrix_t<int> expMat2,
                                   Matrix_t<int> &Mat1, Matrix_t<int> &Mat2) {

  int L1 = id1.size();
  int L2 = id2.size();

  std::vector<size_t> id;

  // ID vector are identical
  if (L1 == L2 && id1 == id2) {
    id = id1;
    Mat1 = expMat1;
    Mat2 = expMat2;
  }
  // ID vectors not identical->merge
  else {

    // merge the two sets
    id = id1;
    Vector_t<size_t> ind2;
    ind2.setZero(id2.size());
    for (auto i = 0; i < id2.size(); i++) {
      int ind = -1;
      for (auto j = 0; j < id.size(); j++) {
        if (id[j] == id2[i]) {
          ind = j;
          break;
        }
      }
      if (ind == -1) {
        id.push_back(id2[i]);
        ind2(i) = id.size() - 1;
      } else {
        ind2(i) = ind;
      }
    }

    // construct the new exponent matrices

    int L = id.size();
    Mat1.resize(expMat1.rows() + L - L1, expMat1.cols());
    Matrix_t<int> tempMa1;

    if ((L - L1) != 0) {
      tempMa1.setZero(L - L1, expMat1.cols());
      Mat1 << expMat1, tempMa1;
    } else {
      Mat1 << expMat1;
    }

    //std::cout << ind2 << std::endl;

    Matrix_t<int> temp;
    temp.setZero(L, expMat2.cols());

    for (auto i = 0; i < expMat2.rows(); i++) {
      temp.row(ind2(i)) = expMat2.row(i);
    }
    Mat2 = temp;
  }
  return id;
}

template <typename Number>
std::ostream &operator<<(std::ostream &os, const Vector_t<Number> &vec) {
  if (vec.rows() != 0) {
    os << "vector: deminsion " << vec.rows() << std::endl;
    for (auto i = 0; i < vec.rows(); i++) {
      os << vec(i, 0) << " ";
    }
    os << std::endl;
  }
  return os;
}

template <typename Number>
std::ostream &operator<<(std::ostream &os, const Matrix_t<Number> &mat) {
  os.precision(16);
  if (mat.rows() != 0) {
    //os << "matrix: rows and cols " << mat.rows() << " " << mat.cols()
    //   << std::endl;

    os << "[";
    for (auto i = 0; i < mat.rows(); i++) {
      for (auto j = 0; j < mat.cols(); j++) {
        os << mat(i, j) << " ";
      }
      os << ";";

      if( i < mat.rows() -1){
        os << std::endl;
      }
    }

    // os << "]" << endl << endl;
    os << "]";

  }
  return os;
}

template <typename Number>
std::ostream &operator<<(std::ostream &os, const IntervalMatrix<Number> &mat) {
  if (inf(mat).rows() != 0) {
    os << "interval matrix: inf" << std::endl;
    os << inf(mat) << std::endl;
    os << "interval matrix: sup" << std::endl;
    os << sup(mat) << std::endl;
  }
  return os;
}

//returns the radius of an interval
template <typename Number> Matrix_t<Number> rad(IntervalMatrix<Number> M) {
  Matrix_t<Number> result(M.rows(), M.cols());
  for (auto i = 0; i < M.rows(); i++) {
    for (auto j = 0; j < M.cols(); j++) {
      result(i, j) = 0.5 * (M(i, j).rightBound() - M(i, j).leftBound());
    }
  }
  return result;
}

template <typename Number> std::vector<Number> findZero(Matrix_t<Number> M) {
  std::vector<Number> res;
  Matrix_t<Number> zero;
  zero.setZero(M.rows(), M.cols());
  Matrix_t<bool> ifZero = (M.array() == zero.array());

  for (auto j = 0; j < M.cols(); j++) {
    for (auto i = 0; i < M.rows(); i++) {
      if(ifZero(i,j) == 1){
        res.push_back(j * M.rows() + i);
      }
    }
  }
  return res;
}

template <typename Number> std::vector<Number> findOnes(Matrix_t<Number> M) {
  std::vector<Number> res;
  Matrix_t<Number> ones;
  ones.setOnes(M.rows(), M.cols());
  Matrix_t<bool> ifZero = (M.array() == ones.array());

  for (auto j = 0; j < M.cols(); j++) {
    for (auto i = 0; i < M.rows(); i++) {
      if(ifZero(i,j) == 1){
        res.push_back(j * M.rows() + i);
      }
    }
  }
  return res;
}

template <typename Number> std::vector<Number> findPositive(Matrix_t<Number> M) {
  std::vector<Number> res;
  Matrix_t<Number> zero;
  zero.setZero(M.rows(), M.cols());
  Matrix_t<bool> ifZero = (M.array() > zero.array());

  for (auto j = 0; j < M.cols(); j++) {
    for (auto i = 0; i < M.rows(); i++) {
      if(ifZero(i,j) == 1){
        res.push_back(j * M.rows() + i);
      }
    }
  }
  return res;
}

template <typename Number> std::vector<Number> findNumber(Matrix_t<Number> M, Number k) {
  std::vector<Number> res;
  Matrix_t<Number> ones;
  ones.setOnes(M.rows(), M.cols());
  Matrix_t<Number> sortNumber = k * ones;
  Matrix_t<bool> ifZero = (M.array() == sortNumber.array());

  for (auto j = 0; j < M.cols(); j++) {
    for (auto i = 0; i < M.rows(); i++) {
      if(ifZero(i,j) == 1){
        res.push_back(j * M.rows() + i);
      }
    }
  }
  return res;
}

template <typename Number> std::vector<Number> Unique(vector<Number> a, vector<Number> b) {
  vector<Number> res;
  set<Number> set;
  for(auto i : a){
    if(!set.count(i)){
      set.insert(i);
    }
  }
  for(auto i : b){
    if(!set.count(i)){
      set.insert(i);
    }
  }
  for(auto i : set){
    res.push_back(i);
  }
  return res;
}
//C = SETDIFF(A,B) for vectors A and B, returns the values in A that are not in B with no repetitions.
// template <typename Number> vector<Number> setdiff(vector<Number> a, vector<Number> b){
//   int i = 0;
//   int j = 0;
//   set<Number> set;
//   vector<Number> res;
//   for(auto i : a){
//     set.insert(i);
//   }

//   while(i < a.size() && j < b.size()){
//     if(a[i] == b[j]){
//       set.erase(a[i]);
//       i++;
//       j++;
//     }else if(a[i] > b[j]){
//       j++;
//     }else{
//       i++;
//     }
//   }

//   for(auto i : set){
//     res.push_back(i);
//   }
//   return res;
// }

template <typename Number> vector<Number> setdiff(vector<Number> a, vector<Number> b){
  vector<Number> res = a;
  typename std::vector<Number>::reverse_iterator r_itA = res.rbegin();
  typename std::vector<Number>::reverse_iterator r_itB = b.rbegin();
  while(r_itA != res.rend() && r_itB != b.rend()){
    if(*r_itA == *r_itB){
      res.erase((++r_itA).base());
      ++r_itB;
    }else if(*r_itA > *r_itB){
      ++r_itA;
    }else{
      ++r_itB;
    }
  }
  return res;
}

template <typename Number> Number Volume(IntervalMatrix<Number> obj){
  Number res;
  if(obj.cols() != 0){
    Matrix_t<Number> radtemp = 2 * rad(obj);
    Matrix_t<Number> temp = radtemp.colwise().prod();
    res = temp(0,0);  
  }else{
    res = 0;
  }
  return res;
}

template <typename Number>
	std::vector<size_t> sort_indexes_ascend(const std::vector<Number> v)
	{
		std::vector<size_t> idx(v.size());
		std::iota(idx.begin(), idx.end(), 0);
		std::sort(idx.begin(), idx.end(), [&v](size_t i1, size_t i2) { return v[i1] < v[i2]; });
		return idx;
	}
template <typename Number>
	Matrix_t<Number> cov(const Matrix_t<Number> input)
  {
    //输入为input，输出为input的协方差矩阵covMat
    Matrix_t<Number> covMat;

    Matrix_t<Number> meanVec = input.colwise().mean();
    RowVector_t<Number> meanVecRow(RowVector_t<Number>::Map(meanVec.data(),input.cols()));

    Matrix_t<Number> zeroMeanMat = input;
    zeroMeanMat.rowwise()-=meanVecRow;
    if(input.rows()==1)
        covMat = (zeroMeanMat.adjoint()*zeroMeanMat)/double(input.rows());
    else
        covMat = (zeroMeanMat.adjoint()*zeroMeanMat)/double(input.rows()-1);

    return covMat;
  }
  template <typename Number>
	bool ismembertol(Matrix_t<Number> const & A, Matrix_t<Number> const & B){
    vector<size_t> idxA = sortrows(A);
    vector<size_t> idxB = sortrows(B);

    Matrix_t<Number> Asort(A.rows(), A.cols());
    Matrix_t<Number> Bsort(B.rows(), B.cols());

    for(int i = 0; i < A.rows(); i++){
      Asort.row(i) = A.row(idxA[i]);
    } 
    for(int i = 0; i < B.rows(); i++){
      Bsort.row(i) = B.row(idxB[i]);
    }
    // cout << Asort << endl;
    // cout << Bsort << endl;
    int i = 0;
    int j = 0;

    while(i < Asort.rows() && j < Bsort.rows()){
      if(abs(Asort(i,0) - Bsort(j,0)) <= 1e-12){
        // cout << "GoOD 1" << endl;
        if(abs(Asort(i,1) - Bsort(j,1)) <= 1e-12){
          // cout <<  "GoOD 2" << endl;
          return 1;
        }else if(Asort(i, 1) > Bsort(j, 1)){
          j++;
        }else{
          i++;
        }
      }else if(Asort(i, 0) > Bsort(j, 0)){
        j++;
      }else{
        i++;
      }
    }
    return 0; 
  }

  //计算阶乘
/*
注意：0！= 1
*/
int factorial(int n)
{
    int fc=1;
    for(int i=1;i<=n;++i) fc *= i;
    return fc;
}

//计算组合数
/*从n个不同元素中取出m(m≤n)个元素的所有组合的个数，叫做n个不同元素中取出m个元素的组合数。用符号c(n,m) 表示。
组合数公式：c(n,m)=n!/(m! * (n-m)!)
性质：c(n,m) = c(n,m-n)
递推公式：c(n,m) = c(n-1,m-1) + c(n-1,m)
*/
int combo(int n,int m)
{
    int com=factorial(n)/(factorial(m)*factorial(n-m));
    return com;
}

// template <typename Number> Matrix_t<Number> isZero(IntervalMatrix<Number> M) {
//   Matrix_t<Number> result(M.rows(), M.cols());
//   for (auto i = 0; i < M.rows(); i++) {
//     for (auto j = 0; j < M.cols(); j++) {
//       result(i, j) = 0.5 * (M(i, j).rightBound() - M(i, j).leftBound());
//     }
//   }
//   return result;
// }
    template <typename Number, typename vecNumber>
	  void RemoveColumn_Vec(Matrix_t<Number> & matrix, const vector<vecNumber>& vec){
      std::vector<vecNumber> Num(matrix.cols());
			std::iota(Num.begin(), Num.end(), 0);
			std::vector<vecNumber> newCol = setdiff(Num, vec);
      Matrix_t<Number> newMatrix(matrix.rows(), newCol.size());

      for(Eigen::Index i = 0; i < newCol.size(); i++){
        newMatrix.col(i) = matrix.col(newCol[i]);
      }

      matrix = newMatrix;
    }
    template <typename Number, typename vecNumber>
	  void RemoveRow_Vec(Matrix_t<Number> & matrix, const vector<vecNumber>& vec){
      std::vector<vecNumber> Num(matrix.rows());
			std::iota(Num.begin(), Num.end(), 0);
			std::vector<vecNumber> newRow = setdiff(Num, vec);
      Matrix_t<Number> newMatrix(newRow.size(), matrix.cols());

      for(Eigen::Index i = 0; i < newRow.size(); i++){
        newMatrix.row(i) = matrix.row(newRow[i]);
      }

      matrix = newMatrix;
    }

  template <typename Number, typename vecNumber>
	void RemoveColumn(Matrix_t<Number> & matrix, vecNumber colToRemove) {
 		vecNumber numRows = matrix.rows();
  	vecNumber numCols = matrix.cols() - 1;
 
  		if( colToRemove < numCols ) {
    		matrix.block(0, colToRemove, numRows, numCols - colToRemove) =
    	  	matrix.block(0, colToRemove + 1, numRows, numCols - colToRemove);
  		}
	
  		matrix.conservativeResize(numRows,numCols);
	}

	template <typename Number, typename vecNumber>
	void RemoveRow(Matrix_t<Number> & matrix, vecNumber rowToRemove) {
		vecNumber numRows = matrix.rows() - 1;
		vecNumber numCols = matrix.cols();
		
		if( rowToRemove < numRows ) {
		matrix.block(rowToRemove,0,numRows-rowToRemove,numCols) =
			matrix.block(rowToRemove+1,0,numRows-rowToRemove,numCols);
		}
		
		matrix.conservativeResize(numRows,numCols);
	}

    // template <typename Number, typename vecNumber>
	  // void RemoveColumn_Vec(Matrix_t<Number> & matrix, vector<vecNumber>& vec){

		//   for (typename std::vector<vecNumber>::reverse_iterator r_it = vec.rbegin(); r_it != vec.rend(); ++r_it)
		//   {
		// 	  RemoveColumn(matrix, *r_it);
		//   }
    // }
    // template <typename Number, typename vecNumber>
	  // void RemoveRow_Vec(Matrix_t<Number> & matrix, vector<vecNumber>& vec){

		//   for (typename std::vector<vecNumber>::reverse_iterator r_it = vec.rbegin(); r_it != vec.rend(); ++r_it)
		//   {
		// 	  RemoveRow(matrix, *r_it);
		//   }
    // }
  template <typename Number>
  void comb(int n, vector<Number>& picked, int toPick, vector<vector<Number>>& com){
      if (toPick == 0) { com.push_back(picked); return; }
      int smallest = picked.empty() ? 0 : picked.back() + 1;
      for (int next = smallest; next < n; ++next){
          picked.push_back(next);
          comb(n, picked, toPick-1, com);
          picked.pop_back();
                                  // 关键！！！
      }
  }

  template <typename Number>
	bool isIncludePlanes(const vector<Number>& comb, const vector<vector<Number>>& planes){
		for(int i = 0; i < planes.size(); i++){
			int comb_index = 0;
			int planes_index = 0;
			while(comb_index < comb.size() && planes_index < planes[i].size()){
				if(comb[comb_index] == planes[i][planes_index]){
					comb_index ++;
					planes_index ++;
				}else if(comb[comb_index] > planes[i][planes_index]){
					planes_index ++;
				}else{
					break;
				}
			}
			if(comb_index == comb.size()){
				return 1;
			}
		}
		return 0;
	}

  template <typename Number>
  void matrix2glpkForm(const Matrix_t<Number>& matrix, int* ia, int* ja, Number* ar){

    int rows = matrix.rows();
    int cols = matrix.cols();
    
    for(int i = 0; i < rows; i++){
      for(int j = 0; j < cols; j++){
        ia[1 + i*cols + j] = i + 1;
        ja[1 + i*cols + j] = j + 1;
        ar[1 + i*cols + j] = matrix(i, j);
      }
    }
  }

  template <typename Number>
  double vector2NormSquare(const Vector_t<Number>& vector){

    double res = 0;

    for(int i = 0; i < vector.rows(); i++){
      res += vector(i) * vector(i);
    }

    return res;
  }

  /**
 * @brief: sort generators of a zonotope generator matrix in descending order.
 */
  template <typename Number>
  Matrix_t<Number> sortGenDes(const Matrix_t<Number>& A) {
    
    Matrix_t<Number> res(A.rows(), A.cols());

    std::vector<int> idx(A.cols());
		std::iota(idx.begin(), idx.end(), 0);

    std::vector<double> norm2(A.cols(), 0);

    for(int i = 0; i < norm2.size(); i++){
      for(int j = 0; j < A.rows(); j++){
        norm2[i] += A(j,i) * A(j,i);
      }
    }
    
		std::sort(idx.begin(), idx.end(), [&norm2](int i1, int i2) { return norm2[i1] > norm2[i2]; });
		
    for(int i = 0; i < A.cols(); i++){
      res.col(i) = A.col(idx[i]);
    }

    return res;
  }

  /**
 * @brief: 根据与目标 generator cos 值从小到大对 zonotope 的 generator 进行排序
 */
  template <typename Number>
  Matrix_t<Number> sortGenCos(const Matrix_t<Number>& A, const Vector_t<Number>& g, vector<int>& idx) {
    
    Matrix_t<Number> res(A.rows(), A.cols());

    // std::vector<int> idx(A.cols());
		std::iota(idx.begin(), idx.end(), 0);

    std::vector<double> norm2(A.cols(), 0);
    double gnorm2 = 0;

    std::vector<double> cos(A.cols());

    for(int i = 0; i < g.cols(); i++){
      gnorm2 += g(i) * g(i);
    }

    gnorm2 = sqrt(gnorm2);

    for(int i = 0; i < norm2.size(); i++){
      for(int j = 0; j < A.rows(); j++){
        norm2[i] += A(j,i) * A(j,i);
      }
      norm2[i] = sqrt(norm2[i]);
    }

    for(int i = 0; i < A.cols(); i++){

      cos[i] = abs(A.col(i).dot(g))/(gnorm2 * norm2[i]);
    }
    
		std::sort(idx.begin(), idx.end(), [&cos](int i1, int i2) { return cos[i1] > cos[i2]; });
		
    for(int i = 0; i < A.cols(); i++){
      res.col(i) = A.col(idx[i]);
    }

    return res;
  }
} // namespace reachSolver
#endif // POLYZONOTOPE_H