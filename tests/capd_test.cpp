#include <iostream>
#include "capd/capdlib.h"
using namespace capd;
using namespace std;
using capd::autodiff::Node;
// ####################################################
/*
 * This is a map we will evaluate and differentiate (dimIn = 3, dimOut = 2, noParams = 5)
 * @param in is an array of independent variables
 * @param out is an array of dependent variables, i.e. out = f(in)
 * @param params - parameters of the map.
 */

double k  = 0.015;
double k2 = 0.01;
double g  = 9.81;

void _f(Node/* t*/, Node in[], int /*dimIn*/, Node out[], int/* dimOut*/, Node params[], int noParams){
	out[0] = k * sqrt(2 * g) * (-sqrt(in[0]));
}
// ####################################################
int main(){
int dimIn			   = 1;
int dimOut			   = 1;
int noParam			   = 0;
int MaxDerivativeOrder = 3;
  IMap f(_f,dimIn,dimOut,noParam,MaxDerivativeOrder);
  // set parameter values

  IVector x(dimIn);
  x[0].setLeftBound(3.79304);
	x[0].setRightBound(4.20672);
  // declare an object for storing jets
  IJet jet(dimOut,dimIn,MaxDerivativeOrder);
  f(x,jet);
  // print the result
  // NOTE: indices of variables must form a nondecreasing sequence!
  cout << "df0_dx0_dx0_dx0 = " << jet(0,0,0,0) << endl;
  cout << jet.toString();
}
