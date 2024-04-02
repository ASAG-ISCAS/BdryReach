//////////////////////////////////////////////////////////////////////////////
//   Package:          CAPD

/////////////////////////////////////////////////////////////////////////////
//
/// @file cuubeTest.cpp
///
/// @author Daniel Wilczak
//
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) CAPD group
//
// This file constitutes a part of the CAPD library,
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.uj.edu.pl/ for details.

//#define BOOST_TEST_MODULE cubeTest
#include "compare.h"
BOOST_AUTO_TEST_SUITE(cubeSuite)

std::vector<double> computeCubeDer(MapType::VectorType & u){
  double x = u[0].leftBound();
  double y = u[1].leftBound();

  // code generated by the following Mathematica code
  // W[n_,m_]:=D[(Sin[x]*Cos[y])^3,{x,n},{y,m}]/(n!m!)//FullSimplify
  // Table[Table[W[m-n,n]//CForm,{n,0,m}],{m,0,5}]//Flatten
  capd::rounding::DoubleRounding::roundNearest();
  double r[] = {Power(Cos(y),3)*Power(Sin(x),3),3*Cos(x)*Power(Cos(y),3)*Power(Sin(x),2),-3*Power(Cos(y),2)*Power(Sin(x),3)*Sin(y),(-3*Power(Cos(y),3)*(Sin(x) - 3*Sin(3*x)))/8.,-9*Cos(x)*Power(Cos(y),2)*Power(Sin(x),2)*Sin(y),(-3*(Cos(y) + 3*Cos(3*y))*Power(Sin(x),3))/8.,-((Cos(x) - 9*Cos(3*x))*Power(Cos(y),3))/8.,(9*Power(Cos(y),2)*(Sin(x) - 3*Sin(3*x))*Sin(y))/8.,(-9*Cos(x)*(Cos(y) + 3*Cos(3*y))*Power(Sin(x),2))/8.,(Power(Sin(x),3)*(Sin(y) + 9*Sin(3*y)))/8.,(Power(Cos(y),3)*(Sin(x) - 27*Sin(3*x)))/32.,(3*(Cos(x) - 9*Cos(3*x))*Power(Cos(y),2)*Sin(y))/8.,(9*(Cos(y) + 3*Cos(3*y))*(Sin(x) - 3*Sin(3*x)))/64.,(3*Cos(x)*Power(Sin(x),2)*(Sin(y) + 9*Sin(3*y)))/8.,((Cos(y) + 27*Cos(3*y))*Power(Sin(x),3))/32.,((Cos(x) - 81*Cos(3*x))*Power(Cos(y),3))/160.,(-3*Power(Cos(y),2)*(Sin(x) - 27*Sin(3*x))*Sin(y))/32.,(3*(Cos(x) - 9*Cos(3*x))*(Cos(y) + 3*Cos(3*y)))/64.,(-3*(Sin(x) - 3*Sin(3*x))*(Sin(y) + 9*Sin(3*y)))/64.,(3*Cos(x)*(Cos(y) + 27*Cos(3*y))*Power(Sin(x),2))/32.,-(Power(Sin(x),3)*(Sin(y) + 81*Sin(3*y)))/160.};
  return std::vector<double> (r,r+sizeof(r)/sizeof(double));
}

BOOST_AUTO_TEST_CASE(xcube)
{
  std::string txt = "var:x,y;fun:(sin(x)*cos(y))^3;",
              msg = "Function \"" + txt + "\"  x = " ;
  MapType f(txt,5);
  VectorType x(2);
  JetType df(1,2,5);

  x[0] = .5; x[1] = 0.25;
  std::vector<double> expected = computeCubeDer(x);
  f(x,df);
  compareResults(expected, df, msg+"(0.5,0.25)");

  MapType g("var:x,y;fun:-(sin(-x)*cos(-y))^3;",5);
  g(x,df);
  compareResults(expected, df, msg+"(0.5,0.25)");

  x[0] = -0.75; x[1] = -0.1;
  expected = computeCubeDer(x);
  f(x,df);
  compareResults(expected, df, msg+"(-0.75,-0.1)");

  x[0] = 0.0; x[1]= 0.0;
  expected = computeCubeDer(x);
  f(x,df);
  compareResults(expected, df, msg+"(0.0,0.0)");

  x[0] = 1.0; x[1] = 0.0;
  expected = computeCubeDer(x);
  f(x,df);
  compareResults(expected, df, msg+"(1.0,0.0)");

  x[0] = 0.0; x[1] = 1.0;
  expected = computeCubeDer(x);
  f(x,df);
  compareResults(expected, df, msg+"(0.0,1.0)");
}

using capd::autodiff::Node;

void _f(Node /*t*/, Node in[], int /*dimIn*/, Node out[], int /*dimOut*/, Node /*params*/[], int /*noParams*/)
{
  out[0] = (sin(in[0])*cos(in[1]))^3;
}

BOOST_AUTO_TEST_CASE(xcubenode)
{
  std::string msg = "Function \"(sin(x)*cos(y))^3\"  (x,y) = " ;
  MapType f(_f,2,1,0,5);
  VectorType x(2);
  JetType df(1,2,5);

  x[0] = .125; x[1]=.75;
  std::vector<double> expected = computeCubeDer(x);
  f(x,df);
  compareResults(expected, df, msg+"(.125,.75)");
  return;
  x[0] = -0.125; x[1] = 0.6;
  expected = computeCubeDer(x);
  f(x,df);
  compareResults(expected, df, msg+"(-0.125,0.6)");

  x[0] = 0.0; x[1] = 0.0;
  expected = computeCubeDer(x);
  f(x,df);
  compareResults(expected, df, msg+"(0.0,0.0)");

  x[0] = 1.0; x[1] = 0.0;
  expected = computeCubeDer(x);
  f(x,df);
  compareResults(expected, df, msg+"(1.0,0.0)");

  x[0] = 0.0; x[1] = 1.0;
  expected = computeCubeDer(x);
  f(x,df);
  compareResults(expected, df, msg+"(0.0,1.0)");
}
BOOST_AUTO_TEST_SUITE_END()