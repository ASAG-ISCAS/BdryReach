//////////////////////////////////////////////////////////////////////////////
//   Package:          CAPD

/////////////////////////////////////////////////////////////////////////////
//
/// @file sqrTest.cpp
///
/// @author Daniel Wilczak
//
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) CAPD group
//
// This file constitutes a part of the CAPD library,
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.uj.edu.pl/ for details.

//#define BOOST_TEST_MODULE sqrTest
#include "compare.h"
BOOST_AUTO_TEST_SUITE(sqrSuite)

std::vector<double> computeSqrDer(MapType::VectorType & u){
  double x = u[0].leftBound();
  double y = u[1].leftBound();

  // code generated by the following Mathematica code
  // W[n_,m_]:=D[(Sin[x]*Cos[y])^2,{x,n},{y,m}]/(n!m!)//FullSimplify
  // Table[Table[W[m-n,n]//CForm,{n,0,m}],{m,0,5}]//Flatten
  capd::rounding::DoubleRounding::roundNearest();
  double r[] = {Power(Cos(y),2)*Power(Sin(x),2),Power(Cos(y),2)*Sin(2*x),-2*Cos(y)*Power(Sin(x),2)*Sin(y),Cos(2*x)*Power(Cos(y),2),-4*Cos(x)*Cos(y)*Sin(x)*Sin(y),-(Cos(2*y)*Power(Sin(x),2)),(-4*Cos(x)*Power(Cos(y),2)*Sin(x))/3.,-(Cos(2*x)*Sin(2*y)),-(Cos(2*y)*Sin(2*x)),(4*Cos(y)*Power(Sin(x),2)*Sin(y))/3.,-(Cos(2*x)*Power(Cos(y),2))/3.,(8*Cos(x)*Cos(y)*Sin(x)*Sin(y))/3.,-(Cos(2*x)*Cos(2*y)),(8*Cos(x)*Cos(y)*Sin(x)*Sin(y))/3.,(Cos(2*y)*Power(Sin(x),2))/3.,(4*Cos(x)*Power(Cos(y),2)*Sin(x))/15.,(Cos(2*x)*Sin(2*y))/3.,(2*Cos(2*y)*Sin(2*x))/3.,(2*Cos(2*x)*Sin(2*y))/3.,(Cos(2*y)*Sin(2*x))/3.,(-4*Cos(y)*Power(Sin(x),2)*Sin(y))/15.};
  return std::vector<double> (r,r+sizeof(r)/sizeof(double));
}

BOOST_AUTO_TEST_CASE(xsqr)
{
  std::string txt = "var:x,y;fun:sqr(sin(x)*cos(y));",
              msg = "Function \"" + txt + "\"  x = " ;
  MapType f(txt,5);
  VectorType x(2);
  JetType df(1,2,5);

  x[0] = .5; x[1] = 0.25;
  std::vector<double> expected = computeSqrDer(x);
  f(x,df);
  compareResults(expected, df, msg+"(0.5,0.25)");

  MapType g("var:x,y;fun:sqr(-sin(x)*cos(-y));",5);
  g(x,df);
  compareResults(expected, df, msg+"(0.5,0.25)");

  x[0] = -0.75; x[1] = -0.1;
  expected = computeSqrDer(x);
  f(x,df);
  compareResults(expected, df, msg+"(-0.75,-0.1)");

  x[0] = 0.0; x[1]= 0.0;
  expected = computeSqrDer(x);
  f(x,df);
  compareResults(expected, df, msg+"(0.0,0.0)");

  x[0] = 1.0; x[1] = 0.0;
  expected = computeSqrDer(x);
  f(x,df);
  compareResults(expected, df, msg+"(1.0,0.0)");

  x[0] = 0.0; x[1] = 1.0;
  expected = computeSqrDer(x);
  f(x,df);
  compareResults(expected, df, msg+"(0.0,1.0)");
}

using capd::autodiff::Node;

void _f(Node /*t*/, Node in[], int /*dimIn*/, Node out[], int /*dimOut*/, Node /*params*/[], int /*noParams*/)
{
  out[0] = sqr(sin(in[0])*cos(in[1]));
}

BOOST_AUTO_TEST_CASE(xsqrnode)
{
  std::string msg = "Function \"sqr(sin(x)*cos(y))\"  (x,y) = " ;
  MapType f(_f,2,1,0,5);
  VectorType x(2);
  JetType df(1,2,5);

  x[0] = .1; x[1]=1;
  std::vector<double> expected = computeSqrDer(x);
  f(x,df);
  compareResults(expected, df, msg+"(.1,1)");

  x[0] = -0.143; x[1] = 0.6;
  expected = computeSqrDer(x);
  f(x,df);
  compareResults(expected, df, msg+"(-0.143,0.6)");

  x[0] = 0.0; x[1] = 0.0;
  expected = computeSqrDer(x);
  f(x,df);
  compareResults(expected, df, msg+"(0.0,0.0)");

  x[0] = 1.0; x[1] = 0.0;
  expected = computeSqrDer(x);
  f(x,df);
  compareResults(expected, df, msg+"(1.0,0.0)");

  x[0] = 0.0; x[1] = 1.0;
  expected = computeSqrDer(x);
  f(x,df);
  compareResults(expected, df, msg+"(0.0,1.0)");
}
BOOST_AUTO_TEST_SUITE_END()