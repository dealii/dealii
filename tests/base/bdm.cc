//-----------------------------------------------------------------------------
//    bdm.cc,v 1.1 2004/01/04 19:07:59 guido Exp
//    Version: 
//
//    Copyright (C) 2000, 2001, 2002, 2003 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//-----------------------------------------------------------------------------

#include "../tests.h"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>

#include <base/logstream.h>
#include <base/table.h>
#include <base/polynomial.h>
#include <base/polynomials_bdm.h>
#include <base/quadrature_lib.h>

using namespace Polynomials;

template <int dim>
void moments (PolynomialsBDM<dim>& poly)
{
  std::vector<Polynomial<double> > legendre(2);
  for (unsigned int i=0;i<legendre.size();++i)
    legendre[i] = Monomial<double>(i);

  QGauss<1> qface(2);

  Table<2,double> integrals (poly.n(), poly.n());

  std::vector<Tensor<1,dim> > values(poly.n());
  std::vector<Tensor<2,dim> > grads;
  std::vector<Tensor<3,dim> > grad_grads;
  values.resize(poly.n());

  for (unsigned int face=0;face<2*dim;++face)
    for (unsigned int k=0;k<qface.n_quadrature_points;++k)
      {
	const double w = qface.weight(k);
	const double x = qface.point(k)(0);
	Point<dim> p;
	switch (face)
	  {
	    case 2:
	      p(1) = 1.;
	    case 0:
	      p(0) = x;
	      break;
	    case 1:
	      p(0) = 1.;
	    case 3:
	      p(1) = x;
	      break;	      
	  }
	std::cerr << p << std::endl;
	
	poly.compute (p, values, grads, grad_grads);
	for (unsigned int i=0;i<poly.n();++i)
	  {
	    integrals(2*face,i) += w * values[i][1-face%2] ;
	    integrals(2*face+1,i) += w * values[i][1-face%2] * legendre[1].value(x);
	  }
      }

  std::cout.setf(std::ios::fixed);
  for (unsigned int i=0;i<integrals.n_rows();++i)
    {
      for (unsigned int j=0;j<poly.n();++j)
	std::cout << std::setw(13) << integrals(i,j);
      std::cout << std::endl;
    }
}



int main()
{
  PolynomialsBDM<2> bdm1(1);

  moments(bdm1);
}
