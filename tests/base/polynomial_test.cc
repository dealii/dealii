//----------------------------  polynomial_test.cc  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2000, 2001, 2002 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  polynomial_test.cc  ---------------------------

#include <iostream>
#include <fstream>
#include <cmath>

#include <base/logstream.h>
#include <base/tensor_product_polynomials.h>
#include <base/polynomial_space.h>

//using std;

extern "C"
void abort()
{}


template<int dim, class POLY>
void check_poly(const Point<dim>& x,
		const POLY& p)
{
  const unsigned int n = p.n();
  vector<double> values (n);
  vector<Tensor<1,dim> > gradients(n);
  vector<Tensor<2,dim> > second(n);
  
  p.compute (x, values, gradients, second);
  
  for (unsigned int k=0;k<n;++k)
    {
      values[k] *= pow(10, dim);
      gradients[k] *= pow(10, dim);
      
      deallog << 'P' << k << "\t= " << values[k]
	      << "\tgradient\t";
      for (unsigned int d=0;d<dim;++d)
	deallog << gradients[k][d] << '\t';
      deallog << "\t2nd\t";
      for (unsigned int d1=0;d1<dim;++d1)
	for (unsigned int d2=0;d2<dim;++d2)
	  deallog << second[k][d1][d2] << '\t';
      deallog << endl;
    }
  deallog << endl;
}


template <int dim>
void
check_tensor (const vector<Polynomial<double> >& v,
	      const Point<dim>& x)
{
  deallog.push("Tensor");
  TensorProductPolynomials<dim> p(v);
  check_poly (x, p);
  deallog.pop();
}


template <int dim>
void
check_poly (const vector<Polynomial<double> >& v,
	    const Point<dim>& x)
{
  deallog.push("Polyno");
  PolynomialSpace<dim> p(v);
  check_poly (x, p);
  deallog.pop();
}


void
check_dimensions (const vector<Polynomial<double> >& p)
{
  deallog.push("1d");
  check_tensor(p, Point<1>(.5));
  check_poly(p, Point<1>(.5));
  deallog.pop();
  deallog.push("2d");
  check_tensor(p, Point<2>(.5, .2));
  check_poly(p, Point<2>(.5, .2));
  deallog.pop();
  deallog.push("3d");
  check_tensor(p, Point<3>(.5, .2, .3));
  check_poly(p, Point<3>(.5, .2, .3));
  deallog.pop();
}

int main()
{
  std::ofstream logfile("polynomial_test.output");
  logfile.precision(2);
  deallog.attach(logfile);
  deallog.depth_console(0);

  vector<Polynomial<double> > p(3);
  for (unsigned int i=0;i<p.size();++i)
    p[i] = LagrangeEquidistant(p.size(), i);

  check_dimensions(p);

  for (unsigned int i=0;i<p.size();++i)
    p[i] = Legendre<double>(i);

  check_dimensions(p);
}
