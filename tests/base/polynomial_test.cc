//----------------------------  polynomial_test.cc  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2000, 2001, 2002, 2003 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  polynomial_test.cc  ---------------------------


// just output a lot of information about various classes implementing
// polynomials, to make sure that all changes we make to these classes
// do not change the results of these classes.


#include <iostream>
#include <fstream>
#include <cmath>

#include <base/logstream.h>
#include <base/tensor_product_polynomials.h>
#include <base/polynomial_space.h>


using namespace Polynomials;


template<int dim, class POLY>
void check_poly(const Point<dim>& x,
		const POLY& p)
{
  const unsigned int n = p.n();
  std::vector<double> values (n);
  std::vector<Tensor<1,dim> > gradients(n);
  std::vector<Tensor<2,dim> > second(n);
  
  p.compute (x, values, gradients, second);
  
  for (unsigned int k=0;k<n;++k)
    {
				       // Check if compute_value is ok
      double val = p.compute_value(k,x);
      if (val != values[k])
	deallog << 'P' << k << ": values differ " << val << " != "
		<< values[k] << std::endl;

				       // Check if compute_grad is ok
      Tensor<1,dim> grad = p.compute_grad(k,x);
      if (grad != gradients[k])
	deallog << 'P' << k << ": gradients differ " << grad << " != "
		<< gradients[k] << std::endl;
      
				       // Check if compute_grad_grad is ok
      Tensor<2,dim> grad2 = p.compute_grad_grad(k,x);
      if (grad2 != second[k])
	deallog << 'P' << k << ": second derivatives differ " << grad2 << " != "
		<< second[k] << std::endl;
      
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
      deallog << std::endl;
    }
  deallog << std::endl;
}


template <int dim>
void
check_tensor (const std::vector<Polynomial<double> >& v,
	      const Point<dim>& x)
{
  deallog.push("Tensor");
  TensorProductPolynomials<dim> p(v);
  check_poly (x, p);
  deallog.pop();
}


template <int dim>
void
check_poly (const std::vector<Polynomial<double> >& v,
	    const Point<dim>& x)
{
  deallog.push("Polyno");
  PolynomialSpace<dim> p(v);
  check_poly (x, p);
  deallog.pop();
}


void
check_dimensions (const std::vector<Polynomial<double> >& p)
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

  deallog.push("Lagrange");
  std::vector<Polynomial<double> > p;
  for (unsigned int i=0;i<3;++i)
    p.push_back (LagrangeEquidistant(3, i));

  check_dimensions(p);

  deallog.pop();
  deallog.push("Legendre");
  
  p.clear ();
  for (unsigned int i=0;i<3;++i)
    p.push_back (Legendre<double>(i));

  check_dimensions(p);

  deallog.pop();
}
