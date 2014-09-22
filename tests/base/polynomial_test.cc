// ---------------------------------------------------------------------
//
// Copyright (C) 2000 - 2013 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------



// just output a lot of information about various classes implementing
// polynomials, to make sure that all changes we make to these classes
// do not change the results of these classes.


#include "../tests.h"
#include <iomanip>
#include <fstream>
#include <cmath>

#include <deal.II/base/logstream.h>
#include <deal.II/base/tensor_product_polynomials.h>
#include <deal.II/base/polynomial_space.h>


using namespace Polynomials;


template<int dim, class POLY>
void check_poly(const Point<dim> &x,
                const POLY &p)
{
  const unsigned int n = p.n();
  std::vector<double> values (n);
  std::vector<Tensor<1,dim> > gradients(n);
  std::vector<Tensor<2,dim> > second(n);

  p.compute (x, values, gradients, second);

  for (unsigned int k=0; k<n; ++k)
    {
      // first make sure the
      // individual functions work in
      // a consistent way

      // Check if compute_value is ok
      double val = p.compute_value(k,x);
      if (std::fabs(val - values[k]) > 5.0E-15)
        deallog << 'P' << k << ": values differ " << val << " != "
                << values[k] << std::endl;

      // Check if compute_grad is ok
      Tensor<1,dim> grad = p.compute_grad(k,x);
      if ((grad-gradients[k]) * (grad-gradients[k]) > 5e-15*5e-15)
        deallog << 'P' << k << ": gradients differ " << grad << " != "
                << gradients[k] << std::endl;

      // Check if compute_grad_grad is ok
      Tensor<2,dim> grad2 = p.compute_grad_grad(k,x);
      Tensor<2,dim> diff = grad2-second[k];
      double s = 0;
      for (unsigned int i=0; i<dim; ++i)
        for (unsigned int j=0; j<dim; ++j)
          s += diff[i][j]*diff[i][j];
      if (s > 5e-15*5e-15)
        deallog << 'P' << k << ": second derivatives differ " << grad2 << " != "
                << second[k] << std::endl;


      // finally output values,
      // gradients, etc, to make sure
      // that they are not only
      // consistent, but also
      // correct. Multiply them
      // somewhat to make them
      // significant despite our
      // two-post-dot-digits limit
      values[k] *= std::pow(10., dim);
      gradients[k] *= std::pow(10., dim);

      deallog << 'P' << k << "\t= " << values[k]
              << "\tgradient\t";
      for (unsigned int d=0; d<dim; ++d)
        deallog << gradients[k][d] << '\t';
      deallog << "\t2nd\t";
      for (unsigned int d1=0; d1<dim; ++d1)
        for (unsigned int d2=0; d2<dim; ++d2)
          deallog << second[k][d1][d2] << '\t';
      deallog << std::endl;
    }
  deallog << std::endl;
}


template <int dim>
void
check_tensor (const std::vector<Polynomial<double> > &v,
              const Point<dim> &x)
{
  deallog.push("Tensor");
  TensorProductPolynomials<dim> p(v);
  check_poly (x, p);
  deallog.pop();
}


template <int dim>
void
check_poly (const std::vector<Polynomial<double> > &v,
            const Point<dim> &x)
{
  deallog.push("Polyno");
  PolynomialSpace<dim> p(v);
  p.output_indices(deallog);
  check_poly (x, p);
  deallog.pop();
}


void
check_dimensions (const std::vector<Polynomial<double> > &p)
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
  std::ofstream logfile("output");
  deallog << std::setprecision(2);
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  deallog.push("Lagrange");
  std::vector<Polynomial<double> > p;
  for (unsigned int i=0; i<3; ++i)
    p.push_back (LagrangeEquidistant(3, i));

  check_dimensions(p);

  deallog.pop();
  deallog.push("Legendre");

  p.clear ();
  for (unsigned int i=0; i<3; ++i)
    p.push_back (Legendre(i));

  check_dimensions(p);

  deallog.pop();
  deallog.push("Hierarchical");

  p.clear ();
  for (unsigned int i=0; i<3; ++i)
    p.push_back (Hierarchical(i));

  check_dimensions(p);

  deallog.pop();
}
