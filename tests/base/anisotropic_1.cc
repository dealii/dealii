// ---------------------------------------------------------------------
//
// Copyright (C) 2003 - 2013 by the deal.II authors
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



// check that TensorProducPolynomials and AnisotropicPolynomials
// compute the same thing if the basis polynomials for the anisotropic
// class happen to be the same for all coordinate directions


#include "../tests.h"
#include <iomanip>
#include <fstream>
#include <cmath>

#include <deal.II/base/logstream.h>
#include <deal.II/base/tensor_product_polynomials.h>


using namespace Polynomials;


template<int dim, class POLY1, class POLY2>
void check_poly(const Point<dim> &x,
                const POLY1      &p,
                const POLY2      &q)
{
  const unsigned int n = p.n();

  std::vector<double> values1 (n), values2 (n);
  std::vector<Tensor<1,dim> > gradients1(n), gradients2(n);
  std::vector<Tensor<2,dim> > second1(n), second2(n);

  p.compute (x, values1, gradients1, second1);
  q.compute (x, values2, gradients2, second2);

  for (unsigned int k=0; k<n; ++k)
    {
      // Check if compute_value is ok
      double val1 = p.compute_value(k,x);
      if (std::fabs(val1 - values1[k]) > 5.0E-16)
        deallog << 'P' << k << ": values differ " << val1 << " != "
                << values1[k] << std::endl;
      double val2 = q.compute_value(k,x);
      if (std::fabs(val2 - values2[k]) > 5.0E-16)
        deallog << 'Q' << k << ": values differ " << val2 << " != "
                << values2[k] << std::endl;
      if (std::fabs(val2 - val1) > 5.0E-16)
        deallog << "PQ" << k << ": values differ " << val1 << " != "
                << val2 << std::endl;

      // Check if compute_grad is ok
      Tensor<1,dim> grad1 = p.compute_grad(k,x);
      if (grad1 != gradients1[k])
        deallog << 'P' << k << ": gradients differ " << grad1 << " != "
                << gradients1[k] << std::endl;
      Tensor<1,dim> grad2 = q.compute_grad(k,x);
      if (grad2 != gradients2[k])
        deallog << 'Q' << k << ": gradients differ " << grad1 << " != "
                << gradients2[k] << std::endl;
      if (grad2 != grad1)
        deallog << "PQ" << k << ": gradients differ " << grad1 << " != "
                << grad2 << std::endl;

      // Check if compute_grad_grad is ok
      Tensor<2,dim> grad_grad1 = p.compute_grad_grad(k,x);
      if (grad_grad1 != second1[k])
        deallog << 'P' << k << ": second derivatives differ " << grad_grad1 << " != "
                << second1[k] << std::endl;
      Tensor<2,dim> grad_grad2 = q.compute_grad_grad(k,x);
      if (grad_grad2 != second2[k])
        deallog << 'Q' << k << ": second derivatives differ " << grad_grad2 << " != "
                << second2[k] << std::endl;
      if (grad_grad2 != grad_grad1)
        deallog << "PQ" << k << ": second derivatives differ " << grad_grad1 << " != "
                << grad_grad2 << std::endl;


      // finally output values,
      // gradients, etc, to make sure
      // that they are not only
      // consistent, but also
      // correct. Multiply them
      // somewhat to make them
      // significant despite our
      // two-post-dot-digits limit
      values1[k] *= std::pow(10., dim);
      gradients1[k] *= std::pow(10., dim);

      deallog << 'P' << k << "\t= " << values1[k]
              << "\tgradient\t";
      for (unsigned int d=0; d<dim; ++d)
        deallog << gradients1[k][d] << '\t';
      deallog << "\t2nd\t";
      for (unsigned int d1=0; d1<dim; ++d1)
        for (unsigned int d2=0; d2<dim; ++d2)
          deallog << second1[k][d1][d2] << '\t';
      deallog << std::endl;
    }
  deallog << std::endl;
}


template <int dim>
void
check_tensor (const std::vector<Polynomial<double> > &v,
              const Point<dim> &x)
{
  TensorProductPolynomials<dim> p(v);

  std::vector<std::vector<Polynomial<double> > > pols (dim, v);
  AnisotropicPolynomials<dim> q(pols);

  check_poly (x, p, q);
}


void
check_dimensions (const std::vector<Polynomial<double> > &p)
{
  deallog.push("1d");
  check_tensor(p, Point<1>(.5));
  deallog.pop();

  deallog.push("2d");
  check_tensor(p, Point<2>(.5, .2));
  deallog.pop();

  deallog.push("3d");
  check_tensor(p, Point<3>(.5, .2, .3));
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
