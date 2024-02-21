// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2000 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// just output a lot of information about various classes implementing
// polynomials, to make sure that all changes we make to these classes
// do not change the results of these classes.


#include <deal.II/base/polynomial_space.h>
#include <deal.II/base/tensor_product_polynomials.h>

#include "../tests.h"


using namespace Polynomials;


template <int dim, class PolynomialType>
void
check_poly(const Point<dim> &x, const PolynomialType &p)
{
  const unsigned int          n   = p.n();
  const double                eps = 1.0e-14;
  std::vector<double>         values(n);
  std::vector<Tensor<1, dim>> gradients(n);
  std::vector<Tensor<2, dim>> second(n);
  std::vector<Tensor<3, dim>> third(n);
  std::vector<Tensor<4, dim>> fourth(n);

  p.evaluate(x, values, gradients, second, third, fourth);

  for (unsigned int k = 0; k < n; ++k)
    {
      // first make sure the
      // individual functions work in
      // a consistent way

      // Check if compute_value is ok
      double val = p.compute_value(k, x);
      if (std::fabs(val - values[k]) > eps)
        deallog << 'P' << k << ": values differ " << val << " != " << values[k]
                << std::endl;

      // Check if compute_grad is ok
      Tensor<1, dim> grad = p.template compute_derivative<1>(k, x);
      if ((grad - gradients[k]) * (grad - gradients[k]) > eps * eps)
        deallog << 'P' << k << ": gradients differ " << grad
                << " != " << gradients[k] << std::endl;

      // Check if compute_grad_grad is ok
      Tensor<2, dim> grad2 = p.template compute_derivative<2>(k, x);
      Tensor<2, dim> diff  = grad2 - second[k];

      if (diff.norm_square() > eps * eps)
        deallog << 'P' << k << ": second derivatives differ " << grad2
                << " != " << second[k] << std::endl;

      // Check if third derivative is ok
      Tensor<3, dim> grad3 = p.template compute_derivative<3>(k, x);
      if ((grad3 - third[k]).norm_square() > 5e-15 * 5e-15)
        deallog << 'P' << k << ": third derivatives differ " << grad3
                << " != " << third[k] << std::endl;

      if (diff.norm_square() > eps * eps)
        deallog << 'P' << k << ": second derivatives differ " << grad2
                << " != " << second[k] << std::endl;

      // Check if third derivative is ok
      Tensor<4, dim> grad4 = p.template compute_derivative<4>(k, x);
      if ((grad3 - third[k]).norm_square() > eps * eps)
        deallog << 'P' << k << ": fourth derivatives differ " << grad4
                << " != " << fourth[k] << std::endl;


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

      deallog << 'P' << k << "\t= " << values[k] << "\tgradient\t";
      for (unsigned int d = 0; d < dim; ++d)
        deallog << gradients[k][d] << '\t';
      deallog << "\t2nd\t";
      for (unsigned int d1 = 0; d1 < dim; ++d1)
        for (unsigned int d2 = 0; d2 < dim; ++d2)
          deallog << second[k][d1][d2] << '\t';
      deallog << "\t3rd\t";
      for (unsigned int d1 = 0; d1 < dim; ++d1)
        for (unsigned int d2 = 0; d2 < dim; ++d2)
          for (unsigned int d3 = 0; d3 < dim; ++d3)
            deallog << third[k][d1][d2][d3] << '\t';
      deallog << "\t4th\t";
      for (unsigned int d1 = 0; d1 < dim; ++d1)
        for (unsigned int d2 = 0; d2 < dim; ++d2)
          for (unsigned int d3 = 0; d3 < dim; ++d3)
            for (unsigned int d4 = 0; d4 < dim; ++d4)
              deallog << fourth[k][d1][d2][d3][d4] << '\t';
      deallog << std::endl;
    }
  deallog << std::endl;
}


template <int dim>
void
check_tensor(const std::vector<Polynomial<double>> &v, const Point<dim> &x)
{
  deallog.push("Tensor");
  TensorProductPolynomials<dim> p(v);
  check_poly(x, p);
  deallog.pop();
}


template <int dim>
void
check_poly(const std::vector<Polynomial<double>> &v, const Point<dim> &x)
{
  deallog.push("Polyno");
  PolynomialSpace<dim> p(v);
  p.output_indices(deallog);
  check_poly(x, p);
  deallog.pop();
}


void
check_dimensions(const std::vector<Polynomial<double>> &p)
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

int
main()
{
  initlog();
  deallog << std::setprecision(2);

  deallog.push("Lagrange");
  std::vector<Polynomial<double>> p;
  for (unsigned int i = 0; i < 3; ++i)
    p.push_back(LagrangeEquidistant(3, i));

  check_dimensions(p);

  deallog.pop();
  deallog.push("Legendre");

  p.clear();
  for (unsigned int i = 0; i < 3; ++i)
    p.push_back(Legendre(i));

  check_dimensions(p);

  deallog.pop();
  deallog.push("Hierarchical");

  p.clear();
  for (unsigned int i = 0; i < 3; ++i)
    p.push_back(Hierarchical(i));

  check_dimensions(p);

  deallog.pop();
}
