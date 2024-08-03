// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2003 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// check that TensorProducPolynomials and AnisotropicPolynomials
// compute the same thing if the basis polynomials for the anisotropic
// class happen to be the same for all coordinate directions


#include <deal.II/base/tensor_product_polynomials.h>

#include <deal.II/fe/fe_tools.h>

#include "../tests.h"


using namespace Polynomials;


template <int dim, class PolynomialType1, class PolynomialType2>
void
check_poly(const Point<dim>      &x,
           const PolynomialType1 &p,
           const PolynomialType2 &q)
{
  const unsigned int n = p.n();

  std::vector<double>         values1(n), values2(n);
  std::vector<Tensor<1, dim>> gradients1(n), gradients2(n);
  std::vector<Tensor<2, dim>> second1(n), second2(n);
  std::vector<Tensor<3, dim>> third1(n), third2(n);
  std::vector<Tensor<4, dim>> fourth1(n), fourth2(n);

  p.evaluate(x, values1, gradients1, second1, third1, fourth1);
  q.evaluate(x, values2, gradients2, second2, third2, fourth2);

  for (unsigned int k = 0; k < n; ++k)
    {
      // Check if compute_value is ok
      double val1 = p.compute_value(k, x);
      if (std::fabs(val1 - values1[k]) > 5.0E-12)
        deallog << 'P' << k << ": values differ " << val1
                << " != " << values1[k] << std::endl;
      double val2 = q.compute_value(k, x);
      if (std::fabs(val2 - values2[k]) > 5.0E-12)
        deallog << 'Q' << k << ": values differ " << val2
                << " != " << values2[k] << std::endl;
      if (std::fabs(val2 - val1) > 5.0E-12)
        deallog << "PQ" << k << ": values differ " << val1 << " != " << val2
                << std::endl;

      // Check if compute_grad is ok
      Tensor<1, dim> grad1 = p.template compute_derivative<1>(k, x);
      if ((grad1 - gradients1[k]).norm() > 5.0E-12)
        deallog << 'P' << k << ": gradients differ " << grad1
                << " != " << gradients1[k] << std::endl;
      Tensor<1, dim> grad2 = q.template compute_derivative<1>(k, x);
      if ((grad2 - gradients2[k]).norm() > 5.0E-12)
        deallog << 'Q' << k << ": gradients differ " << grad1
                << " != " << gradients2[k] << std::endl;
      if ((grad2 - grad1).norm() > 5.0E-12)
        deallog << "PQ" << k << ": gradients differ " << grad1
                << " != " << grad2 << std::endl;

      // Check if compute_grad_grad is ok
      Tensor<2, dim> grad_grad1 = p.template compute_derivative<2>(k, x);
      if ((grad_grad1 - second1[k]).norm() > 5.0E-12)
        deallog << 'P' << k << ": second derivatives differ " << grad_grad1
                << " != " << second1[k] << std::endl;
      Tensor<2, dim> grad_grad2 = q.template compute_derivative<2>(k, x);
      if ((grad_grad2 - second2[k]).norm() > 5.0E-12)
        deallog << 'Q' << k << ": second derivatives differ " << grad_grad2
                << " != " << second2[k] << std::endl;
      if ((grad_grad2 - grad_grad1).norm() > 5.0E-12)
        deallog << "PQ" << k << ": second derivatives differ " << grad_grad1
                << " != " << grad_grad2 << std::endl;

      // Check if third derivative is ok
      Tensor<3, dim> third_derivative1 = p.template compute_derivative<3>(k, x);
      if ((third_derivative1 - third1[k]).norm() > 5.0E-12)
        deallog << 'P' << k << ": third derivatives differ "
                << third_derivative1 << " != " << third1[k] << std::endl;
      Tensor<3, dim> third_derivative2 = q.template compute_derivative<3>(k, x);
      if ((third_derivative2 - third2[k]).norm() > 5.0E-12)
        deallog << 'Q' << k << ": third derivatives differ "
                << third_derivative2 << " != " << third2[k] << std::endl;
      if ((third_derivative2 - third_derivative1).norm() > 5.0E-12)
        deallog << "PQ" << k << ": third derivatives differ "
                << third_derivative1 << " != " << third_derivative2
                << std::endl;

      // Check if third derivative is ok
      Tensor<4, dim> fourth_derivative1 =
        p.template compute_derivative<4>(k, x);
      if ((fourth_derivative1 - fourth1[k]).norm() > 5.0E-12)
        deallog << 'P' << k << ": fourth derivatives differ "
                << fourth_derivative1 << " != " << fourth1[k] << std::endl;
      Tensor<4, dim> fourth_derivative2 =
        q.template compute_derivative<4>(k, x);
      if ((fourth_derivative2 - fourth2[k]).norm() > 5.0E-12)
        deallog << 'Q' << k << ": fourth derivatives differ "
                << fourth_derivative2 << " != " << fourth2[k] << std::endl;
      if ((fourth_derivative2 - fourth_derivative1).norm() > 5.0E-12)
        deallog << "PQ" << k << ": fourth derivatives differ "
                << fourth_derivative1 << " != " << fourth_derivative2
                << std::endl;


      // finally output values, gradients, etc, to make sure that they are
      // not only consistent, but also correct.
      values1[k] *= std::pow(10., dim);
      gradients1[k] *= std::pow(10., dim);

      deallog << 'P' << k << "\t= " << values1[k] << "\tgradient\t";
      for (unsigned int d = 0; d < dim; ++d)
        deallog << gradients1[k][d] << '\t';
      deallog << "\t2nd\t";
      for (unsigned int d1 = 0; d1 < dim; ++d1)
        for (unsigned int d2 = 0; d2 < dim; ++d2)
          deallog << second1[k][d1][d2] << '\t';
      deallog << "\t3rd\t";
      for (unsigned int d1 = 0; d1 < dim; ++d1)
        for (unsigned int d2 = 0; d2 < dim; ++d2)
          for (unsigned int d3 = 0; d3 < dim; ++d3)
            deallog << third1[k][d1][d2][d3] << '\t';
      deallog << "\t4th\t";
      for (unsigned int d1 = 0; d1 < dim; ++d1)
        for (unsigned int d2 = 0; d2 < dim; ++d2)
          for (unsigned int d3 = 0; d3 < dim; ++d3)
            for (unsigned int d4 = 0; d4 < dim; ++d4)
              deallog << fourth1[k][d1][d2][d3][d4] << '\t';
      deallog << std::endl;
    }
  deallog << std::endl;
}


template <int dim>
void
check_tensor(const std::vector<Polynomial<double>> &v, const Point<dim> &x)
{
  TensorProductPolynomials<dim> p(v);

  std::vector<std::vector<Polynomial<double>>> pols(dim, v);
  AnisotropicPolynomials<dim>                  q(pols);

  check_poly(x, p, q);

  const std::vector<unsigned int> renumber =
    FETools::hierarchic_to_lexicographic_numbering<dim>(v.size() - 1);

  p.set_numbering(renumber);
  q.set_numbering(renumber);

  check_poly(x, p, q);
}


void
check_dimensions(const std::vector<Polynomial<double>> &p)
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

int
main()
{
  initlog();
  deallog.precision(10);

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
