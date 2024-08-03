// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2023 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// Tests point-wise evaluation of higher order derivatives with
// evaluate_tensor_product_higher_derivatives for a scalar function on FE_Q

#include <deal.II/base/function_lib.h>
#include <deal.II/base/point.h>
#include <deal.II/base/polynomial.h>
#include <deal.II/base/quadrature_lib.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_tools.h>

#include <deal.II/lac/vector.h>

#include <deal.II/matrix_free/tensor_product_point_kernels.h>

#include "../tests.h"

template <int dim, int derivative_order>
void
test(const unsigned int degree)
{
  FE_Q<dim> fe(degree);

  // go through all monomials of degree 'derivative_order + 1'
  std::vector<std::array<int, 3>> exponents;
  for (int i2 = 0; i2 <= (dim > 2 ? derivative_order : 0); ++i2)
    for (int i1 = 0; i1 <= (dim > 1 ? derivative_order : 0); ++i1)
      for (int i0 = 0; i0 <= derivative_order; ++i0)
        if (i0 + i1 + i2 == derivative_order)
          exponents.push_back(std::array<int, 3>{{i0, i1, i2}});

  deallog << "Evaluate derivative " << derivative_order << " in " << dim
          << "d with polynomial degree " << degree << std::endl;

  const std::vector<unsigned int> renumbering =
    FETools::lexicographic_to_hierarchic_numbering<dim>(degree);
  const std::vector<Polynomials::Polynomial<double>> polynomials =
    Polynomials::generate_complete_Lagrange_basis(
      QGaussLobatto<1>(degree + 1).get_points());

  Point<dim> p;
  for (unsigned int d = 0; d < dim; ++d)
    p[d] = 1.;

  std::vector<double> function_values(fe.dofs_per_cell);

  for (const std::array<int, 3> &exponent : exponents)
    {
      Tensor<1, dim> exp;
      for (unsigned int d = 0; d < dim; ++d)
        exp[d] = exponent[d];
      Functions::Monomial<dim> monomial(exp);
      for (unsigned int i = 0; i < fe.dofs_per_cell; ++i)
        function_values[i] = monomial.value(fe.get_unit_support_points()[i]);

      deallog << "Monomial [";
      for (unsigned int d = 0; d < dim; ++d)
        deallog << exponent[d] << (d == dim - 1 ? "" : " ");
      deallog << "]: ";
      const auto derivative =
        internal::evaluate_tensor_product_higher_derivatives<derivative_order>(
          polynomials, make_const_array_view(function_values), p, renumbering);

      for (unsigned int d = 0; d < derivative.dimension; ++d)
        deallog << (std::abs(derivative[d]) < 1e-11 ? 0. : derivative[d])
                << " ";
      deallog << std::endl;
    }
  deallog << std::endl;
}



int
main()
{
  initlog();
  deallog << std::setprecision(3);

  // test 2nd derivatives
  test<1, 2>(2);
  test<2, 2>(2);
  test<3, 2>(2);
  test<3, 2>(3);

  // test 3rd derivatives
  test<1, 3>(3);
  test<2, 3>(3);
  test<3, 3>(3);
  test<3, 3>(2);
  test<3, 3>(4);

  // test 4th derivatives
  test<1, 4>(4);
  test<2, 4>(4);
  test<3, 4>(4);
}
