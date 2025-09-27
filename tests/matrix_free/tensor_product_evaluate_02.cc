// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2020 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// Tests point-wise evaluation of functions with
// evaluate_tensor_product_value_and_gradient for a vector-valued function
// represented by a linear function and comparing to the analytical value

#include <deal.II/base/point.h>
#include <deal.II/base/polynomial.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/symmetric_tensor.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_tools.h>

#include <deal.II/lac/vector.h>

#include <deal.II/matrix_free/tensor_product_point_kernels.h>

#include "../tests.h"

template <int dim>
void
test(const unsigned int degree)
{
  FE_Q<dim> fe(degree);

  Tensor<2, dim> matrix = unit_symmetric_tensor<dim>();
  for (unsigned int d = 0; d < dim; ++d)
    for (unsigned int e = 0; e < dim; ++e)
      matrix[d][e] += 0.1 * (d + 1 + 2 * e);
  Tensor<1, dim> offset;
  for (unsigned int d = 0; d < dim; ++d)
    offset[d] = 2 + 0.7 * d;

  std::vector<Tensor<1, dim>> coefficients(fe.dofs_per_cell);
  for (unsigned int i = 0; i < fe.dofs_per_cell; ++i)
    coefficients[i] = matrix * fe.get_unit_support_points()[i] + offset;

  const std::vector<Polynomials::Polynomial<double>> polynomials =
    Polynomials::generate_complete_Lagrange_basis(
      QGaussLobatto<1>(degree + 1).get_points());

  const std::vector<unsigned int> renumbering =
    FETools::lexicographic_to_hierarchic_numbering<dim>(degree);

  const std::vector<Point<dim>> evaluation_points =
    dim == 3 ? QGauss<dim>(2).get_points() :
               QIterated<dim>(QTrapezoid<1>(), 3).get_points();

  deallog << "Evaluate in " << dim << "d with polynomial degree " << degree
          << std::endl;
  for (const auto &p : evaluation_points)
    {
      const auto val = internal::evaluate_tensor_product_value_and_gradient(
        polynomials,
        make_const_array_view(coefficients),
        p,
        false,
        renumbering);
      deallog << "Value " << val.first << " vs " << matrix * p + offset
              << std::endl;
      deallog << "Gradient " << val.second << " vs " << transpose(matrix)
              << std::endl;
    }

  if (degree == 1)
    {
      deallog << "Evaluate d-linear shortcut in " << dim << 'd' << std::endl;
      for (const auto &p : evaluation_points)
        {
          const auto val = internal::evaluate_tensor_product_value_and_gradient(
            polynomials,
            make_const_array_view(coefficients),
            p,
            true,
            renumbering);
          deallog << "Value " << val.first << " vs " << matrix * p + offset
                  << std::endl;
          deallog << "Gradient " << val.second << " vs " << transpose(matrix)
                  << std::endl;
        }
    }
}



int
main()
{
  initlog();
  deallog << std::setprecision(9);

  for (unsigned int degree = 1; degree < 5; ++degree)
    test<1>(degree);
  for (unsigned int degree = 1; degree < 5; ++degree)
    test<2>(degree);
  for (unsigned int degree = 1; degree < 5; ++degree)
    test<3>(degree);
}
