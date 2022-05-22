// ---------------------------------------------------------------------
//
// Copyright (C) 2020 - 2022 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE.md at
// the top level directory of deal.II.
//
// ---------------------------------------------------------------------



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

#include <deal.II/matrix_free/tensor_product_kernels.h>

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
               QIterated<dim>(QTrapez<1>(), 3).get_points();

  deallog << "Evaluate in " << dim << "d with polynomial degree " << degree
          << std::endl;
  for (const auto &p : evaluation_points)
    {
      const auto val = internal::evaluate_tensor_product_value_and_gradient(
        polynomials, coefficients, p, false, renumbering);
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
            polynomials, coefficients, p, true, renumbering);
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
