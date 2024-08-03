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
// evaluate_tensor_product_hessian for a scalar function on FE_DGQ

#include <deal.II/base/point.h>
#include <deal.II/base/polynomial.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/symmetric_tensor.h>
#include <deal.II/base/vectorization.h>

#include <deal.II/fe/fe_dgq.h>

#include <deal.II/lac/vector.h>

#include <deal.II/matrix_free/tensor_product_point_kernels.h>

#include "../tests.h"

template <int dim>
void
test(const unsigned int degree)
{
  FE_DGQ<dim> fe(degree);

  // choose a symmetric matrix and then construct f(x) = x^T A x
  SymmetricTensor<2, dim> matrix = unit_symmetric_tensor<dim>();
  for (unsigned int d = 0; d < dim; ++d)
    for (unsigned int e = 0; e < dim; ++e)
      matrix[d][e] += 0.1 * (d + 1 + 2 * e);

  std::vector<double> coefficients(fe.dofs_per_cell);
  for (unsigned int i = 0; i < fe.dofs_per_cell; ++i)
    coefficients[i] = fe.get_unit_support_points()[i] *
                      (matrix * fe.get_unit_support_points()[i]);

  const std::vector<Polynomials::Polynomial<double>> polynomials =
    Polynomials::generate_complete_Lagrange_basis(
      QGaussLobatto<1>(degree + 1).get_points());

  const std::vector<Point<dim>> evaluation_points =
    dim == 3 ? QGauss<dim>(2).get_points() :
               QIterated<dim>(QTrapezoid<1>(), 3).get_points();

  deallog << "Evaluate in " << dim << "d with polynomial degree " << degree
          << std::endl;
  for (const auto &p : evaluation_points)
    {
      Point<dim, VectorizedArray<double>> p_vec;
      for (unsigned int v = 0; v < VectorizedArray<double>::size(); ++v)
        for (unsigned int d = 0; d < dim; ++d)
          p_vec[d][v] = p[d] + 0.01 * v;

      const auto hess = internal::evaluate_tensor_product_hessian(
        polynomials, make_const_array_view(coefficients), p_vec);

      std::cout << hess << "    " << matrix << std::endl;

      double error = 0;
      for (unsigned int v = 0; v < VectorizedArray<double>::size(); ++v)
        for (unsigned int d = 0; d < dim; ++d)
          for (unsigned int e = 0; e < dim; ++e)
            error += std::abs(hess[d][e][v] - 2. * matrix[e][d]);

      deallog << "Hessian error " << error << std::endl;
    }
}



int
main()
{
  initlog();
  deallog << std::setprecision(9);

  // We need at least degree 2 to correctly represent the Hessian of a
  // quadratic function
  for (unsigned int degree = 2; degree < 5; ++degree)
    test<1>(degree);
  for (unsigned int degree = 2; degree < 5; ++degree)
    test<2>(degree);
  for (unsigned int degree = 2; degree < 5; ++degree)
    test<3>(degree);
}
