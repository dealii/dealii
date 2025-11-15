// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


/*
 * Test FE_DGP and FE_DGQ with arbitrary polynomial spaces.
 */

#include <deal.II/base/quadrature_lib.h>

#include <deal.II/fe/fe_dgp.h>
#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_q1.h>

#include <deal.II/grid/grid_generator.h>

#include "../tests.h"


std::vector<Polynomials::Polynomial<double>>
create_orthonormal_basis(const unsigned int fe_degree)
{
  std::vector<Polynomials::Polynomial<double>> polynomials;

  std::vector<std::vector<double>> coefficients = {
    {{1.000000000000000e+00}},
    {{3.464101615137754e+00, -1.732050807568877e+00}},
    {{1.341640786499874e+01, -1.341640786499874e+01, 2.236067977499790e+00}},
    {{5.291502622129181e+01,
      -7.937253933193772e+01,
      3.174901573277509e+01,
      -2.645751311064591e+00}},
    {{2.100000000000000e+02,
      -4.200000000000000e+02,
      2.700000000000000e+02,
      -6.000000000000000e+01,
      3.000000000000000e+00}}};

  for (unsigned int c = 0; c <= fe_degree; ++c)
    {
      std::reverse(coefficients[c].begin(), coefficients[c].end());
      polynomials.emplace_back(coefficients[c]);
    }

  return polynomials;
}


template <int dim>
void
test(const FiniteElement<dim> &fe)
{
  MappingQ1<dim> mapping;
  QGauss<dim>    quadrature(fe.tensor_degree() + 1);

  FEValues<dim> fe_values(mapping,
                          fe,
                          quadrature,
                          update_JxW_values | update_values);

  Triangulation<dim> tria;
  GridGenerator::hyper_cube(tria);

  fe_values.reinit(tria.begin());

  FullMatrix<double> matrix(fe.n_dofs_per_cell(), fe.n_dofs_per_cell());

  for (const unsigned int q_index : fe_values.quadrature_point_indices())
    for (const unsigned int i : fe_values.dof_indices())
      for (const unsigned int j : fe_values.dof_indices())
        matrix[i][j] += fe_values.shape_value(i, q_index) *
                        fe_values.shape_value(j, q_index) *
                        fe_values.JxW(q_index);

  for (unsigned int i = 0; i < matrix.m(); ++i)
    for (unsigned int j = 0; j < matrix.n(); ++j)
      if (matrix[i][j] < 1e-8)
        matrix[i][j] = 0.0;

  matrix.print_formatted(deallog.get_file_stream());
}



int
main()
{
  initlog();

  const unsigned int dim       = 2;
  const unsigned int fe_degree = 3;

  const auto polynomials = create_orthonormal_basis(fe_degree);

  test(FE_DGP<dim>(polynomials));
  test(FE_DGQArbitraryPolynomialSpace<dim>(polynomials));
}
