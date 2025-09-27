// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2022 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// Test TensorProductMatrixSymmetricSum for zero (constrained) rows and columns.
// We consider a single cell with DBC applied to face 2*(dim-1).

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_q1.h>

#include <deal.II/grid/grid_generator.h>

#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/tensor_product_matrix.h>
#include <deal.II/lac/vector.h>

#include "../tests.h"

#include "../testmatrix.h"

template <int dim, typename Number>
std::tuple<FullMatrix<Number>, FullMatrix<Number>>
compute_reference_matrices(const unsigned int fe_degree)
{
  MappingQ1<dim> mapping;
  QGauss<dim>    quadrature(fe_degree + 1);
  FE_DGQ<dim>    fe(fe_degree);

  Triangulation<dim> tria;
  GridGenerator::hyper_cube(tria);

  DoFHandler<dim> dof_handler(tria);
  dof_handler.distribute_dofs(fe);

  const unsigned int n_dofs = fe.n_dofs_per_cell();

  FullMatrix<Number> mass_matrix_reference(n_dofs, n_dofs);
  FullMatrix<Number> derivative_matrix_reference(n_dofs, n_dofs);

  FEValues<dim> fe_values(mapping,
                          fe,
                          quadrature,
                          update_values | update_gradients | update_JxW_values);

  fe_values.reinit(dof_handler.begin());

  for (const unsigned int q_index : fe_values.quadrature_point_indices())
    for (const unsigned int i : fe_values.dof_indices())
      for (const unsigned int j : fe_values.dof_indices())
        {
          mass_matrix_reference(i, j) +=
            (fe_values.shape_value(i, q_index) *
             fe_values.shape_value(j, q_index) * fe_values.JxW(q_index));

          derivative_matrix_reference(i, j) +=
            (fe_values.shape_grad(i, q_index) *
             fe_values.shape_grad(j, q_index) * fe_values.JxW(q_index));
        }



  return std::tuple<FullMatrix<Number>, FullMatrix<Number>>{
    mass_matrix_reference, derivative_matrix_reference};
}

template <typename Number>
void
print(const FullMatrix<Number> &matrix, const std::string label)
{
  deallog << label << std::endl;
  matrix.print_formatted(deallog.get_file_stream(), 10, true, 15);
  deallog << std::endl << std::endl;
}


void
do_test(const bool zero_out_constraints)
{
  using Number                 = double;
  const unsigned int dim       = 2;
  const unsigned int n_rows_1d = 4;

  // compute 2D stiffness matrix
  const auto reference_matrices_2D =
    compute_reference_matrices<dim, Number>(n_rows_1d - 1);
  auto K_2D = std::get<1>(reference_matrices_2D);

  // ... and apply DBC on face 2*(dim-1)
  for (unsigned int j = 0; j < Utilities::pow(n_rows_1d, dim); ++j)
    for (unsigned int i = 0; i < Utilities::pow(n_rows_1d, dim - 1); ++i)
      {
        K_2D[i][j] = 0.0;
        K_2D[j][i] = 0.0;
      }

  if (zero_out_constraints == false)
    for (unsigned int i = 0; i < Utilities::pow(n_rows_1d, dim - 1); ++i)
      K_2D[i][i] = 1.0;

  // ... print
  print(K_2D, "K_2D:");

  if (zero_out_constraints)
    for (unsigned int i = 0; i < Utilities::pow(n_rows_1d, dim - 1); ++i)
      K_2D[i][i] = 1.0;

  K_2D.gauss_jordan();

  if (zero_out_constraints)
    for (unsigned int i = 0; i < Utilities::pow(n_rows_1d, dim - 1); ++i)
      K_2D[i][i] = 0.0;

  // ... print
  print(K_2D, "K_2D^-1:");

  // compute 1D stiffness and mass matrix
  const auto reference_matrices_1D =
    compute_reference_matrices<1, Number>(n_rows_1d - 1);

  // ... setup FDM
  std::array<Table<2, Number>, dim> mass_matrix;
  std::array<Table<2, Number>, dim> derivative_matrix;

  for (unsigned int d = 0; d < dim; ++d)
    {
      mass_matrix[d]       = std::get<0>(reference_matrices_1D);
      derivative_matrix[d] = std::get<1>(reference_matrices_1D);

      // apply constraints
      if ((d + 1) == dim)
        {
          for (unsigned int i = 0; i < n_rows_1d; ++i)
            {
              mass_matrix[d][i][0]       = 0.0;
              mass_matrix[d][0][i]       = 0.0;
              derivative_matrix[d][i][0] = 0.0;
              derivative_matrix[d][0][i] = 0.0;
            }

          if (zero_out_constraints == false)
            {
              mass_matrix[d][0][0]       = 1.0;
              derivative_matrix[d][0][0] = 1.0;
            }
        }
    }

  // ... print matrix
  TensorProductMatrixSymmetricSum<dim, Number, -1> fdm;
  fdm.reinit(mass_matrix, derivative_matrix);

  FullMatrix<Number> matrix(fdm.m(), fdm.m());

  for (unsigned int i = 0; i < fdm.m(); ++i)
    {
      Vector<Number> dst(fdm.m());
      Vector<Number> src(fdm.n());

      src[i] = 1.0;

      fdm.vmult(make_array_view(dst), make_array_view(src));

      for (unsigned int j = 0; j < fdm.m(); ++j)
        matrix[j][i] = dst[j];
    }

  print(matrix, "K_FDM:");

  // ... print inverse matrix
  for (unsigned int i = 0; i < fdm.m(); ++i)
    {
      Vector<Number> dst(fdm.m());
      Vector<Number> src(fdm.n());

      src[i] = 1.0;

      fdm.apply_inverse(make_array_view(dst), make_array_view(src));

      for (unsigned int j = 0; j < fdm.m(); ++j)
        matrix[j][i] = dst[j];
    }

  print(matrix, "K_FDM^-1:");
}


int
main()
{
  initlog();

  do_test(true);

  return 0;
}
