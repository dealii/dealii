// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2018 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// Test that setting the kernel size in RelaxationBlock actually works. This
// test uses SVD-based smoothing so it can only be run with LAPACK.


#include <deal.II/base/function_lib.h>
#include <deal.II/base/quadrature_lib.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/relaxation_block.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/vector.h>

#include <deal.II/numerics/matrix_tools.h>

#include "../tests.h"


template <int dim>
void
make_stokes_matrix(const DoFHandler<dim> &dof_handler,
                   const Quadrature<dim> &quadrature_formula,
                   SparseMatrix<double>  &system_matrix)
{
  const FiniteElement<dim> &fe     = dof_handler.get_fe();
  const unsigned int        degree = fe.degree;
  system_matrix                    = 0;

  AffineConstraints<double> constraints;
  constraints.close();
  FEValues<dim>      fe_values(fe,
                          quadrature_formula,
                          update_values | update_quadrature_points |
                            update_JxW_values | update_gradients);
  const unsigned int dofs_per_cell = fe.dofs_per_cell;
  const unsigned int n_q_points    = quadrature_formula.size();
  FullMatrix<double> local_matrix(dofs_per_cell, dofs_per_cell);
  std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);
  const FEValuesExtractors::Vector     velocities(0);
  const FEValuesExtractors::Scalar     pressure(dim);
  std::vector<SymmetricTensor<2, dim>> symgrad_phi_u(dofs_per_cell);
  std::vector<double>                  div_phi_u(dofs_per_cell);
  std::vector<double>                  phi_p(dofs_per_cell);

  typename DoFHandler<dim>::active_cell_iterator cell =
                                                   dof_handler.begin_active(),
                                                 endc = dof_handler.end();
  for (; cell != endc; ++cell)
    {
      fe_values.reinit(cell);
      local_matrix = 0;
      for (unsigned int q = 0; q < n_q_points; ++q)
        {
          for (unsigned int k = 0; k < dofs_per_cell; ++k)
            {
              symgrad_phi_u[k] = fe_values[velocities].symmetric_gradient(k, q);
              div_phi_u[k]     = fe_values[velocities].divergence(k, q);
              phi_p[k]         = fe_values[pressure].value(k, q);
            }
          for (unsigned int i = 0; i < dofs_per_cell; ++i)
            for (unsigned int j = 0; j < dofs_per_cell; ++j)
              local_matrix(i, j) +=
                (2 * (symgrad_phi_u[i] * symgrad_phi_u[j]) -
                 div_phi_u[i] * phi_p[j] + phi_p[i] * div_phi_u[j]) *
                fe_values.JxW(q);
        }
      cell->get_dof_indices(local_dof_indices);
      constraints.distribute_local_to_global(local_matrix,
                                             local_dof_indices,
                                             system_matrix);
    }
}


template <int dim>
void
check()
{
  Triangulation<dim> tr(Triangulation<dim>::limit_level_difference_at_vertices);
  GridGenerator::hyper_cube(tr, -1, 1);
  tr.refine_global(1);

  FESystem<dim>   element(FESystem<dim>(FE_Q<dim>(2), dim), 1, FE_Q<dim>(1), 1);
  DoFHandler<dim> dof(tr);
  dof.distribute_dofs(element);
  dof.distribute_mg_dofs();

  QGauss<dim> quadrature(element.degree + 1);

  SparsityPattern sparsity(dof.n_dofs(), dof.n_dofs());
  DoFTools::make_sparsity_pattern(dof, sparsity);
  sparsity.compress();

  SparseMatrix<double> matrix;
  matrix.reinit(sparsity);

  make_stokes_matrix(dof, quadrature, matrix);

  BlockMask exclude_boundary_dofs(std::vector<bool>{true, false});

  using Smoother = RelaxationBlockJacobi<SparseMatrix<double>>;
  {
    Smoother::AdditionalData smoother_data;
    Smoother                 smoother;

    DoFTools::make_vertex_patches(smoother_data.block_list,
                                  dof,
                                  tr.n_levels() - 1,
                                  exclude_boundary_dofs);
    smoother_data.block_list.compress();
    smoother_data.inversion = PreconditionBlockBase<double>::svd;
    smoother_data.threshold = 1.e-8;

    smoother.initialize(matrix, smoother_data);
    smoother.log_statistics();
  }
  {
    Smoother::AdditionalData smoother_data;
    Smoother                 smoother;

    DoFTools::make_vertex_patches(smoother_data.block_list,
                                  dof,
                                  tr.n_levels() - 1,
                                  exclude_boundary_dofs);
    smoother_data.block_list.compress();
    smoother_data.inversion   = PreconditionBlockBase<double>::svd;
    smoother_data.kernel_size = 1;

    smoother.initialize(matrix, smoother_data);
    smoother.log_statistics();
  }
}



int
main()
{
  initlog();

  deallog.push("1d");
  check<1>();
  deallog.pop();
  deallog.push("2d");
  check<2>();
  deallog.pop();
  deallog.push("3d");
  check<3>();
  deallog.pop();
}
