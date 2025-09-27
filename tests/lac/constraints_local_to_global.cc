// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2012 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// this function tests the correctness of the implementation of
// AffineConstraints<double>::distribute_local_to_global for FullMatrix by
// comparing the results with a sparse matrix. As a test case, we use a square
// mesh that is refined once globally and then the first cell is refined
// adaptively.

#include <deal.II/base/function.h>

#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_refinement.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>

#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/sparse_matrix.h>

#include <deal.II/numerics/vector_tools.h>

#include <complex>
#include <iostream>

#include "../tests.h"

template <int dim>
void
test()
{
  Triangulation<dim> tria;
  GridGenerator::hyper_cube(tria);
  tria.begin()->face(0)->set_boundary_id(1);
  tria.refine_global(1);
  tria.begin_active()->set_refine_flag();
  tria.execute_coarsening_and_refinement();

  FE_Q<dim>       fe(1);
  DoFHandler<dim> dof(tria);
  dof.distribute_dofs(fe);

  AffineConstraints<double> constraints;
  DoFTools::make_hanging_node_constraints(dof, constraints);
  VectorTools::interpolate_boundary_values(dof,
                                           1,
                                           Functions::ZeroFunction<dim>(),
                                           constraints);
  constraints.close();

  SparsityPattern sparsity;
  {
    DynamicSparsityPattern csp(dof.n_dofs(), dof.n_dofs());
    DoFTools::make_sparsity_pattern(dof, csp, constraints, false);
    sparsity.copy_from(csp);
  }
  SparseMatrix<double> sparse(sparsity);
  FullMatrix<double>   full(dof.n_dofs(), dof.n_dofs());

  FullMatrix<double> local_mat(fe.dofs_per_cell, fe.dofs_per_cell);
  std::vector<types::global_dof_index> local_dof_indices(fe.dofs_per_cell);

  // loop over cells, fill local matrix with
  // random values, insert both into sparse and
  // full matrix. Make some random entries equal
  // to zero
  typename DoFHandler<dim>::active_cell_iterator cell = dof.begin_active(),
                                                 endc = dof.end();
  unsigned int counter                                = 0;
  for (; cell != endc; ++cell)
    {
      for (unsigned int i = 0; i < fe.dofs_per_cell; ++i)
        for (unsigned int j = 0; j < fe.dofs_per_cell; ++j, ++counter)
          if (counter % 42 == 0)
            local_mat(i, j) = 0;
          else
            local_mat(i, j) = random_value<double>();
      cell->get_dof_indices(local_dof_indices);
      constraints.distribute_local_to_global(local_mat,
                                             local_dof_indices,
                                             sparse);
      constraints.distribute_local_to_global(local_mat,
                                             local_dof_indices,
                                             full);
    }

  // now check that the entries are indeed the
  // same by copying the sparse matrix into a
  // full matrix and checking the Frobenius norm
  // of the difference matrix
  FullMatrix<double> ref;
  ref.copy_from(sparse);
  full.add(-1., ref);
  deallog << "Difference between full and sparse matrix: "
          << full.frobenius_norm() << std::endl;
}


int
main()
{
  initlog();
  deallog << std::setprecision(2);
  deallog.get_file_stream() << std::setprecision(2);

  test<2>();
}
