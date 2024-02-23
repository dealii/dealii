// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2004 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// simply check what happens when condensing matrices. This test was written
// when I changed a few things in the algorithm. By simply looping over all
// entries of the sparse matrix, we also check that things went right during
// compression of the sparsity pattern.

#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/sparsity_pattern.h>
#include <deal.II/lac/vector.h>

#include "../tests.h"


template <int dim>
void
test()
{
  deallog << dim << 'D' << std::endl;

  Triangulation<dim> triangulation;
  GridGenerator::hyper_cube(triangulation);

  // refine once, then refine first cell to
  // create hanging nodes
  triangulation.refine_global(1);
  triangulation.begin_active()->set_refine_flag();
  triangulation.execute_coarsening_and_refinement();
  deallog << "Number of cells: " << triangulation.n_active_cells() << std::endl;

  // set up a DoFHandler and compute hanging
  // node constraints for a Q2 element
  FE_Q<dim>       fe(2);
  DoFHandler<dim> dof_handler(triangulation);
  dof_handler.distribute_dofs(fe);
  deallog << "Number of dofs: " << dof_handler.n_dofs() << std::endl;

  AffineConstraints<double> constraints;
  DoFTools::make_hanging_node_constraints(dof_handler, constraints);
  constraints.close();
  deallog << "Number of constraints: " << constraints.n_constraints()
          << std::endl;

  // then set up a sparsity pattern and a
  // matrix on top of it
  SparsityPattern sparsity(dof_handler.n_dofs(),
                           dof_handler.n_dofs(),
                           dof_handler.max_couplings_between_dofs());
  DoFTools::make_sparsity_pattern(dof_handler, sparsity);
  constraints.condense(sparsity);
  SparseMatrix<double> A(sparsity);

  // then fill the matrix by setting up
  // bogus matrix entries
  std::vector<types::global_dof_index> local_dofs(fe.dofs_per_cell);
  FullMatrix<double> local_matrix(fe.dofs_per_cell, fe.dofs_per_cell);
  for (typename DoFHandler<dim>::active_cell_iterator cell =
         dof_handler.begin_active();
       cell != dof_handler.end();
       ++cell)
    {
      cell->get_dof_indices(local_dofs);
      local_matrix = 0;
      for (unsigned int i = 0; i < fe.dofs_per_cell; ++i)
        for (unsigned int j = 0; j < fe.dofs_per_cell; ++j)
          local_matrix(i, j) =
            (i + 1.) * (j + 1.) * (local_dofs[i] + 1.) * (local_dofs[j] + 1.);

      // copy local to global
      for (unsigned int i = 0; i < fe.dofs_per_cell; ++i)
        for (unsigned int j = 0; j < fe.dofs_per_cell; ++j)
          A.add(local_dofs[i], local_dofs[j], local_matrix(i, j));
    }

  // now condense away constraints from A
  constraints.condense(A);

  // and output what we have
  for (SparseMatrix<double>::const_iterator i = A.begin(); i != A.end(); ++i)
    deallog << i->row() << ' ' << i->column() << ' ' << i->value() << std::endl;
}



int
main()
{
  initlog();

  try
    {
      test<1>();
      test<2>();
      test<3>();
    }
  catch (const std::exception &exc)
    {
      deallog << std::endl
              << std::endl
              << "----------------------------------------------------"
              << std::endl;
      deallog << "Exception on processing: " << std::endl
              << exc.what() << std::endl
              << "Aborting!" << std::endl
              << "----------------------------------------------------"
              << std::endl;

      return 1;
    }
  catch (...)
    {
      deallog << std::endl
              << std::endl
              << "----------------------------------------------------"
              << std::endl;
      deallog << "Unknown exception!" << std::endl
              << "Aborting!" << std::endl
              << "----------------------------------------------------"
              << std::endl;
      return 1;
    };
}
