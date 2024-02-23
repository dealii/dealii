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



// this test, extracted from dof_constraints_09, used to fail with aborts

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
#include <deal.II/lac/block_sparse_matrix.h>
#include <deal.II/lac/block_sparsity_pattern.h>
#include <deal.II/lac/block_vector.h>

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
  triangulation.execute_coarsening_and_refinement();
  deallog << "Number of cells: " << triangulation.n_active_cells() << std::endl;

  // set up a DoFHandler and compute hanging
  // node constraints for a Q2 element
  FE_Q<dim>       fe(1);
  DoFHandler<dim> dof_handler(triangulation);
  dof_handler.distribute_dofs(fe);
  deallog << "Number of dofs: " << dof_handler.n_dofs() << std::endl;

  // then set up a sparsity pattern and a
  // matrix on top of it
  std::vector<unsigned int> block_sizes(2);
  block_sizes[0] = dof_handler.n_dofs() / 3;
  block_sizes[1] = dof_handler.n_dofs() - block_sizes[0];

  BlockSparsityPattern sparsity(2, 2);
  for (unsigned int i = 0; i < 2; ++i)
    for (unsigned int j = 0; j < 2; ++j)
      sparsity.block(i, j).reinit(block_sizes[i],
                                  block_sizes[j],
                                  dof_handler.max_couplings_between_dofs());
  sparsity.collect_sizes();

  DoFTools::make_sparsity_pattern(dof_handler, sparsity);
  sparsity.compress();
  BlockSparseMatrix<double> A(sparsity);

  // and output what we have. first for the
  // individual blocks, and later for all
  // together (which yielded an abort)
  deallog << "Blockwise output" << std::endl;
  for (unsigned int i = 0; i < 2; ++i)
    for (unsigned int j = 0; j < 2; ++j)
      for (SparseMatrix<double>::const_iterator k = A.block(i, j).begin();
           k != A.block(i, j).end();
           ++k)
        deallog << i << ' ' << j << ' ' << k->row() << ' ' << k->column() << ' '
                << k->value() << std::endl;

  deallog << "Global output" << std::endl;
  for (BlockSparseMatrix<double>::const_iterator i = A.begin(); i != A.end();
       ++i)
    deallog << i->block_row() << ' ' << i->block_column() << ' ' << i->row()
            << ' ' << i->column() << ' ' << i->value() << std::endl;
}



int
main()
{
  initlog();

  try
    {
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
