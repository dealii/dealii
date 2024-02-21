// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2017 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// check n_nonzero_elements() for an empty matrix

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>

#include <deal.II/grid/grid_generator.h>

#include <deal.II/lac/sparsity_pattern.h>
#include <deal.II/lac/trilinos_sparse_matrix.h>
#include <deal.II/lac/trilinos_sparsity_pattern.h>

#include "../tests.h"


void
test()
{
  Triangulation<2> tria;
  GridGenerator::hyper_cube(tria);

  FE_Q<2>       fe(1);
  DoFHandler<2> dof_handler(tria);
  dof_handler.distribute_dofs(fe);

  Table<2, DoFTools::Coupling> coupling(1, 1);
  coupling.fill(DoFTools::none);

  // create an empty sparsity pattern
  TrilinosWrappers::SparsityPattern sparsity;
  sparsity.reinit(dof_handler.n_dofs(), dof_handler.n_dofs());
  DoFTools::make_sparsity_pattern(dof_handler,
                                  coupling,
                                  sparsity,
                                  AffineConstraints<double>(),
                                  false,
                                  Utilities::MPI::this_mpi_process(
                                    MPI_COMM_WORLD));
  sparsity.compress();

  // attach a sparse matrix to it
  TrilinosWrappers::SparseMatrix A;
  A.reinit(sparsity);

  // see how many nonzero elements it reports
  deallog << A.n_nonzero_elements() << std::endl;
}



int
main(int argc, char **argv)
{
  initlog();

  Utilities::MPI::MPI_InitFinalize mpi_initialization(
    argc, argv, testing_max_num_threads());

  try
    {
      test();
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
