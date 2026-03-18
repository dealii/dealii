// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
// Copyright (C) 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Detailed license information governing the source code and contributions
// can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
//
// -----------------------------------------------------------------------------



// like n_nonzeros_03, but this time create a BlockDynamicSparsityPattern
// and copy it into Tpetrawrappers::BlockSparsityPattern

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>

#include <deal.II/grid/grid_generator.h>

#include <deal.II/lac/block_sparsity_pattern.h>
#include <deal.II/lac/sparsity_pattern.h>
#include <deal.II/lac/trilinos_tpetra_block_sparse_matrix.h>

#include "../tests.h"


void
test(const Table<2, DoFTools::Coupling> &coupling)
{
  Triangulation<2> tria;
  GridGenerator::hyper_cube(tria);

  FESystem<2>   fe(FE_Q<2>(1), FE_Q<2>(2));
  DoFHandler<2> dof_handler(tria);
  dof_handler.distribute_dofs(fe);

  const IndexSet &locally_owned_total = dof_handler.locally_owned_dofs();
  const IndexSet  relevant_total =
    DoFTools::extract_locally_relevant_dofs(dof_handler);

  const std::vector<types::global_dof_index> dofs_per_block =
    DoFTools::count_dofs_per_fe_block(dof_handler, {{0, 1}});

  std::vector<IndexSet> locally_owned(2), relevant_set(2);
  locally_owned[0] = locally_owned_total.get_view(0, dofs_per_block[0]);
  locally_owned[1] =
    locally_owned_total.get_view(dofs_per_block[0], dof_handler.n_dofs());
  relevant_set[0] = relevant_total.get_view(0, dofs_per_block[0]);
  relevant_set[1] =
    relevant_total.get_view(dofs_per_block[0], dof_handler.n_dofs());

  // create an empty sparsity pattern
  dealii::BlockDynamicSparsityPattern dsp(relevant_set);
  DoFTools::make_sparsity_pattern(dof_handler, coupling, dsp);
  dsp.compress();

  dealii::LinearAlgebra::TpetraWrappers::BlockSparsityPattern sparsity;
  sparsity.copy_from(dsp);
  sparsity.compress();

  // attach a sparse matrix to it
  LinearAlgebra::TpetraWrappers::BlockSparseMatrix<double, MemorySpace::Host> A;
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
      // Test a few different couplings
      Table<2, DoFTools::Coupling> coupling(2, 2);
      coupling.fill(DoFTools::none);
      test(coupling);

      coupling(0, 0) = DoFTools::always;
      test(coupling);

      coupling(0, 0) = DoFTools::none;
      coupling(1, 1) = DoFTools::always;
      test(coupling);

      coupling.fill(DoFTools::always);
      test(coupling);
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
