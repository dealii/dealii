// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// Check Trilinos sparse matrix iterators in a parallel
// context. One used to get errors when walking off the locally-owned
// set of rows of a sparsity pattern.

#include <deal.II/base/index_set.h>
#include <deal.II/base/utilities.h>

#include <deal.II/distributed/grid_refinement.h>
#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_renumbering.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>

#include <deal.II/grid/grid_generator.h>

#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/trilinos_sparse_matrix.h>
#include <deal.II/lac/trilinos_sparsity_pattern.h>
#include <deal.II/lac/trilinos_vector.h>

#include <iostream>

#include "../tests.h"


namespace dealiiLA
{
  using namespace dealii::TrilinosWrappers;
}

int
main(int argc, char *argv[])
{
  initlog();

  dealii::Utilities::MPI::MPI_InitFinalize mpiInitialization(argc, argv, 1);
  const MPI_Comm                          &mpiCommunicator(MPI_COMM_WORLD);

  dealii::parallel::distributed::Triangulation<3> tria(mpiCommunicator);
  dealii::GridGenerator::hyper_cube(tria);
  tria.refine_global(2);

  dealii::DoFHandler<3> dofHandler(tria);
  dealii::FE_Q<3>       fe(1);
  dofHandler.distribute_dofs(fe);

  dealii::DoFRenumbering::component_wise(dofHandler);

  dealii::IndexSet locallyOwnedDofs = dofHandler.locally_owned_dofs();

  dealiiLA::SparsityPattern sparsityPattern;
  sparsityPattern.reinit(locallyOwnedDofs, mpiCommunicator);
  dealii::DoFTools::make_sparsity_pattern(
    dofHandler,
    sparsityPattern,
    dealii::AffineConstraints(),
    false,
    dealii::Utilities::MPI::this_mpi_process(mpiCommunicator));
  sparsityPattern.compress();

  dealiiLA::SparseMatrix matrix(sparsityPattern);

  for (unsigned int i = 0; i < dofHandler.n_dofs(); ++i)
    if (sparsityPattern.row_is_stored_locally(i))
      {
        auto it = matrix.begin(i);

        // matrix.end(i) does not work for i=74,
        // because i=75 is not locally owned
        auto itend = matrix.end(i);
      }

  deallog << "OK" << std::endl;
}
