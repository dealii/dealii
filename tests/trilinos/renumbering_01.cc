// ---------------------------------------------------------------------
//
// Copyright (C) 2009 - 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE.md at
// the top level directory of deal.II.
//
// ---------------------------------------------------------------------

// Test that dof renumbering is also reflected in TrilinosWrappers::SparseMatrix

#include <deal.II/base/conditional_ostream.h>
#include <deal.II/base/utilities.h>

#include <deal.II/distributed/shared_tria.h>
#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_renumbering.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe.h>
#include <deal.II/fe/fe_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>

#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/sparsity_tools.h>
#include <deal.II/lac/trilinos_sparse_matrix.h>
#include <deal.II/lac/trilinos_vector.h>

#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/vector_tools.h>

#include <mpi.h>

#include <iostream>
#include <vector>

#include "../tests.h"


class Test
{
private:
  MPI_Comm mpi_communicator;

  const unsigned int rank;
  const unsigned int n_ranks;

  parallel::distributed::Triangulation<2> triangulation;

  DoFHandler<2> dof_handler;

  FE_Q<2> fe;

  AffineConstraints<double> constraints;

  IndexSet locally_owned_dofs;
  IndexSet locally_relevant_dofs;

  TrilinosWrappers::SparseMatrix system_matrix;

public:
  Test(const bool do_renumber)
    : mpi_communicator(MPI_COMM_WORLD)
    , rank(Utilities::MPI::this_mpi_process(mpi_communicator))
    , n_ranks(Utilities::MPI::n_mpi_processes(mpi_communicator))
    , triangulation(mpi_communicator)
    , dof_handler(triangulation)
    , fe(1)
  {
    deallog << "Start";

    if (do_renumber)
      deallog << " with renumbering" << std::endl;
    else
      deallog << " without renumbering" << std::endl;

    GridGenerator::hyper_cube(triangulation);
    triangulation.refine_global(2);

    dof_handler.distribute_dofs(fe);

    constraints.clear();
    constraints.close();

    if (do_renumber)
      renumber();

    init_structures();

    deallog << "Finished";

    if (do_renumber)
      deallog << " with renumbering" << std::endl;
    else
      deallog << " without renumbering" << std::endl;
  }

  ~Test()
  {
    dof_handler.clear();
  }

private:
  void
  init_structures()
  {
    locally_owned_dofs = dof_handler.locally_owned_dofs();
    locally_owned_dofs.print(deallog);
    DoFTools::extract_locally_relevant_dofs(dof_handler, locally_relevant_dofs);
    locally_relevant_dofs.print(deallog);

    DynamicSparsityPattern sparsity_pattern(locally_relevant_dofs);
    DoFTools::make_sparsity_pattern(dof_handler,
                                    sparsity_pattern,
                                    constraints,
                                    /*keep constrained dofs*/ false);
    SparsityTools::distribute_sparsity_pattern(sparsity_pattern,
                                               locally_owned_dofs,
                                               MPI_COMM_WORLD,
                                               locally_relevant_dofs);

    system_matrix.reinit(locally_owned_dofs,
                         locally_owned_dofs,
                         sparsity_pattern,
                         MPI_COMM_WORLD);
    deallog << "local_range: " << system_matrix.local_range().first << " - "
            << system_matrix.local_range().second << std::endl;
  }

  void
  renumber()
  {
    // DoFRenumbering::Cuthill_McKee(dof_handler);

    locally_owned_dofs = dof_handler.locally_owned_dofs();

    std::vector<types::global_dof_index> new_number(dof_handler.n_dofs());
    for (types::global_dof_index i = 0; i < dof_handler.n_dofs(); i++)
      new_number[i] = dof_handler.n_dofs() - i - 1;

    std::vector<types::global_dof_index> local_new_number;
    for (IndexSet::ElementIterator dof = locally_owned_dofs.begin();
         dof != locally_owned_dofs.end();
         ++dof)
      local_new_number.push_back(new_number[*dof]);

    deallog << "n_dofs = " << dof_handler.n_dofs() << std::endl;
    deallog << "before renumbering:" << std::endl;
    locally_owned_dofs.print(dealii::deallog);
    dof_handler.renumber_dofs(local_new_number);

    deallog << "after renumbering:" << std::endl;
    locally_owned_dofs = dof_handler.locally_owned_dofs();
    locally_owned_dofs.print(dealii::deallog);
  }
};

int
main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi(argc, argv, 1);
  MPILogInitAll                    log;

  Test test1(false);
  MPI_Barrier(MPI_COMM_WORLD);
  Test test2(true);

  return 0;
}
