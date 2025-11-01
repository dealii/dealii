// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2013 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// check AffineConstraints<double>::condense(Vector) for parallel vectors. this
// function used to crash
//
// original test case by Daniel Arndt

#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>

#include <deal.II/grid/grid_generator.h>

#include <deal.II/lac/petsc_vector.h>

#include "../tests.h"


void
test()
{
  MPI_Comm                                mpi_communicator(MPI_COMM_WORLD);
  parallel::distributed::Triangulation<3> triangulation(mpi_communicator);

  FE_Q<3> fe(1);

  DoFHandler<3> dof_handler(triangulation);
  GridGenerator::hyper_cube(triangulation);
  triangulation.refine_global(1);

  dof_handler.distribute_dofs(fe);

  const IndexSet &locally_owned_dofs = dof_handler.locally_owned_dofs();
  const IndexSet  locally_relevant_dofs =
    DoFTools::extract_locally_relevant_dofs(dof_handler);

  PETScWrappers::MPI::Vector force;
  force.reinit(locally_owned_dofs, mpi_communicator);
  Assert(!force.has_ghost_elements(), ExcInternalError());

  AffineConstraints<PetscScalar> constraints(locally_owned_dofs,
                                             locally_relevant_dofs);
  constraints.clear();
  {
    const IndexSet boundary_dofs =
      DoFTools::extract_boundary_dofs(dof_handler, ComponentMask(1, true));

    unsigned int first_nboundary_dof = 0;
    while (boundary_dofs.is_element(first_nboundary_dof))
      first_nboundary_dof++;

    if (locally_relevant_dofs.is_element(first_nboundary_dof))
      {
        constraints.add_line(first_nboundary_dof);
        for (unsigned int i = 0; i < dof_handler.n_dofs(); ++i)
          if (boundary_dofs.is_element(i) == true)
            constraints.add_entry(first_nboundary_dof, i, -1);
      }
  }
  constraints.close();

  constraints.condense(force);
}


int
main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

  unsigned int myid = Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);

  if (myid == 0)
    {
      initlog();
      deallog << std::setprecision(4);

      test();

      deallog << "OK" << std::endl;
    }
  else
    test();


  return 0;
}
