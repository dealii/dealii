// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2010 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// SolutionTransfer locked up when a process has no locally owned
// cells. this was fixed with r24007

#include <deal.II/base/function.h>
#include <deal.II/base/tensor.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>

#include <deal.II/grid/filtered_iterator.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/lac/petsc_vector.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/solution_transfer.h>
#include <deal.II/numerics/vector_tools.h>

#include "../tests.h"


template <int dim>
void
test()
{
  parallel::distributed::Triangulation<2> tria(
    MPI_COMM_WORLD,
    typename Triangulation<2>::MeshSmoothing(
      Triangulation<2>::smoothing_on_refinement |
      Triangulation<2>::smoothing_on_coarsening));

  GridGenerator::hyper_cube(tria, -1.0, 1.0);

  FE_Q<2>       fe(1);
  DoFHandler<2> dh(tria);

  dh.distribute_dofs(fe);

  IndexSet locally_owned_dofs    = dh.locally_owned_dofs();
  IndexSet locally_relevant_dofs = DoFTools::extract_locally_relevant_dofs(dh);

  PETScWrappers::MPI::Vector solution(locally_owned_dofs,
                                      locally_relevant_dofs,
                                      MPI_COMM_WORLD);

  SolutionTransfer<2, PETScWrappers::MPI::Vector> soltrans(dh);

  tria.set_all_refine_flags();
  tria.prepare_coarsening_and_refinement();

  soltrans.prepare_for_coarsening_and_refinement(solution);


  tria.execute_coarsening_and_refinement();

  dh.distribute_dofs(fe);
  locally_owned_dofs    = dh.locally_owned_dofs();
  locally_relevant_dofs = DoFTools::extract_locally_relevant_dofs(dh);

  PETScWrappers::MPI::Vector tmp(locally_owned_dofs, MPI_COMM_WORLD);

  soltrans.interpolate(tmp);

  solution.reinit(locally_owned_dofs, locally_relevant_dofs, MPI_COMM_WORLD);
  solution = tmp;

  // make sure no processor is
  // hanging
  MPI_Barrier(MPI_COMM_WORLD);

  if (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0)
    deallog << "OK" << std::endl;
}


int
main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);


  unsigned int myid = Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);


  deallog.push(Utilities::int_to_string(myid));

  if (myid == 0)
    {
      initlog();

      deallog.push("2d");
      test<2>();
      deallog.pop();
      deallog.push("3d");
      test<3>();
      deallog.pop();
    }
  else
    {
      test<2>();
      test<3>();
    }
}
