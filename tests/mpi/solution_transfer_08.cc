// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
// Copyright (C) 2021 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Detailed license information governing the source code and contributions
// can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
//
// -----------------------------------------------------------------------------



// Verify that SolutionTransfer works on parallel shared Triangulation
// objects with or without artificial cells.
//
// Test 'pure_refinement'.


#include <deal.II/distributed/shared_tria.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>

#include <deal.II/lac/la_parallel_vector.h>

#include <deal.II/numerics/solution_transfer.h>

#include "../tests.h"


template <int dim>
void
test(const bool allow_artificial_cells)
{
  parallel::shared::Triangulation<dim> tria(MPI_COMM_WORLD,
                                            Triangulation<dim>::none,
                                            allow_artificial_cells);
  GridGenerator::subdivided_hyper_cube(tria, 2);

  FE_Q<dim>       fe(1);
  DoFHandler<dim> dh(tria);
  dh.distribute_dofs(fe);

  tria.begin_active()->set_refine_flag();
  tria.prepare_coarsening_and_refinement();

  using VectorType = LinearAlgebra::distributed::Vector<double>;
  VectorType sol_old(dh.locally_owned_dofs(),
                     DoFTools::extract_locally_relevant_dofs(dh),
                     MPI_COMM_WORLD);
  sol_old = 1.;
  sol_old.update_ghost_values();

  SolutionTransfer<dim, VectorType> soltrans(dh);
  soltrans.prepare_for_coarsening_and_refinement(sol_old);

  tria.execute_coarsening_and_refinement();
  dh.distribute_dofs(fe);

  VectorType sol_new(dh.locally_owned_dofs(),
                     DoFTools::extract_locally_relevant_dofs(dh),
                     MPI_COMM_WORLD);
  soltrans.interpolate(sol_new);
  for (const auto i : sol_new.locally_owned_elements())
    AssertThrow(sol_new[i] == 1., ExcInternalError());

  deallog << "OK" << std::endl;
}


int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  MPILogInitAll                    log;

  deallog.push("1d");
  test<1>(false);
  test<1>(true);
  deallog.pop();
  deallog.push("2d");
  test<2>(false);
  test<2>(true);
  deallog.pop();
  deallog.push("3d");
  test<3>(false);
  test<3>(true);
  deallog.pop();
}
