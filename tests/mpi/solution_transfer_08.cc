// ---------------------------------------------------------------------
//
// Copyright (C) 2021 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------



// Verify that SolutionTransfer works on parallel shared Triangulation
// objects with or without artificial cells.
//
// Test 'pure_refinement'.


#include <deal.II/distributed/shared_tria.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>

#include <deal.II/lac/vector.h>

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

  Vector<double> sol_old(dh.n_dofs());
  sol_old = 1.;

  SolutionTransfer<dim> soltrans(dh);
  soltrans.prepare_for_pure_refinement();

  tria.execute_coarsening_and_refinement();
  dh.distribute_dofs(fe);

  Vector<double> sol_new(dh.n_dofs());
  soltrans.refine_interpolate(sol_old, sol_new);
  for (unsigned int i = 0; i < sol_new.size(); ++i)
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
