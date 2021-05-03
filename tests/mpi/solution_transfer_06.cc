// ---------------------------------------------------------------------
//
// Copyright (C) 2020 by the deal.II authors
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



// Test parallel::distributed::SolutionTransfer for FE_Nothing
// (see also https://github.com/dealii/dealii/issues/10570).

#include <deal.II/base/function.h>

#include <deal.II/distributed/solution_transfer.h>
#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_nothing.h>
#include <deal.II/fe/fe_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_refinement.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/lac/la_parallel_vector.h>

#include <iostream>
#include <vector>

#include "../tests.h"


template <int dim>
void
transfer(const MPI_Comm &comm)
{
  AssertDimension(Utilities::MPI::n_mpi_processes(comm), 1);

  parallel::distributed::Triangulation<dim> tria(comm);
  GridGenerator::hyper_cube(tria);
  tria.refine_global(1);

  hp::FECollection<dim> fe;
  fe.push_back(FE_Q<dim>(1));
  fe.push_back(FE_Nothing<dim>());

  // create a DoFHandler on which we have both cells with FE_Q as well as
  // FE_Nothing
  DoFHandler<dim> dof_handler(tria);
  dof_handler.begin(0)->child(0)->set_active_fe_index(1);

  IndexSet locally_relevant_dofs;
  DoFTools::extract_locally_relevant_dofs(dof_handler, locally_relevant_dofs);

  LinearAlgebra::distributed::Vector<double> solution(
    dof_handler.locally_owned_dofs(), locally_relevant_dofs, comm);
  AffineConstraints<double> cm;
  cm.close();

  dof_handler.distribute_dofs(fe);
  solution.reinit(dof_handler.n_dofs());

  for (unsigned int i = 0; i < solution.size(); ++i)
    solution(i) = i;

  parallel::distributed::
    SolutionTransfer<dim, LinearAlgebra::distributed::Vector<double>>
      soltrans(dof_handler);

  for (const auto &cell : tria.active_cell_iterators())
    cell->set_refine_flag();

  LinearAlgebra::distributed::Vector<double> old_solution = solution;

  tria.prepare_coarsening_and_refinement();
  soltrans.prepare_for_coarsening_and_refinement(old_solution);
  tria.execute_coarsening_and_refinement();
  dof_handler.distribute_dofs(fe);
  solution.reinit(dof_handler.n_dofs());

  soltrans.interpolate(solution);

  deallog << "OK" << std::endl;
}


int
main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

  const MPI_Comm comm = MPI_COMM_WORLD;
  initlog();

  deallog.push("2D solution transfer");
  transfer<2>(comm);
  deallog.pop();

  deallog.push("3D solution transfer");
  transfer<3>(comm);
  deallog.pop();
}
