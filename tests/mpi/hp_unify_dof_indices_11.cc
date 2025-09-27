// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2020 - 2021 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// have a 2x1 coarse mesh (or 2x1x1) and verify DoF indices in the hp-
// case with a FECollection that contains two finite elements that do
// not dominate each other. Here, a (FE_Q(1) x FE_Q(2)) and a
// (FE_Q(2) x FE_Q(1)) element on two separate subdomains face each
// other. the hp-code will unify DoF indices on boundaries between all
// subdomains.


#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>

#include <deal.II/grid/filtered_iterator.h>
#include <deal.II/grid/grid_generator.h>

#include <deal.II/hp/fe_collection.h>

#include "../tests.h"

#include "../test_grids.h"
#include "hp_unify_dof_indices.h"


template <int dim>
void
test()
{
  parallel::distributed::Triangulation<dim> triangulation(
    MPI_COMM_WORLD, Triangulation<dim>::limit_level_difference_at_vertices);
  TestGrids::hyper_line(triangulation, 2);
  Assert(triangulation.n_active_cells() == 2, ExcInternalError());

  hp::FECollection<dim> fe;
  fe.push_back(FESystem<dim>(FE_Q<dim>(1), 1, FE_Q<dim>(2), 1));
  fe.push_back(FESystem<dim>(FE_Q<dim>(2), 1, FE_Q<dim>(1), 1));

  DoFHandler<dim> dof_handler(triangulation);
  for (const auto &cell : dof_handler.active_cell_iterators() |
                            IteratorFilters::LocallyOwnedCell())
    {
      if (cell->id().to_string() == "0_0:")
        cell->set_active_fe_index(0);
      if (cell->id().to_string() == "1_0:")
        cell->set_active_fe_index(1);
    }
  dof_handler.distribute_dofs(fe);

  log_dof_diagnostics(dof_handler);
}


int
main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  MPILogInitAll                    log;

  deallog.push("2d");
  test<2>();
  deallog.pop();

  deallog.push("3d");
  test<3>();
  deallog.pop();
}
