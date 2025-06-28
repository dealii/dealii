// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2013 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// check that cell->active_cell_index() works as advertised

#include <deal.II/base/geometry_info.h>
#include <deal.II/base/quadrature_lib.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria_accessor.h>

#include "../tests.h"



template <int dim>
void
check(const parallel::distributed::Triangulation<dim> &tria)
{
  unsigned int index = 0;
  for (typename Triangulation<dim>::active_cell_iterator cell =
         tria.begin_active();
       cell != tria.end();
       ++cell, ++index)
    AssertThrow(cell->active_cell_index() == index, ExcInternalError());

  AssertThrow(index == tria.n_active_cells(), ExcInternalError());
  AssertThrow(index >= tria.n_locally_owned_active_cells(), ExcInternalError());
}



template <int dim>
void
check()
{
  parallel::distributed::Triangulation<dim> tria(MPI_COMM_WORLD);
  GridGenerator::hyper_cube(tria);
  tria.refine_global(3);

  DoFHandler<dim> dof_handler(tria);

  check(tria);

  if (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0)
    deallog << "OK for " << dim << 'd' << std::endl;
}


int
main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  MPILogInitAll                    log;

  check<2>();
  check<3>();
}
