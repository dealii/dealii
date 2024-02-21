// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2010 - 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// extracted from count_dofs_per_block that crashed at one point


#include <deal.II/base/tensor.h>
#include <deal.II/base/utilities.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/intergrid_map.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include <numeric>

#include "../tests.h"


template <int dim>
void
test()
{
  parallel::distributed::Triangulation<dim> triangulation(
    MPI_COMM_WORLD, Triangulation<dim>::limit_level_difference_at_vertices);

  GridGenerator::hyper_cube(triangulation);
  triangulation.refine_global(2);

  for (unsigned int i = 0; i < 2; ++i)
    {
      // refine one-fifth of cells randomly
      std::vector<bool> flags(triangulation.n_active_cells(), false);
      for (unsigned int k = 0; k < flags.size() / 5 + 1; ++k)
        flags[Testing::rand() % flags.size()] = true;
      // make sure there's at least one that
      // will be refined
      flags[0] = true;

      // refine triangulation
      unsigned int index = 0;
      for (typename Triangulation<dim>::active_cell_iterator cell =
             triangulation.begin_active();
           cell != triangulation.end();
           ++cell, ++index)
        if (flags[index])
          cell->set_refine_flag();
      AssertThrow(index == triangulation.n_active_cells(), ExcInternalError());

      // flag all other cells for coarsening
      // (this should ensure that at least
      // some of them will actually be
      // coarsened)
      index = 0;
      for (typename Triangulation<dim>::active_cell_iterator cell =
             triangulation.begin_active();
           cell != triangulation.end();
           ++cell, ++index)
        if (!flags[index])
          cell->set_coarsen_flag();

      triangulation.execute_coarsening_and_refinement();
    }
}


int
main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

  unsigned int myid = Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);
  if (myid == 0)
    {
      initlog();

      deallog.push("2d");
      test<2>();
      deallog.pop();

      deallog << "OK" << std::endl;
    }
  else
    {
      test<2>();
    }
}
