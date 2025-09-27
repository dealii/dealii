// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2010 - 2019 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// Tests that functions connected to 'pre_distributed_refinement'
// signals work on valid 'refine' and 'coarsen' flags


#include <deal.II/base/logstream.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_q.h>

#include <deal.II/grid/grid_generator.h>

#include <functional>

#include "../tests.h"


// check that cells flagged for coarsening cannot have siblings
// flagged for refinement
template <int dim>
void
verify_callback(const parallel::distributed::Triangulation<dim> &tria)
{
  for (const auto &cell : tria.active_cell_iterators())
    if (cell->coarsen_flag_set())
      {
        const auto parent = cell->parent();
        for (unsigned int c = 0; c < parent->n_children(); ++c)
          Assert(parent->child(c)->refine_flag_set() == false,
                 ExcInternalError());
      }
}



template <int dim>
void
test()
{
  parallel::distributed::Triangulation<dim> triangulation(
    MPI_COMM_WORLD, Triangulation<dim>::limit_level_difference_at_vertices);
  triangulation.signals.pre_distributed_refinement.connect(
    std::bind(verify_callback<dim>, std::cref(triangulation)));

  GridGenerator::hyper_cube(triangulation);
  triangulation.refine_global(2);

  DoFHandler<dim> dof_handler(triangulation);
  FE_Q<dim>       fe(1);
  dof_handler.distribute_dofs(fe);

  const unsigned int n_refinements[] = {0, 4, 3, 2};
  for (unsigned int i = 0; i < n_refinements[dim]; ++i)
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

  deallog << "OK" << std::endl;
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
