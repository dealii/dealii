// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2018 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// Cell Weights Test
// -----------------
// Create a larger 8x8(x8) grid, on which all cells are associated with
// a Q1 element, besides the very first one which has a Q7 element.
// We choose a cell weighting algorithm based on the number of degrees
// of freedom and check if load is balanced as expected after
// repartitioning the triangulation. The expected accumulated weight on
// each processor should correlate to the sum of all degrees of
// freedom on all cells of the corresponding subdomain.
//
// This test works on a parallel::distributed::Triangulation.
//
// This test runs on a larger domain with a Lagrangian element of higher
// order, compared to the previous one. If we would have used a Q7
// element on the smaller grid, load balancing would fail in such a way
// that the second (last) processor owns the whole domain -- p4est wants
// to 'cut' its tree on a parent branch that does not exist in this case.


#include <deal.II/distributed/cell_weights.h>
#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_q.h>

#include <deal.II/grid/grid_generator.h>

#include "../tests.h"


template <int dim>
void
test()
{
  parallel::distributed::Triangulation<dim> tria(MPI_COMM_WORLD);

  GridGenerator::hyper_cube(tria);
  tria.refine_global(3);

  // Apply ndof cell weights.
  hp::FECollection<dim> fe_collection;
  fe_collection.push_back(FE_Q<dim>(1));
  fe_collection.push_back(FE_Q<dim>(7));

  DoFHandler<dim> dh(tria);

  // default: active_fe_index = 0
  for (auto &cell : dh.active_cell_iterators())
    if (cell->is_locally_owned())
      if (cell->id().to_string() == "0_3:000")
        cell->set_active_fe_index(1);

  dh.distribute_dofs(fe_collection);

  deallog << "Number of cells before repartitioning: "
          << tria.n_locally_owned_active_cells() << std::endl;
  {
    unsigned int dof_counter = 0;
    for (auto &cell : dh.active_cell_iterators())
      if (cell->is_locally_owned())
        dof_counter += cell->get_fe().dofs_per_cell;
    deallog << "  Cumulative dofs per cell: " << dof_counter << std::endl;
  }


  const parallel::CellWeights<dim> cell_weights(
    dh, parallel::CellWeights<dim>::ndofs_weighting({1, 1}));

  tria.repartition();


  deallog << "Number of cells after repartitioning: "
          << tria.n_locally_owned_active_cells() << std::endl;
  {
    unsigned int dof_counter = 0;
    for (auto &cell : dh.active_cell_iterators())
      if (cell->is_locally_owned())
        dof_counter += cell->get_fe().dofs_per_cell;
    deallog << "  Cumulative dofs per cell: " << dof_counter << std::endl;
  }

  if constexpr (running_in_debug_mode())
    {
      parallel::distributed::Triangulation<dim> other_tria(MPI_COMM_WORLD);
      GridGenerator::hyper_cube(other_tria);
      other_tria.refine_global(3);

      dh.reinit(other_tria);
      dh.distribute_dofs(fe_collection);

      try
        {
          tria.repartition();
        }
      catch (const ExceptionBase &e)
        {
          deallog << e.get_exc_name() << std::endl;
        }
    }
  else
    {
      deallog
        << "ExcMessage(\"Triangulation associated with the DoFHandler has changed!\")"
        << std::endl;
    }

  // make sure no processor is hanging
  MPI_Barrier(MPI_COMM_WORLD);

  deallog << "OK" << std::endl;
}


int
main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  MPILogInitAll                    log;

  deal_II_exceptions::disable_abort_on_exception();

  deallog.push("2d");
  test<2>();
  deallog.pop();
  deallog.push("3d");
  test<3>();
  deallog.pop();
}
