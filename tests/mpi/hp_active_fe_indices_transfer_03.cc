// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2018 - 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// active FE indices transfer on repartitioning


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
  const unsigned int myid = Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);

  // ------ setup ------
  parallel::distributed::Triangulation<dim> tria(MPI_COMM_WORLD);

  GridGenerator::subdivided_hyper_cube(tria, 2);
  tria.refine_global(1);
  deallog << "cells before: " << tria.n_global_active_cells() << std::endl;

  // prepare FECollection with arbitrary number of entries
  hp::FECollection<dim> fe_collection;
  for (unsigned int i = 0; i < Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD);
       ++i)
    fe_collection.push_back(FE_Q<dim>(i + 1));

  DoFHandler<dim> dh(tria);

  for (auto &cell : dh.active_cell_iterators())
    if (cell->is_locally_owned())
      {
        // set active FE index
        if (!(cell->is_artificial()))
          cell->set_active_fe_index(myid);

        deallog << "cellid=" << cell->id()
                << " fe_index=" << cell->active_fe_index() << std::endl;
      }

  dh.distribute_dofs(fe_collection);

  // ----- transfer -----
  const parallel::CellWeights<dim> cell_weights(
    dh, parallel::CellWeights<dim>::ndofs_weighting({100000, 1}));

  tria.repartition();

  deallog << "cells after: " << tria.n_global_active_cells() << std::endl;

  // ------ verify ------
  // check if all children adopted the correct id
  for (auto &cell : dh.active_cell_iterators())
    if (!cell->is_artificial())
      {
        deallog << "cellid=" << cell->id()
                << " fe_index=" << cell->active_fe_index();
        if (cell->is_ghost())
          deallog << " ghost";
        deallog << std::endl;
      }

  // for further calculations, distribute dofs, i.e.
  // dh.distribute_dofs(fe_collection);

  // make sure no processor is hanging
  MPI_Barrier(MPI_COMM_WORLD);

  deallog << "OK" << std::endl;
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
