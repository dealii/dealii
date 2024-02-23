// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2019 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// p::d::CellDataTransfer test by refinement


#include <deal.II/distributed/cell_data_transfer.templates.h>
#include <deal.II/distributed/tria.h>

#include <deal.II/grid/grid_generator.h>

#include <string>

#include "../tests.h"


template <int dim, int spacedim>
void
test()
{
  // ------ setup ------
  parallel::distributed::Triangulation<dim, spacedim> tria(MPI_COMM_WORLD);

  GridGenerator::subdivided_hyper_cube(tria, 2);
  tria.refine_global(1);
  deallog << "cells before: " << tria.n_global_active_cells() << std::endl;

  // ----- prepare -----
  // set refinement/coarsening flags
  for (auto &cell : tria.active_cell_iterators())
    if (cell->is_locally_owned())
      {
        if (cell->id().to_string() == "0_1:0")
          cell->set_refine_flag();
        else if (cell->parent()->id().to_string() ==
                 ((dim == 2) ? "3_0:" : "7_0:"))
          cell->set_coarsen_flag();
      }

  // ----- gather -----
  // store parent id of all cells
  std::vector<unsigned int> cell_ids(tria.n_active_cells());
  for (auto &cell : tria.active_cell_iterators())
    if (cell->is_locally_owned())
      {
        const std::string  parent_cellid = cell->parent()->id().to_string();
        const unsigned int parent_coarse_cell_id =
          (unsigned int)std::stoul(parent_cellid);
        cell_ids[cell->active_cell_index()] = parent_coarse_cell_id;

        deallog << "cellid=" << cell->id()
                << " parentid=" << cell_ids[cell->active_cell_index()];
        if (cell->coarsen_flag_set())
          deallog << " coarsening";
        else if (cell->refine_flag_set())
          deallog << " refining";
        deallog << std::endl;
      }

  // ----- transfer -----
  parallel::distributed::
    CellDataTransfer<dim, spacedim, std::vector<unsigned int>>
      cell_data_transfer(tria);

  cell_data_transfer.prepare_for_coarsening_and_refinement(cell_ids);
  tria.execute_coarsening_and_refinement();
  deallog << "cells after: " << tria.n_global_active_cells() << std::endl;

  cell_ids.resize(tria.n_active_cells());
  cell_data_transfer.unpack(cell_ids);

  // ------ verify ------
  // check if all children adopted the correct id
  for (auto &cell : tria.active_cell_iterators())
    if (cell->is_locally_owned())
      deallog << "cellid=" << cell->id()
              << " parentid=" << cell_ids[(cell->active_cell_index())]
              << std::endl;

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
  test<2, 2>();
  deallog.pop();
  deallog.push("3d");
  test<3, 3>();
  deallog.pop();
}
