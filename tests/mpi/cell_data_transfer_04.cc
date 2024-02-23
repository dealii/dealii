// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// p::d::CellDataTransfer serialization with a different number of cpus


#include <deal.II/distributed/cell_data_transfer.templates.h>
#include <deal.II/distributed/tria.h>

#include <deal.II/grid/grid_generator.h>

#include <string>

#include "../tests.h"


template <int dim, int spacedim>
void
test()
{
  unsigned int myid    = Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);
  MPI_Comm     com_all = MPI_COMM_WORLD;
  MPI_Comm     com_small;

  // split the communicator in proc 0,1,2 and 3,4
  MPI_Comm_split(com_all, (myid < 3) ? 0 : 1, myid, &com_small);

  // write with small com
  if (myid < 3)
    {
      deallog << "writing with " << Utilities::MPI::n_mpi_processes(com_small)
              << std::endl;

      // ------ setup ------
      parallel::distributed::Triangulation<dim, spacedim> tria(com_small);
      GridGenerator::subdivided_hyper_cube(tria, 2);
      tria.refine_global(1);

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
                    << " parentid=" << cell_ids[cell->active_cell_index()]
                    << std::endl;
          }

      // ----- transfer -----
      parallel::distributed::
        CellDataTransfer<dim, spacedim, std::vector<unsigned int>>
          cell_data_transfer(tria);

      cell_data_transfer.prepare_for_coarsening_and_refinement(cell_ids);
      tria.save("file");
    }

  // make sure no processor is hanging
  MPI_Barrier(MPI_COMM_WORLD);

  {
    deallog << "reading with " << Utilities::MPI::n_mpi_processes(com_all)
            << std::endl;

    // ------ setup ------
    parallel::distributed::Triangulation<dim, spacedim> tria(com_all);
    GridGenerator::subdivided_hyper_cube(tria, 2);
    // triangulation has to be initialized with correct coarse cells

    // ----- transfer -----
    tria.load("file");

    parallel::distributed::
      CellDataTransfer<dim, spacedim, std::vector<unsigned int>>
        cell_data_transfer(tria);

    std::vector<unsigned int> cell_ids(tria.n_active_cells());
    cell_data_transfer.deserialize(cell_ids);

    // ------ verify ------
    // check if all children adopted the correct id
    for (auto &cell : tria.active_cell_iterators())
      if (cell->is_locally_owned())
        deallog << "cellid=" << cell->id()
                << " parentid=" << cell_ids[(cell->active_cell_index())]
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

  deallog.push("2d");
  test<2, 2>();
  deallog.pop();
  deallog.push("3d");
  test<3, 3>();
  deallog.pop();
}
