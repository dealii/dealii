// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2015 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



#include <deal.II/distributed/tria.h>

#include <deal.II/grid/cell_id.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_tools.h>

#include <algorithm>

#include "../tests.h"

template <int dim>
void
test()
{
  const MPI_Comm mpi_communicator = MPI_COMM_WORLD;
  deallog << "dim = " << dim << std::endl;

  parallel::distributed::Triangulation<dim> tria(mpi_communicator);
  GridGenerator::hyper_cube(tria);
  tria.refine_global(2);

  using cell_iterator =
    typename parallel::distributed::Triangulation<dim>::active_cell_iterator;

  // Mark a small block at the corner of the hypercube
  std::vector<cell_iterator> ghost_cells_tria;
  cell_iterator              cell = tria.begin_active(), endc = tria.end();
  for (; cell != endc; ++cell)
    {
      if (cell->is_ghost() == true)
        ghost_cells_tria.push_back(cell);
    }
  std::sort(ghost_cells_tria.begin(), ghost_cells_tria.end());

  // Compute a halo layer around the locally owned cells.
  // These should all be ghost cells
  std::vector<cell_iterator> ghost_cell_halo_layer =
    GridTools::compute_ghost_cell_halo_layer(tria);
  AssertThrow(ghost_cell_halo_layer.size() > 0,
              ExcMessage("Ghost cell halo layer found."));
  AssertThrow(ghost_cell_halo_layer.size() == ghost_cells_tria.size(),
              ExcMessage("Ghost cell halo layer wrong size."));
  std::sort(ghost_cell_halo_layer.begin(), ghost_cell_halo_layer.end());

  for (unsigned int proc = 0;
       proc < Utilities::MPI::n_mpi_processes(mpi_communicator);
       ++proc)
    {
      if (proc == Utilities::MPI::this_mpi_process(mpi_communicator))
        {
          for (typename std::vector<cell_iterator>::const_iterator
                 it_1 = ghost_cells_tria.begin(),
                 it_2 = ghost_cell_halo_layer.begin();
               it_1 != ghost_cells_tria.end() &&
               it_2 != ghost_cell_halo_layer.end();
               ++it_1, ++it_2)
            {
              const cell_iterator &cell_1 = *it_1;
              const cell_iterator &cell_2 = *it_2;
              AssertThrow(cell_1->is_ghost() == true,
                          ExcMessage("Cell is not a ghost cell!"));
              AssertThrow(cell_2->is_ghost() == true,
                          ExcMessage("Halo cell is not a ghost cell!"));
              deallog << "Ghost " << cell_1->level() << ' ' << cell_1->index()
                      << ' ' << cell_1->id() << ' ' << cell_1->id().to_string()
                      << ' ' << "Halo " << cell_2->level() << ' '
                      << cell_2->index() << ' ' << cell_2->id() << ' '
                      << cell_2->id().to_string() << std::endl;
              AssertThrow(cell_2 == cell_1,
                          ExcMessage(
                            "Halo cell is not identical to ghost cell."));
            }
        }
    }
}


int
main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  MPILogInitAll                    log;

  test<2>();
  test<3>();

  return 0;
}
