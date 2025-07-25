// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2021 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found incd
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

// Test for Gmsh::read_partitioned_msh() on a fully distributed
// triangulation such that each MPI process only reads its own
// partitioned mesh.

#include <deal.II/base/mpi.h>
#include <deal.II/base/utilities.h>

#include <deal.II/distributed/fully_distributed_tria.h>

#include <deal.II/gmsh/utilities.h>

#include <deal.II/grid/grid_out.h>

#include <fstream>

#include "../tests.h"

using namespace dealii;

/**
 * @brief Main driver function
 *
 * Initializes MPI, creates a fully distributed triangulation, reads a
 * partitioned Gmsh mesh for each rank, and reports ownership stats.
 */
int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc,
                                                      argv,
                                                      /*argc_mpi=*/1);
  MPILogInitAll                    mpilog;
  const MPI_Comm                   mpi_comm = MPI_COMM_WORLD;

  int rank = 0;
  int size = 0;
  MPI_Comm_rank(mpi_comm, &rank);
  MPI_Comm_size(mpi_comm, &size);

  // Create fully distributed triangulation
  parallel::fullydistributed::Triangulation<2, 2> tria(mpi_comm);

  // Read the partitioned Gmsh mesh for this rank
  Gmsh::read_partitioned_msh(tria,
                             mpi_comm,
                             SOURCE_DIR "/../grid/grids/unit-square");

  // Output ownership information
  deallog << "Rank " << rank << " owns " << tria.n_active_cells() << " cells "
          << "and " << tria.n_vertices() << " vertices." << std::endl;

  unsigned int n_locally_owned_cells = 0;
  unsigned int n_ghost_cells         = 0;

  // Count locally owned and ghost cells
  for (const auto &cell : tria.active_cell_iterators())
    {
      if (cell->is_locally_owned())
        ++n_locally_owned_cells;
      else if (cell->is_ghost())
        ++n_ghost_cells;
    }

  // Report detailed statistics
  deallog << "Rank " << rank << " owns " << n_locally_owned_cells
          << " cells and has " << n_ghost_cells << " ghost cells." << std::endl;

  // Final output
  std::cout << "Rank " << rank << " finished deallog output" << std::endl;
  std::cout << "Rank " << rank << " exiting main()" << std::endl;

  return 0;
}
