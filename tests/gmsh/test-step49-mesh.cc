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
// partitioned mesh and provide information about physical groups and
// boundary_ids

#include <deal.II/base/mpi.h>
#include <deal.II/base/utilities.h>

#include <deal.II/distributed/fully_distributed_tria.h>

#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/grid_out.h>

#include <fstream>

#include "../tests.h"

using namespace dealii;

// Main driver function.
// Initializes MPI, creates a fully distributed triangulation,
// reads a partitioned Gmsh mesh for each rank, and reports ownership stats.

int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc,
                                                      argv,
                                                      /*argc_mpi=*/1);
  MPILogInitAll                    mpilog;
  const MPI_Comm                   mpi_comm = MPI_COMM_WORLD;

  const unsigned int rank = Utilities::MPI::this_mpi_process(mpi_comm);
  const unsigned int size = Utilities::MPI::n_mpi_processes(mpi_comm);


  // Create fully distributed triangulation
  parallel::fullydistributed::Triangulation<2, 2> tria(mpi_comm);

  // Use GridIn to read the mesh
  GridIn<2, 2> grid_in;
  grid_in.attach_triangulation(tria);
  grid_in.read_partitioned_msh(SOURCE_DIR "/../grid/grids/mesh-from-step49",
                               "msh");


  // Output ownership information
  deallog << "Rank " << rank << " owns " << tria.n_active_cells() << " cells "
          << "and " << tria.n_vertices() << " vertices." << std::endl;


  // Count locally owned and ghost cells
  unsigned int n_locally_owned_cells = 0;
  unsigned int n_ghost_cells         = 0;
  for (const auto &cell : tria.active_cell_iterators())
    {
      if (cell->is_locally_owned())
        ++n_locally_owned_cells;
      else if (cell->is_ghost())
        ++n_ghost_cells;
    }

  // Collect boundary IDs from faces
  std::map<types::boundary_id, unsigned int> boundary_id_counts;
  for (const auto &cell : tria.active_cell_iterators())
    if (cell->is_locally_owned())
      for (unsigned int f = 0; f < GeometryInfo<2>::faces_per_cell; ++f)
        if (cell->face(f)->at_boundary())
          ++boundary_id_counts[cell->face(f)->boundary_id()];

  // Collect material (physical surface) IDs from cells
  std::map<types::material_id, unsigned int> material_id_counts;
  for (const auto &cell : tria.active_cell_iterators())
    if (cell->is_locally_owned())
      ++material_id_counts[cell->material_id()];

  // Report detailed statistics
  deallog << "Rank " << rank << " owns " << n_locally_owned_cells
          << " cells and has " << n_ghost_cells << " ghost cells." << std::endl;

  deallog << "Rank " << rank << " boundary ID counts:" << std::endl;
  for (const auto &p : boundary_id_counts)
    deallog << "  ID " << p.first << ": " << p.second << " faces" << std::endl;

  deallog << "Rank " << rank << " material ID counts:" << std::endl;
  for (const auto &p : material_id_counts)
    deallog << "  ID " << p.first << ": " << p.second << " cells" << std::endl;

  // Final output
  std::cout << "Rank " << rank << " finished deallog output" << std::endl;
  std::cout << "Rank " << rank << " exiting main()" << std::endl;

  return 0;
}
