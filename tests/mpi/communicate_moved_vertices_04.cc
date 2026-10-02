// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
// Copyright (C) 2026 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Detailed license information governing the source code and contributions
// can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
//
// -----------------------------------------------------------------------------

// Test Triangulation::communicate_locally_moved_vertices() on a 3x1 mesh
// with cells owned by ranks 2, 1, 0. The vertices shared by the middle and
// the right cell are owned by rank 0, but the middle cell is a ghost cell on
// rank 2 that does not neighbor any cell of rank 0.

#include <deal.II/distributed/fully_distributed_tria.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_description.h>

#include "../tests.h"


int
main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  MPILogInitAll                    log;

  constexpr int  dim  = 2;
  const MPI_Comm comm = MPI_COMM_WORLD;

  Triangulation<dim> base_triangulation;
  GridGenerator::subdivided_hyper_rectangle(base_triangulation,
                                            {3, 1},
                                            Point<dim>(0, 0),
                                            Point<dim>(3, 1));
  for (const auto &cell : base_triangulation.active_cell_iterators())
    cell->set_subdomain_id(2 - cell->index());

  parallel::fullydistributed::Triangulation<dim> triangulation(comm);
  triangulation.create_triangulation(
    TriangulationDescription::Utilities::create_description_from_triangulation(
      base_triangulation, comm));

  const std::vector<bool> locally_owned_vertices =
    GridTools::get_locally_owned_vertices(triangulation);

  std::vector<bool> vertex_moved(triangulation.n_vertices(), false);
  for (const auto &cell : triangulation.active_cell_iterators())
    if (cell->is_locally_owned())
      for (const auto v : cell->vertex_indices())
        {
          const auto index = cell->vertex_index(v);
          if (locally_owned_vertices[index] && !vertex_moved[index])
            {
              cell->vertex(v)[1] += 1.;
              vertex_moved[index] = true;
            }
        }

  triangulation.communicate_locally_moved_vertices(vertex_moved);

  for (const auto &cell : triangulation.active_cell_iterators())
    if (!cell->is_artificial())
      {
        deallog << "cell " << cell->id() << " (owner " << cell->subdomain_id()
                << "):";
        for (const auto v : cell->vertex_indices())
          deallog << " (" << cell->vertex(v) << ')';
        deallog << std::endl;
      }
}
