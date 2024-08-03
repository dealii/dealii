// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2015 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



#include <deal.II/distributed/fully_distributed_tria.h>

#include <deal.II/grid/cell_id.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_tools.h>

#include <map>
#include <vector>

#include "../tests.h"

// test compute_local_to_global_vertex_index_map() for p::fd::T

void
test()
{
  const MPI_Comm            mpi_comm = MPI_COMM_WORLD;
  const types::subdomain_id rank = Utilities::MPI::this_mpi_process(mpi_comm);

  Triangulation<3> serial_tria;
  GridGenerator::subdivided_hyper_cube_with_simplices(
    serial_tria, 4, 0.0, 1.0, true);
  parallel::fullydistributed::Triangulation<3> tria(mpi_comm);
  tria.copy_triangulation(serial_tria);

  const std::map<unsigned int, types::global_vertex_index> local_to_global_id =
    GridTools::compute_local_to_global_vertex_index_map(tria);

  deallog << "Number of vertices = " << tria.n_vertices() << std::endl;
  deallog << "Number of id pairs = " << local_to_global_id.size() << std::endl;

  for (const auto &cell : tria.active_cell_iterators())
    if (cell->is_locally_owned())
      for (const auto v_no : cell->vertex_indices())
        Assert(local_to_global_id.find(cell->vertex_index(v_no)) !=
                 local_to_global_id.end(),
               ExcMessage("Every local vertex should be in the map."));

  // Verify that we have a consistent and global ordering by checking physical
  // coordinates.
  std::map<types::global_vertex_index, Point<3>> vertices;
  for (const auto &pair : local_to_global_id)
    {
      AssertIndexRange(pair.first, tria.n_vertices());
      vertices[pair.second] = tria.get_vertices()[pair.first];
    }

  const auto all_vertices = Utilities::MPI::all_gather(mpi_comm, vertices);
  for (types::subdomain_id other_rank = 0; other_rank < all_vertices.size();
       ++other_rank)
    {
      unsigned int n_matches = 0;
      for (const auto &pair : all_vertices[other_rank])
        {
          const auto it = vertices.find(pair.first);
          if (it != vertices.end())
            {
              AssertThrow(it->second == pair.second,
                          ExcMessage("global vertices should match."));
              n_matches += (it->second == pair.second);
            }
        }
      deallog << "Found " << n_matches << " matching vertices between ranks "
              << rank << " and " << other_rank << std::endl;
    }
}


int
main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  MPILogInitAll                    log;

  test();

  deallog << "OK" << std::endl;

  return 0;
}
