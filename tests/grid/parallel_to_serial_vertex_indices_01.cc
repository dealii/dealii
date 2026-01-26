// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2015 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// GridTools::parallel_to_serial_vertex_indices() test


#include <deal.II/distributed/fully_distributed_tria.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria.h>

#include <algorithm>
#include <cstdio>

#include "../tests.h"


template <int dim, int spacedim>
void
test()
{
  deallog << "Test " << dim << "d in " << spacedim << "d" << std::endl;

  Triangulation<dim, spacedim> serial_tria;
  GridGenerator::subdivided_hyper_cube(serial_tria, 5 - dim);

  parallel::fullydistributed::Triangulation<dim, spacedim> parallel_tria(
    MPI_COMM_WORLD);
  parallel_tria.copy_triangulation(serial_tria);

  const auto vertex_indices =
    GridTools::parallel_to_serial_vertex_indices(serial_tria, parallel_tria);

  const auto locally_owned_vertices =
    GridTools::get_locally_owned_vertices(parallel_tria);
  const auto n_locally_owned_vertices =
    std::count(locally_owned_vertices.begin(),
               locally_owned_vertices.end(),
               true);

  deallog << " serial vertices: " << serial_tria.n_vertices()
          << ", parallel vertices: " << parallel_tria.n_vertices()
          << " (owned vertices: " << n_locally_owned_vertices << ")"
          << std::endl;

  const auto &serial_vertices   = serial_tria.get_vertices();
  const auto &parallel_vertices = parallel_tria.get_vertices();

  unsigned int n_checked_vertices = 0;
  for (unsigned int i = 0; i < vertex_indices.size(); ++i)
    {
      auto serial_vertex_index = vertex_indices[i];
      if (serial_vertex_index != numbers::invalid_unsigned_int)
        {
          const auto &v_parallel = parallel_vertices[i];
          const auto &v_serial   = serial_vertices[serial_vertex_index];
          ++n_checked_vertices;
          if (v_parallel.distance(v_serial) > 1e-12)
            {
              deallog << "Error: vertex " << n_checked_vertices - 1
                      << " differs: parallel " << v_parallel << ", serial "
                      << v_serial << std::endl;
            }
        }
    }

  deallog << "Checked " << n_checked_vertices << " locally owned vertices."
          << std::endl;
}



int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_init(argc, argv);

  MPILogInitAll mpi_log_init_all;

  test<1, 1>();
  test<1, 2>();
  test<1, 3>();

  test<2, 2>();
  test<2, 3>();

  test<3, 3>();

  return 0;
}
