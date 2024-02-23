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

// Test insert_global_particles. Make sure we don't lose particles
// along the way. Test the case where all particles are owned by a
// single mpi process, and where we add one property for particle.

#include <deal.II/base/mpi.h>
#include <deal.II/base/utilities.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/fe/mapping_q.h>

#include <deal.II/grid/filtered_iterator.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/grid_tools_cache.h>

#include <deal.II/particles/particle_handler.h>

#include <unistd.h>

#include <iostream>

#include "../tests.h"

template <int dim, int spacedim>
void
test()
{
  parallel::distributed::Triangulation<dim, spacedim> tr(MPI_COMM_WORLD);

  GridGenerator::hyper_cube(tr);
  tr.refine_global(2);

  MappingQ<dim, spacedim> mapping(1);

  Particles::ParticleHandler<dim, spacedim> particle_handler(tr, mapping, 1);

  const unsigned int n_points = 10;
  const unsigned int my_cpu = Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);
  const unsigned int n_cpus = Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD);
  Testing::srand(my_cpu + 1);

  // Distribute the local points to the processor that owns them
  // on the triangulation
  auto my_bounding_box = GridTools::compute_mesh_predicate_bounding_box(
    tr, IteratorFilters::LocallyOwnedCell());

  auto global_bounding_boxes =
    Utilities::MPI::all_gather(MPI_COMM_WORLD, my_bounding_box);

  std::vector<Point<spacedim>>     points(n_points);
  std::vector<std::vector<double>> properties(n_points, {1. * my_cpu});

  for (auto &p : points)
    p = random_point<spacedim>(0, 0.1);

  auto cpu_to_index =
    particle_handler.insert_global_particles(points,
                                             global_bounding_boxes,
                                             properties);

  for (auto p : particle_handler)
    {
      deallog << "Particle : " << p.get_id() << ", property: "
              << static_cast<unsigned int>(p.get_properties()[0]) << std::endl;
    }
}



int
main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

  MPILogInitAll all;

  deallog.push("2d/2d");
  test<2, 2>();
  deallog.pop();

  // [TODO]: There is a bug that prevents this from working in dimensions <2,3>
  // The construction of the bounding boxes is not robust in that case
  // and there is no easy and obvious fix for the moment. Keep this here
  // as a reminder that we still need to address that.

  // deallog.push("2d/3d");
  // test<2, 3>();
  // deallog.pop();
  deallog.push("3d/3d");
  test<3, 3>();
  deallog.pop();
}
