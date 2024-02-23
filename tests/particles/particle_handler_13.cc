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
// single mpi process, and where we add two properties per particle.

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

  Particles::ParticleHandler<dim, spacedim> particle_handler(tr, mapping, 2);

  const unsigned int n_points = 3;
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
  std::vector<std::vector<double>> properties(n_points,
                                              {my_cpu + 1000.0,
                                               my_cpu + 1000.0});

  for (unsigned int i = 0; i < n_points; ++i)
    properties[i][1] = i + 100;

  for (auto &p : points)
    p = random_point<spacedim>();

  auto cpu_to_index =
    particle_handler.insert_global_particles(points,
                                             global_bounding_boxes,
                                             properties);

  for (const auto &c : cpu_to_index)
    {
      deallog << "From cpu: " << c.first << " I got : ";
      c.second.print(deallog);
    }

  if (cpu_to_index.find(my_cpu) != cpu_to_index.end())
    cpu_to_index.erase(cpu_to_index.find(my_cpu));
  auto received = Utilities::MPI::some_to_some(MPI_COMM_WORLD, cpu_to_index);


  for (const auto &c : received)
    {
      deallog << "To cpu : " << c.first << " I sent : ";
      c.second.print(deallog);
    }


  for (auto p : particle_handler)
    {
      deallog << "Particle : " << p.get_id() << ", properties: "
              << static_cast<unsigned int>(p.get_properties()[0]) << " - "
              << static_cast<unsigned int>(p.get_properties()[1]) << std::endl;
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
}
