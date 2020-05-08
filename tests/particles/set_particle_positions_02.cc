// ---------------------------------------------------------------------
//
// Copyright (C) 2017 - 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE.md at
// the top level directory of deal.II.
//
// ---------------------------------------------------------------------

// Test insert_global_particles. Make sure we don't lose particles
// along the way

#include <deal.II/base/mpi.h>
#include <deal.II/base/utilities.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/fe/mapping_q.h>

#include <deal.II/grid/filtered_iterator.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>

#include <deal.II/particles/particle_handler.h>

#include "../tests.h"

template <int dim, int spacedim>
void
test()
{
  parallel::distributed::Triangulation<dim, spacedim> tr(MPI_COMM_WORLD);

  GridGenerator::hyper_cube(tr);
  tr.refine_global(2);

  MappingQ<dim, spacedim> mapping(1);

  Particles::ParticleHandler<dim, spacedim> particle_handler(tr, mapping);

  const unsigned int n_points = 2;
  const unsigned int my_cpu = Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);
  const unsigned int n_cpus = Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD);
  Testing::srand(my_cpu + 1);

  // Create the bounding boxes for the global insertion
  auto my_bounding_box = GridTools::compute_mesh_predicate_bounding_box(
    tr, IteratorFilters::LocallyOwnedCell());

  auto global_bounding_boxes =
    Utilities::MPI::all_gather(MPI_COMM_WORLD, my_bounding_box);

  std::vector<Point<spacedim>> points(n_points);
  for (auto &p : points)
    p = random_point<spacedim>(0.1, 0.2);

  auto cpu_to_index =
    particle_handler.insert_global_particles(points, global_bounding_boxes);

  for (const auto &particle : particle_handler)
    {
      deallog << "Particle id : " << particle.get_id();
      deallog << " Particle location: " << particle.get_location() << std::endl;
    }

  Tensor<1, spacedim> displacement;
  displacement[0] = 0.1;
  std::vector<Point<spacedim>> new_particle_positions(
    particle_handler.n_locally_owned_particles());

  std::vector<Point<spacedim>> old_particle_positions(
    particle_handler.n_locally_owned_particles());
  particle_handler.get_particle_positions(old_particle_positions);

  for (unsigned int i = 0; i < old_particle_positions.size(); ++i)
    new_particle_positions[i] = old_particle_positions[i] + displacement;

  particle_handler.set_particle_positions(new_particle_positions, false);

  deallog << "Displacement" << std::endl;
  for (const auto &particle : particle_handler)
    {
      deallog << "Particle id : " << particle.get_id();
      deallog << " Particle location: " << particle.get_location() << std::endl;
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
  deallog << "---" << std::endl;
  deallog.push("3d/3d");
  test<3, 3>();
  deallog.pop();
}
