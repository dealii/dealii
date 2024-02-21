// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2019 - 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

// Test insert_global_particles and set_particles_positions using a
// Function<spacedim>

#include <deal.II/base/mpi.h>
#include <deal.II/base/utilities.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/fe/mapping_q.h>

#include <deal.II/grid/filtered_iterator.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>

#include <deal.II/particles/particle_handler.h>

#include "../tests.h"

template <int spacedim>
class Positions : public Function<spacedim>
{
public:
  Positions(unsigned int n_components)
    : Function<spacedim>(n_components)
  {}

  virtual double
  value(const Point<spacedim> &p, unsigned int component = 0) const
  {
    if (component == 0)
      return p[component] + 0.2;
    else
      return p[component];
  }

private:
};

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

  Positions<spacedim> new_positions(spacedim);

  particle_handler.set_particle_positions(new_positions, false);

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
