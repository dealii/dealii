// ---------------------------------------------------------------------
//
// Copyright (C) 2017 - 2019 by the deal.II authors
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



// check the particle sort within the particle handler class in serial models.

#include <deal.II/distributed/tria.h>

#include <deal.II/fe/mapping_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>

#include <deal.II/particles/particle_handler.h>

#include "../tests.h"

template <int dim, int spacedim>
void
test()
{
  {
    parallel::distributed::Triangulation<dim, spacedim> tr(MPI_COMM_WORLD);

    GridGenerator::hyper_cube(tr);
    tr.refine_global(1);
    MappingQ<dim, spacedim> mapping(1);

    Particles::ParticleHandler<dim, spacedim> particle_handler(tr, mapping);

    std::vector<Point<spacedim>> position(2);
    std::vector<Point<dim>>      reference_position(2);

    for (unsigned int i = 0; i < dim; ++i)
      {
        position[0](i) = 0.25;
        position[1](i) = 0.75;
      }

    Particles::Particle<dim, spacedim> particle1(position[0],
                                                 reference_position[0],
                                                 0);
    Particles::Particle<dim, spacedim> particle2(position[1],
                                                 reference_position[1],
                                                 1);

    typename Triangulation<dim, spacedim>::active_cell_iterator cell1(&tr,
                                                                      1,
                                                                      0);
    typename Triangulation<dim, spacedim>::active_cell_iterator cell2(&tr,
                                                                      1,
                                                                      0);

    particle_handler.insert_particle(particle1, cell1);
    particle_handler.insert_particle(particle2, cell2);

    for (const auto &particle : particle_handler)
      deallog << "Before sort particle id " << particle.get_id()
              << " is in cell " << particle.get_surrounding_cell(tr)
              << std::endl;

    particle_handler.sort_particles_into_subdomains_and_cells();

    for (const auto &particle : particle_handler)
      deallog << "After sort particle id " << particle.get_id()
              << " is in cell " << particle.get_surrounding_cell(tr)
              << std::endl;

    // Move all points up by 0.5. This will change cell for particle 1, and will
    // move particle 2 out of the domain. Note that we need to change the
    // coordinate dim-1 despite having a spacedim point.
    Point<spacedim> shift;
    shift(dim - 1) = 0.5;
    for (auto &particle : particle_handler)
      particle.set_location(particle.get_location() + shift);

    particle_handler.sort_particles_into_subdomains_and_cells();
    for (const auto &particle : particle_handler)
      deallog << "After shift particle id " << particle.get_id()
              << " is in cell " << particle.get_surrounding_cell(tr)
              << std::endl;
  }

  deallog << "OK" << std::endl;
}



int
main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

  initlog();

  deallog.push("2d/2d");
  test<2, 2>();
  deallog.pop();
  deallog.push("2d/3d");
  test<2, 3>();
  deallog.pop();
  deallog.push("3d/3d");
  test<3, 3>();
  deallog.pop();
}
