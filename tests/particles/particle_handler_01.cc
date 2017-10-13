// ---------------------------------------------------------------------
//
// Copyright (C) 2017 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------



// check the creation and destruction of particle within the particle handler class.

#include "../tests.h"
#include <deal.II/particles/particle_handler.h>

#include <deal.II/distributed/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/fe/mapping_q.h>

template <int dim, int spacedim>
void test ()
{
  {
    parallel::distributed::Triangulation<dim,spacedim> tr(MPI_COMM_WORLD);

    GridGenerator::hyper_cube(tr);
    MappingQ<dim,spacedim> mapping(1);

    Particles::ParticleHandler<dim,spacedim> particle_handler(tr,mapping);

    Particles::Particle<dim,spacedim> particle;

    Point<spacedim> position;
    position(0) = 0.3;
    if (spacedim>1)
      position(1) = 0.5;
    if (spacedim>2)
      position(2) = 0.7;

    Point<dim> reference_position;
    reference_position(0) = 0.2;
    if (dim>1)
      reference_position(1) = 0.4;
    if (dim>2)
      reference_position(1) = 0.6;

    particle.set_location(position);
    particle.set_reference_location(reference_position);
    deallog << "Particle location: " << particle.get_location() << std::endl;


    std::pair<typename parallel::distributed::Triangulation<dim, spacedim>::active_cell_iterator, Point<dim> > cell_position =
      GridTools::find_active_cell_around_point(mapping, tr, particle.get_location());

    particle_handler.insert_particle(particle,cell_position.first);
    particle_handler.update_cache();

    deallog << "Particle number: " << particle_handler.n_global_particles() << std::endl;

    for (auto particle = particle_handler.begin(); particle != particle_handler.end(); ++particle)
      {
        deallog << "Particle location: " << particle->get_location() << std::endl;
        deallog << "Particle reference location: " << particle->get_reference_location() << std::endl;
      }
  }

  deallog << "OK" << std::endl;
}



int main (int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization (argc, argv, 1);

  initlog();

  test<2,2>();
  test<2,3>();

  test<3,3>();
}
