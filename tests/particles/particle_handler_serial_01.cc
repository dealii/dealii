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



// check the creation and destruction of particle within the particle handler
// class for a serial triangulation.

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
  Triangulation<dim, spacedim> tr;

  GridGenerator::hyper_cube(tr);
  MappingQ<dim, spacedim> mapping(1);

  Particles::ParticleHandler<dim, spacedim> particle_handler(tr, mapping);


  Point<spacedim> position;
  position(0) = 0.3;
  if (spacedim > 1)
    position(1) = 0.5;
  if (spacedim > 2)
    position(2) = 0.7;

  Point<dim> reference_position;
  reference_position(0) = 0.2;
  if (dim > 1)
    reference_position(1) = 0.4;
  if (dim > 2)
    reference_position(2) = 0.6;

  Particles::Particle<dim, spacedim> particle(position, reference_position, 7);
  deallog << "Particle location: " << particle.get_location() << std::endl;


  std::pair<typename parallel::distributed::Triangulation<dim, spacedim>::
              active_cell_iterator,
            Point<dim>>
    cell_position =
      GridTools::find_active_cell_around_point(mapping,
                                               tr,
                                               particle.get_location());

  particle_handler.insert_particle(particle, cell_position.first);
  particle_handler.update_cached_numbers();

  deallog << "Particle number: " << particle_handler.n_global_particles()
          << std::endl;

  for (const auto &particle : particle_handler)
    {
      deallog << "Particle location: " << particle.get_location() << std::endl;
      deallog << "Particle reference location: "
              << particle.get_reference_location() << std::endl;
    }

  deallog << "OK" << std::endl;
}



int
main(int argc, char *argv[])
{
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
