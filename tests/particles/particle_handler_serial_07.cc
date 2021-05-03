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
  {
    Triangulation<dim, spacedim> tr;

    GridGenerator::hyper_cube(tr);
    tr.refine_global(2);

    MappingQ<dim, spacedim> mapping(1);

    Particles::ParticleHandler<dim, spacedim> particle_handler(tr, mapping);

    std::vector<Point<spacedim>> points(10);
    for (unsigned int i = 0; i < 10; ++i)
      {
        const double coordinate = static_cast<double>(i) / 10.0;
        for (unsigned int j = 0; j < spacedim; ++j)
          points[i][j] = 0.05 + coordinate;
      }

    particle_handler.insert_particles(points);
    particle_handler.update_cached_numbers();

    deallog << "Particle number: " << particle_handler.n_global_particles()
            << std::endl;

    for (auto particle = particle_handler.begin();
         particle != particle_handler.end();
         ++particle)
      deallog << "Particle id " << particle->get_id() << " is in cell "
              << particle->get_surrounding_cell(tr) << std::endl
              << "     at location " << particle->get_location() << std::endl
              << "     at reference location "
              << particle->get_reference_location() << std::endl;
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
