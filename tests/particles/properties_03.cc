// ---------------------------------------------------------------------
//
// Copyright (C) 2020 by the deal.II authors
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



// When inserting multiple particles at once, we forgot to set their
// property pool pointer. This test differs from properties_01 by
// testing another function that forgot to set the property pool.

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/grid/grid_generator.h>

#include <deal.II/particles/particle_handler.h>
#include <deal.II/particles/utilities.h>

#include "../tests.h"

template <int dim>
void
test()
{
  // number of particles to place
  const unsigned int      n_particles = 10;
  std::vector<Point<dim>> particle_positions(n_particles);

  for (unsigned int i = 0; i < n_particles; ++i)
    {
      for (unsigned int d = 0; d < dim; ++d)
        particle_positions[i][d] =
          static_cast<double>(i) / static_cast<double>(n_particles);
    }

  const MappingQ1<dim> mapping;
  Triangulation<dim>   triangulation;
  GridGenerator::hyper_cube(triangulation, 0.0, 1.0);
  triangulation.refine_global(1);

  // create the particle handler
  const unsigned int              n_property_components = 1;
  Particles::ParticleHandler<dim> particle_handler(triangulation,
                                                   mapping,
                                                   n_property_components);

  particle_handler.insert_particles(particle_positions);

  // Now try to set the interpolated quantities as particle
  // properties. This used to fail before because we had forgotten to
  // the set the property_pool pointer of the particles that were
  // bulk-inserted before.
  for (auto &particle : particle_handler)
    {
      particle.set_properties(
        std::vector<double>(1, particle.get_location()[0]));
    }

  // Finally output these properties
  for (const auto particle : particle_handler)
    {
      deallog << "Particle " << particle.get_id() << ": ";
      for (const auto p : particle.get_properties())
        deallog << p << ' ';
      deallog << std::endl;
    }
}



int
main()
{
  initlog();
  test<1>();
  test<2>();
  test<3>();
}
