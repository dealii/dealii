// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2017 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// Like particle_iterator_01, but changes the location through the
// particle iterator.

#include <deal.II/base/array_view.h>

#include <deal.II/fe/mapping_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>

#include <deal.II/particles/particle.h>
#include <deal.II/particles/particle_handler.h>
#include <deal.II/particles/particle_iterator.h>

#include "../tests.h"


template <int dim>
void
test()
{
  {
    Triangulation<dim> tr;

    GridGenerator::hyper_cube(tr);
    MappingQ<dim>      mapping(1);
    const unsigned int n_properties_per_particle = 3;

    Particles::ParticleHandler<dim> particle_handler(tr,
                                                     mapping,
                                                     n_properties_per_particle);

    Point<dim> position;
    position[0] = 0.3;

    if (dim > 1)
      position[1] = 0.5;

    Point<dim> reference_position;
    reference_position[0] = 0.2;

    if (dim > 1)
      reference_position[1] = 0.4;

    const types::particle_index index(7);

    std::vector<double> properties = {0.15, 0.45, 0.75};

    Particles::Particle<dim> particle(position, reference_position, index);

    // Insert two identical particles
    Particles::ParticleIterator<dim> particle_it =
      particle_handler.insert_particle(particle, tr.begin());
    particle_it->set_properties(
      ArrayView<double>(&properties[0], properties.size()));

    particle_it = particle_handler.insert_particle(particle, tr.begin());
    particle_it->set_properties(
      ArrayView<double>(&properties[0], properties.size()));

    // Modify the location of the second particle
    particle_it->get_location()[0] = 0.05;

    for (const auto &particle : particle_handler)
      {
        deallog << "Particle position: " << particle.get_location() << std::endl
                << "Particle properties: "
                << std::vector<double>(particle.get_properties().begin(),
                                       particle.get_properties().end())
                << std::endl;
      }
  }

  deallog << "OK" << std::endl;
}



int
main()
{
  initlog();
  test<2>();
}
