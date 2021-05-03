// ---------------------------------------------------------------------
//
// Copyright (C) 2017 - 2018 by the deal.II authors
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



// Like particle_03, but tests the creation and use of a
// particle iterator from the created particle.

#include <deal.II/base/array_view.h>

#include <deal.II/particles/particle.h>
#include <deal.II/particles/particle_iterator.h>

#include "../tests.h"


template <int dim>
void
test()
{
  {
    const unsigned int           n_properties_per_particle = 3;
    Particles::PropertyPool<dim> pool(n_properties_per_particle);

    Point<dim> position;
    position(0) = 0.3;

    if (dim > 1)
      position(1) = 0.5;

    Point<dim> reference_position;
    reference_position(0) = 0.2;

    if (dim > 1)
      reference_position(1) = 0.4;

    const types::particle_index index(7);

    std::vector<double> properties = {0.15, 0.45, 0.75};

    Particles::Particle<dim> particle(position, reference_position, index);
    particle.set_property_pool(pool);
    particle.set_properties(
      ArrayView<double>(&properties[0], properties.size()));

    std::multimap<Particles::internal::LevelInd, Particles::Particle<dim>> map;

    Particles::internal::LevelInd level_index = std::make_pair(0, 0);
    map.insert(std::make_pair(level_index, particle));

    particle.get_properties()[0] = 0.05;
    map.insert(std::make_pair(level_index, particle));

    Particles::ParticleIterator<dim> particle_it(map, map.begin());
    Particles::ParticleIterator<dim> particle_end(map, map.end());

    for (; particle_it != particle_end; ++particle_it)
      {
        deallog << "Particle position: " << (*particle_it).get_location()
                << std::endl
                << "Particle properties: "
                << std::vector<double>(particle_it->get_properties().begin(),
                                       particle_it->get_properties().end())
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
