// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2018 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// Like particle_04, but also attach properties.

#include <deal.II/base/array_view.h>

#include <deal.II/particles/particle.h>

#include "../tests.h"


template <int dim, int spacedim>
void
test()
{
  {
    const unsigned int                     n_properties_per_particle = 3;
    Particles::PropertyPool<dim, spacedim> pool(n_properties_per_particle);

    Point<spacedim> position;

    position[0] = 0.3;
    if (spacedim > 1)
      position[1] = 0.5;
    if (spacedim > 2)
      position[2] = 0.7;

    Point<dim> reference_position;
    reference_position[0] = 0.2;
    if (dim > 1)
      reference_position[1] = 0.4;
    if (dim > 2)
      reference_position[2] = 0.6;

    const types::particle_index index(7);

    std::vector<double> properties = {0.15, 0.45, 0.75};

    Particles::Particle<dim, spacedim> particle(position,
                                                reference_position,
                                                index);
    particle.set_property_pool(pool);
    particle.set_properties(
      ArrayView<double>(&properties[0], properties.size()));

    deallog << "Particle location: " << particle.get_location() << std::endl
            << "Particle reference location: "
            << particle.get_reference_location() << std::endl
            << "Particle index: " << particle.get_id() << std::endl
            << "Particle properties: "
            << std::vector<double>(particle.get_properties().begin(),
                                   particle.get_properties().end())
            << std::endl;

    std::vector<char> data(particle.serialized_size_in_bytes());
    void             *write_pointer = static_cast<void *>(&data.front());

    write_pointer = particle.write_particle_data_to_memory(write_pointer);

    const void *read_pointer = static_cast<const void *>(&data.front());
    const Particles::Particle<dim, spacedim> new_particle(read_pointer, &pool);

    deallog << "Copy particle location: " << new_particle.get_location()
            << std::endl
            << "Copy particle reference location: "
            << new_particle.get_reference_location() << std::endl
            << "Copy particle index: " << new_particle.get_id() << std::endl
            << "Copy particle properties: "
            << std::vector<double>(new_particle.get_properties().begin(),
                                   new_particle.get_properties().end())
            << std::endl;
  }

  deallog << "OK" << std::endl;
}



int
main()
{
  initlog();

  test<1, 1>();
  test<1, 2>();
  test<1, 3>();

  test<2, 2>();
  test<2, 3>();

  test<3, 3>();
}
