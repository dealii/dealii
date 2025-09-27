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



// Like particle_02, but tests particle serialization and deserialization using
// boost archive.

#include <deal.II/base/array_view.h>

#include <deal.II/particles/particle.h>

#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>

#include "../tests.h"

template <int dim, int spacedim>
void
test()
{
  {
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

    Particles::Particle<dim, spacedim> particle(position,
                                                reference_position,
                                                index);

    deallog << "Particle location: " << particle.get_location() << std::endl
            << "Particle reference location: "
            << particle.get_reference_location() << std::endl
            << "Particle index: " << particle.get_id() << std::endl;

    std::stringstream             stream;
    boost::archive::text_oarchive archive(stream);

    archive << particle;

    Particles::Particle<dim, spacedim> new_particle;

    boost::archive::text_iarchive iarchive(stream);
    iarchive >> new_particle;

    deallog << "Copy particle location: " << new_particle.get_location()
            << std::endl
            << "Copy particle reference location: "
            << new_particle.get_reference_location() << std::endl
            << "Copy particle index: " << new_particle.get_id() << std::endl;
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
