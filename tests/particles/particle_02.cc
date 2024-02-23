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



// check Particle constructors, copy, and move operations.

#include <deal.II/particles/particle.h>

#include "../tests.h"


template <int dim>
void
test()
{
  {
    Point<2> position;
    position[0] = 0.3;
    position[1] = 0.5;

    Point<2> reference_position;
    reference_position[0] = 0.2;
    reference_position[1] = 0.4;

    const types::particle_index index(7);

    Particles::Particle<2> particle(position, reference_position, index);

    deallog << "Particle location: " << particle.get_location() << std::endl
            << "Particle reference location: "
            << particle.get_reference_location() << std::endl
            << "Particle index: " << particle.get_id() << std::endl;

    const Particles::Particle<2> copy(particle);

    deallog << "Copy particle location: " << copy.get_location() << std::endl
            << "Copy particle reference location: "
            << copy.get_reference_location() << std::endl
            << "Copy particle index: " << copy.get_id() << std::endl;

    const Particles::Particle<2> moved_particle(std::move(particle));

    deallog << "Moved particle location: " << moved_particle.get_location()
            << std::endl
            << "Moved particle reference location: "
            << moved_particle.get_reference_location() << std::endl
            << "Moved particle index: " << moved_particle.get_id() << std::endl;

    deallog << "Original particle location: " << particle.get_location()
            << std::endl
            << "Original particle reference location: "
            << particle.get_reference_location() << std::endl
            << "Original particle index: " << particle.get_id() << std::endl;
  }

  deallog << "OK" << std::endl;
}



int
main()
{
  initlog();
  test<2>();
}
