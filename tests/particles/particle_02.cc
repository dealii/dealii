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



// check Particle constructors, copy, and move operations.

#include "../tests.h"
#include <deal.II/particles/particle.h>


template <int dim>
void test ()
{
  {
    Point<2> position;
    position(0) = 0.3;
    position(1) = 0.5;

    Point<2> reference_position;
    reference_position(0) = 0.2;
    reference_position(1) = 0.4;

    const types::particle_index index(7);

    Particles::Particle<2> particle(position,reference_position,index);

    deallog << "Particle location: " << particle.get_location() << std::endl
            << "Particle reference location: " << particle.get_reference_location() << std::endl
            << "Particle index: " << particle.get_id() << std::endl;

    const Particles::Particle<2> copy(particle);

    deallog << "Copy particle location: " << copy.get_location() << std::endl
            << "Copy particle reference location: " << copy.get_reference_location() << std::endl
            << "Copy particle index: " << copy.get_id() << std::endl;

    const Particles::Particle<2> moved_particle(std::move(particle));

    deallog << "Moved particle location: " << moved_particle.get_location() << std::endl
            << "Moved particle reference location: " << moved_particle.get_reference_location() << std::endl
            << "Moved particle index: " << moved_particle.get_id() << std::endl;

    deallog << "Original particle location: " << particle.get_location() << std::endl
            << "Original particle reference location: " << particle.get_reference_location() << std::endl
            << "Original particle index: " << particle.get_id() << std::endl;
  }

  deallog << "OK" << std::endl;
}



int main ()
{
  initlog();
  test<2>();
}
