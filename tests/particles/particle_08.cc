// ---------------------------------------------------------------------
//
// Copyright (C) 2019 by the deal.II authors
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



// Test set/get id.

#include <deal.II/particles/particle.h>

#include "../tests.h"


template <int dim>
void
test()
{
  {
    Point<2> position;
    position(0) = 0.3;
    position(1) = 0.5;

    Point<2> reference_position;
    reference_position(0) = 0.2;
    reference_position(1) = 0.4;

    const types::particle_index index(7);

    Particles::Particle<2> particle(position, reference_position, index);

    deallog << "ID: " << particle.get_id() << std::endl;
    particle.set_id(9);
    deallog << "New ID: " << particle.get_id() << std::endl;
  }
}

int
main()
{
  initlog();
  test<2>();
}
