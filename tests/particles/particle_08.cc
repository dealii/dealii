// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2019 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// Test set/get id.

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
