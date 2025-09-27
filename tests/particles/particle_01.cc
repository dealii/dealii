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



// check the creation and destruction of particles

#include <deal.II/particles/particle.h>

#include "../tests.h"


template <int dim>
void
test()
{
  {
    Particles::Particle<dim> particle;

    Point<dim> position;
    position[0] = 1.0;
    particle.set_location(position);

    deallog << "Particle location: " << particle.get_location() << std::endl;
  }

  deallog << "OK" << std::endl;
}



int
main()
{
  initlog();
  test<2>();
  test<3>();
}
