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



// Like particle_03, but test property initialization.

#include <deal.II/base/array_view.h>

#include <deal.II/particles/particle.h>

#include "../tests.h"


template <int dim>
void
test()
{
  {
    const unsigned int           n_properties_per_particle = 3;
    Particles::PropertyPool<dim> pool(n_properties_per_particle);

    Point<2> position;
    position[0] = 0.3;
    position[1] = 0.5;

    Point<2> reference_position;
    reference_position[0] = 0.2;
    reference_position[1] = 0.4;

    const types::particle_index index(7);

    std::vector<double> properties = {0.15, 0.45, 0.75};

    Particles::Particle<2> particle(position, reference_position, index);
    particle.set_property_pool(pool);

    for (unsigned int i = 0; i < n_properties_per_particle; ++i)
      AssertThrow(particle.get_properties()[i] == 0.0, ExcInternalError());
  }

  deallog << "OK" << std::endl;
}



int
main()
{
  initlog();
  test<2>();
}
