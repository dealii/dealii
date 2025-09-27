// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2019 - 2022 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// Test constructing an rtree of particles.

#include <deal.II/numerics/rtree.h>

#include <deal.II/particles/particle.h>

#include "../tests.h"

namespace bgi = boost::geometry::index;

template <int dim, int spacedim>
void
test()
{
  const int                                       n_particles = 10;
  std::vector<Particles::Particle<dim, spacedim>> particles(n_particles);

  unsigned int id = 0;
  for (auto &p : particles)
    {
      p.set_location(random_point<spacedim>());
      p.set_id(id++);
    }

  auto tree = pack_rtree(particles);

  auto p = random_point<spacedim>();
  for (const auto part : tree | bgi::adaptors::queried(bgi::nearest(p, 3)))
    deallog << "Particle " << part.get_id() << " is close to " << p
            << " (location = " << part.get_location() << ')' << std::endl;
}

int
main()
{
  initlog();
  test<2, 2>();
}
