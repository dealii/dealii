// ---------------------------------------------------------------------
//
// Copyright (C) 2019 - 2022 by the deal.II authors
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
