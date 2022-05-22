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

// Test constructing an rtree of particles from a ParticleHandler object.

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>

#include <deal.II/numerics/rtree.h>

#include <deal.II/particles/particle.h>
#include <deal.II/particles/particle_handler.h>

#include "../tests.h"

namespace bgi = boost::geometry::index;

template <int dim, int spacedim>
void
test()
{
  Triangulation<dim, spacedim> tria;
  GridGenerator::hyper_cube(tria);
  tria.refine_global(2);

  Particles::ParticleHandler<dim, spacedim> particle_handler(
    tria, StaticMappingQ1<dim, spacedim>::mapping);

  const int                    n_particles = 10;
  std::vector<Point<spacedim>> particles(n_particles);

  for (auto &p : particles)
    p = random_point<spacedim>();

  particle_handler.insert_particles(particles);

  auto tree = pack_rtree(particle_handler.begin(), particle_handler.end());

  auto p = random_point<spacedim>();
  for (const auto &part : tree | bgi::adaptors::queried(bgi::nearest(p, 3)))
    deallog << "Particle " << part.get_id() << " is close to " << p
            << " (location = " << part.get_location() << ')' << std::endl;
}

int
main()
{
  initlog();
  test<2, 2>();
}
