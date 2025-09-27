// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2020 - 2022 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

// Test constructing an rtree of particles from a ParticleHandler object.

#include <deal.II/fe/mapping_q1.h>

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
