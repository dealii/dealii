// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2021 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

// Test the signal that is called if a particle is lost.

#include <deal.II/fe/mapping_q1.h>

#include <deal.II/grid/grid_generator.h>

#include <deal.II/particles/particle_handler.h>

#include "../tests.h"


int
main()
{
  initlog();

  const int dim = 2;

  Triangulation<dim> tria;
  GridGenerator::hyper_cube(tria);

  MappingQ1<dim> mapping;

  Particles::ParticleHandler<dim> particle_handler(tria, mapping);

  std::vector<Point<dim>> particle_positions;
  particle_positions.emplace_back(0.5, 0.5);
  particle_handler.insert_particles(particle_positions);

  particle_handler.signals.particle_lost.connect(
    [](const auto &particle, const auto &cell) {
      deallog << particle->get_id() << ' ' << cell->id() << std::endl;
    });

  particle_handler.begin()->set_location(Point<dim>{2.0, 2.0});

  particle_handler.sort_particles_into_subdomains_and_cells();
}
