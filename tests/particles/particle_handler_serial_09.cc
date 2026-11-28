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

// This test verifies the correctness of retrieving all particles that belong to
// a specified parent cell using
//   Particles::ParticleHandler::particles_in_descendant_active_cells() and
//   Particles::ParticleHandler::particles_in_active_subtrees_of_level_cells().
//
// The test setup consists of a square domain that is refined three times. A
// fixed number of particles are inserted at random positions in the domain.
//
// For each cell on level 1, the test selects one of its descendant active cells
// on level 3 and queries all particles that belong to the corresponding parent
// cell on level 1. The retrieved particle IDs are printed. In addition, the
// test calls
//   Particles::ParticleHandler::particles_in_active_subtrees_of_level_cells()
// to obtain a map from level-1 cells to the particles they contain, and prints
// the particle IDs for each cell.

#include <deal.II/fe/mapping_q1.h>

#include <deal.II/grid/cell_id.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>

#include <deal.II/particles/particle.h>
#include <deal.II/particles/particle_handler.h>

#include "../tests.h"


template <int dim, int spacedim>
void
test()
{
  using particle_iterator_vector =
    std::vector<typename Particles::ParticleHandler<dim, spacedim>::
                  particle_iterator_range>;

  using cell_iterator = typename Triangulation<dim, spacedim>::cell_iterator;

  Triangulation<dim, spacedim> tria;
  GridGenerator::hyper_cube(tria);
  tria.refine_global(3);

  Particles::ParticleHandler<dim, spacedim> particle_handler(
    tria, StaticMappingQ1<dim, spacedim>::mapping);

  const int                    n_particles = 10;
  std::vector<Point<spacedim>> particles(n_particles);

  for (auto &p : particles)
    p = random_point<spacedim>();

  particle_handler.insert_particles(particles);

  auto log_particle_ids = [](const particle_iterator_vector &particle_ranges) {
    for (const auto &particle_range : particle_ranges)
      for (const auto &particle : particle_range)
        {
          deallog << particle.get_id() << " ";
        }
  };

  for (const auto &level_1_cell_iterator : tria.cell_iterators_on_level(1))
    {
      cell_iterator active_cell = level_1_cell_iterator;
      while (active_cell->has_children())
        active_cell = active_cell->child(0);

      particle_iterator_vector particles_in_level_1_parent_cell =
        particle_handler.particles_in_descendant_active_cells(
          active_cell, active_cell->level() - 1);

      deallog << "Particles (IDs) in parent cell at level 1 (CellId: "
              << level_1_cell_iterator->id() << "): " << std::endl;
      log_particle_ids(particles_in_level_1_parent_cell);
      deallog << std::endl;
    }

  std::map<cell_iterator, particle_iterator_vector>
    all_particles_in_level_1_parent_cell_map =
      particle_handler.particles_in_active_subtrees_of_parent_cells(1);

  deallog << "Cells and Particles on Level 1: " << std::endl;
  for (const auto &[cell, particle_ranges] :
       all_particles_in_level_1_parent_cell_map)
    {
      deallog << "CellId: " << cell->id() << ", Particles (IDs): ";
      log_particle_ids(particle_ranges);
      deallog << std::endl;
    }
}

int
main()
{
  initlog();
  test<2, 2>();
}
