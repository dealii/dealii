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



// check the removal of a single particle works as expected

#include <deal.II/distributed/tria.h>

#include <deal.II/fe/mapping_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>

#include <deal.II/particles/generators.h>
#include <deal.II/particles/particle_handler.h>

#include "../tests.h"

template <int dim, int spacedim>
void
test()
{
  {
    Triangulation<dim, spacedim> tr;

    GridGenerator::hyper_cube(tr);
    tr.refine_global(1);
    MappingQ<dim, spacedim> mapping(1);

    Particles::ParticleHandler<dim, spacedim> particle_handler(tr, mapping);

    std::vector<Point<dim>> particle_reference_locations(3, Point<dim>());

    for (unsigned int i = 0; i < dim; ++i)
      {
        particle_reference_locations[0][i] = 0.25;
        particle_reference_locations[1][i] = 0.5;
        particle_reference_locations[2][i] = 0.75;
      }

    Particles::Generators::regular_reference_locations(
      tr, particle_reference_locations, particle_handler);

    deallog << "Particle number: " << particle_handler.n_global_particles()
            << std::endl;

    for (const auto &cell : tr.active_cell_iterators())
      {
        const unsigned int n_particles_to_remove =
          cell->active_cell_index() % 4;

        for (unsigned int i = 0; i < n_particles_to_remove; ++i)
          {
            const unsigned int particle_index_to_remove =
              Testing::rand() % particle_handler.n_particles_in_cell(cell);
            auto particle_to_remove =
              particle_handler.particles_in_cell(cell).begin();
            std::advance(particle_to_remove, particle_index_to_remove);

            deallog << "Removing particle index: "
                    << particle_to_remove->get_id()
                    << ". Advanced by: " << particle_index_to_remove
                    << std::endl;

            particle_handler.remove_particle(particle_to_remove);
          }

        deallog << "Cell index: " << cell->active_cell_index()
                << ". N Particles: "
                << particle_handler.n_particles_in_cell(cell) << std::endl;
      }

    for (const auto &particle : particle_handler)
      {
        deallog << "Particle id: " << particle.get_id() << std::endl;
      }
  }

  deallog << "OK" << std::endl;
}



int
main(int argc, char *argv[])
{
  initlog();

  deallog.push("2d/2d");
  test<2, 2>();
  deallog.pop();
  deallog.push("3d/3d");
  test<3, 3>();
  deallog.pop();
}
