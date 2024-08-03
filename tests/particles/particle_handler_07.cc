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



// check the creation and destruction of particle within the particle handler
// class.

#include <deal.II/distributed/tria.h>

#include <deal.II/fe/mapping_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>

#include <deal.II/particles/particle_handler.h>

#include "../tests.h"

template <int dim, int spacedim>
void
test()
{
  {
    parallel::distributed::Triangulation<dim, spacedim> tr(MPI_COMM_WORLD);

    GridGenerator::hyper_cube(tr);
    tr.refine_global(2);

    MappingQ<dim, spacedim> mapping(1);

    Particles::ParticleHandler<dim, spacedim> particle_handler(tr, mapping);

    std::vector<Point<spacedim>> points(10);
    for (unsigned int i = 0; i < 10; ++i)
      {
        const double coordinate = static_cast<double>(i) / 10.0;
        for (unsigned int j = 0; j < spacedim; ++j)
          points[i][j] = 0.05 + coordinate;
      }

    particle_handler.insert_particles(points);
    particle_handler.update_cached_numbers();

    deallog << "Particle number: " << particle_handler.n_global_particles()
            << std::endl;

    for (auto particle = particle_handler.begin();
         particle != particle_handler.end();
         ++particle)
      deallog << "Particle id " << particle->get_id() << " is in cell "
              << particle->get_surrounding_cell() << std::endl
              << "     at location " << particle->get_location() << std::endl
              << "     at reference location "
              << particle->get_reference_location() << std::endl;
  }

  deallog << "OK" << std::endl;
}



int
main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

  initlog();

  deallog.push("2d/2d");
  test<2, 2>();
  deallog.pop();
  deallog.push("3d/3d");
  test<3, 3>();
  deallog.pop();
}
