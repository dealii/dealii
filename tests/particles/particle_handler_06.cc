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



// like particle_handler_04, but tests the exchange of ghost particles.

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

    // both processes create a particle handler with a particle in a
    // position that is in a ghost cell to the other process
    Particles::ParticleHandler<dim, spacedim> particle_handler(tr, mapping);

    Point<spacedim> position;
    Point<dim>      reference_position;

    if (Utilities::MPI::this_mpi_process(tr.get_mpi_communicator()) == 0)
      for (unsigned int i = 0; i < dim; ++i)
        position[i] = 0.475;
    else
      for (unsigned int i = 0; i < dim; ++i)
        position[i] = 0.525;

    Particles::Particle<dim, spacedim> particle(
      position,
      reference_position,
      Utilities::MPI::this_mpi_process(tr.get_mpi_communicator()));

    // We give a local random cell hint to check that sorting and
    // transferring ghost particles works.
    typename Triangulation<dim, spacedim>::active_cell_iterator cell =
      tr.begin_active();
    while (!cell->is_locally_owned())
      ++cell;

    particle_handler.insert_particle(particle, cell);

    particle_handler.sort_particles_into_subdomains_and_cells();
    particle_handler.exchange_ghost_particles();

    for (auto particle = particle_handler.begin();
         particle != particle_handler.end();
         ++particle)
      deallog << "Particle id " << particle->get_id()
              << " is local particle on process "
              << Utilities::MPI::this_mpi_process(tr.get_mpi_communicator())
              << std::endl;

    for (auto particle = particle_handler.begin_ghost();
         particle != particle_handler.end_ghost();
         ++particle)
      deallog << "Particle id " << particle->get_id()
              << " is ghost particle on process "
              << Utilities::MPI::this_mpi_process(tr.get_mpi_communicator())
              << std::endl;
  }

  deallog << "OK" << std::endl;
}



int
main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

  MPILogInitAll all;

  deallog.push("2d/2d");
  test<2, 2>();
  deallog.pop();
  deallog.push("2d/3d");
  test<2, 3>();
  deallog.pop();
  deallog.push("3d/3d");
  test<3, 3>();
  deallog.pop();
}
