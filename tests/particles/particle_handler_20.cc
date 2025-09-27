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



// like particle_handler_06, but after exchanging the ghost particles
// modifies the location of the particles and updates them through
// the update_ghosts mechanism. This test also adds a single property
// to the particles. This property is also modified and then updated
// through the update_ghosts mechanism.

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

    Particles::ParticleHandler<dim, spacedim> particle_handler(tr, mapping, 2);

    unsigned int    n_particles = 3;
    Point<spacedim> position;
    Point<dim>      reference_position;

    for (unsigned int p = 0; p < n_particles; ++p)
      {
        if (Utilities::MPI::this_mpi_process(tr.get_mpi_communicator()) == 0)
          {
            for (unsigned int i = 0; i < dim; ++i)
              position[i] = 0.410 + 0.01 * p;

            Particles::Particle<dim, spacedim> particle(
              position,
              reference_position,
              Utilities::MPI::this_mpi_process(tr.get_mpi_communicator()) *
                  n_particles +
                p);
            typename Triangulation<dim, spacedim>::active_cell_iterator cell =
              tr.begin_active();
            particle_handler.insert_particle(particle, cell);
          }
      }

    particle_handler.sort_particles_into_subdomains_and_cells();


    unsigned int counter = 0;
    // Set the properties of the particle to be a unique number
    for (auto particle = particle_handler.begin();
         particle != particle_handler.end();
         ++particle)
      {
        particle->get_properties()[0] =
          1000 +
          100 * Utilities::MPI::this_mpi_process(tr.get_mpi_communicator()) +
          10 * particle->get_id();
        particle->get_properties()[1] =
          2000 +
          100 * Utilities::MPI::this_mpi_process(tr.get_mpi_communicator()) +
          10 * particle->get_id();
        counter++;
      }


    particle_handler.exchange_ghost_particles(true);

    for (auto particle = particle_handler.begin();
         particle != particle_handler.end();
         ++particle)
      deallog << "Particle id : " << particle->get_id()
              << " location : " << particle->get_location()
              << " property : " << particle->get_properties()[0] << " and "
              << particle->get_properties()[1] << " is local on process : "
              << Utilities::MPI::this_mpi_process(tr.get_mpi_communicator())
              << std::endl;

    for (auto particle = particle_handler.begin_ghost();
         particle != particle_handler.end_ghost();
         ++particle)
      deallog << "Particle id : " << particle->get_id()
              << " location : " << particle->get_location()
              << " property : " << particle->get_properties()[0] << " and "
              << particle->get_properties()[1] << " is ghost on process : "
              << Utilities::MPI::this_mpi_process(tr.get_mpi_communicator())
              << std::endl;

    deallog << "Modifying particles positions and properties" << std::endl;

    // Modify the location of a single particle on processor 0 and update the
    // ghosts
    for (auto particle = particle_handler.begin();
         particle != particle_handler.end();
         ++particle)
      {
        auto location = particle->get_location();
        location[0] += 0.1;
        particle->get_properties()[0] += 10000;
        particle->get_properties()[1] += 10000;
        particle->set_location(location);
      }

    // Update the ghost particles
    particle_handler.update_ghost_particles();


    for (auto particle = particle_handler.begin();
         particle != particle_handler.end();
         ++particle)
      deallog << "Particle id : " << particle->get_id()
              << " location : " << particle->get_location()
              << " property : " << particle->get_properties()[0] << " and "
              << particle->get_properties()[1] << " is local on process : "
              << Utilities::MPI::this_mpi_process(tr.get_mpi_communicator())
              << std::endl;

    for (auto particle = particle_handler.begin_ghost();
         particle != particle_handler.end_ghost();
         ++particle)
      deallog << "Particle id : " << particle->get_id()
              << " location : " << particle->get_location()
              << " property : " << particle->get_properties()[0] << " and "
              << particle->get_properties()[1] << " is ghost on process : "
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
