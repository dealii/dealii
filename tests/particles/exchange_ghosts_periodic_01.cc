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

// This test a triangulation which has a periodic boundary along the x axis.
// Particles are generated close to the x- and x+ limit of the domain
// and ghost particles are exchanged. Because the domain is periodic,
// ghost particles should be generated on both processors even if the cells
// do not share vertices because they are connected by the periodicity.
// This test verifies that ghost particles are indeed generated on both cores.

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

    // Create a subdivided hyper rectangle with divisions along the x axis
    // Colorize is set to true to enable setting one set of periodic boundary
    // conditions along id 0 and 1
    std::vector<unsigned int> repetitions({6, 1});
    if (dim == 3)
      repetitions.push_back(1);
    Point<spacedim> left;
    Point<spacedim> right;
    for (int i = 0; i < spacedim; ++i)
      right[i] = 1;

    GridGenerator::subdivided_hyper_rectangle(
      tr, repetitions, left, right, true);
    MappingQ<dim, spacedim> mapping(1);


    // Generate periodic boundary conditions
    std::vector<GridTools::PeriodicFacePair<
      typename Triangulation<dim, spacedim>::cell_iterator>>
      periodicity_vector;
    GridTools::collect_periodic_faces(tr, 0, 1, 0, periodicity_vector);
    tr.add_periodicity(periodicity_vector);

    // Both processes create a particle handler with a particle in a
    // position that is in a ghost cell to the other process connected by a
    // periodic boundary
    Particles::ParticleHandler<dim, spacedim> particle_handler(tr, mapping);

    Point<spacedim> position;
    Point<dim>      reference_position;

    if (Utilities::MPI::this_mpi_process(tr.get_mpi_communicator()) == 0)
      position[0] = 0.001;
    else
      position[0] = 0.999;

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

    // We count the number of ghost particles before they are exchanged
    // The expected result is zero

    unsigned int n_ghost_particles = 0;
    for (auto ghost_particle = particle_handler.begin_ghost();
         ghost_particle != particle_handler.end_ghost();
         ++ghost_particle)
      {
        ++n_ghost_particles;
      }

    deallog << "Number of ghost particles before exchange ghost: "
            << std::to_string(n_ghost_particles) << std::endl;

    particle_handler.exchange_ghost_particles();

    n_ghost_particles = 0;
    for (auto ghost_particle = particle_handler.begin_ghost();
         ghost_particle != particle_handler.end_ghost();
         ++ghost_particle)
      {
        ++n_ghost_particles;
      }

    // We count the number of ghost particles after they are exchanged
    // The expected result is one

    deallog << "Number of ghost particles after exchange ghost: "
            << std::to_string(n_ghost_particles) << std::endl;
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
  deallog.push("3d/3d");
  test<3, 3>();
  deallog.pop();
}
