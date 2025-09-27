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



// like particle_handler_04, but for fully distributed triangulations in
// parallel computations

#include <deal.II/distributed/fully_distributed_tria.h>
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
    const MPI_Comm comm = MPI_COMM_WORLD;

    // create a distributed triangulation
    parallel::distributed::Triangulation<dim, spacedim> tria_pdt(comm);

    GridGenerator::hyper_cube(tria_pdt);
    tria_pdt.refine_global(2);
    MappingQ<dim, spacedim> mapping(1);


    // create a fully distributed triangulation
    parallel::fullydistributed::Triangulation<dim, spacedim> tria_pft(comm);

    // extract relevant information from distributed triangulation
    auto construction_data = TriangulationDescription::Utilities::
      create_description_from_triangulation(tria_pdt, comm);

    // actually create triangulation
    tria_pft.create_triangulation(construction_data);

    //    both processes create a particle handler,
    //      but only the first creates particles
    Particles::ParticleHandler<dim, spacedim> particle_handler(tria_pft,
                                                               mapping);

    if (Utilities::MPI::this_mpi_process(tria_pft.get_mpi_communicator()) == 0)
      {
        std::vector<Point<spacedim>> position(2);
        std::vector<Point<dim>>      reference_position(2);

        for (unsigned int i = 0; i < dim; ++i)
          {
            position[0][i] = 0.125;
            position[1][i] = 0.525;
          }

        Particles::Particle<dim, spacedim> particle1(position[0],
                                                     reference_position[0],
                                                     0);
        Particles::Particle<dim, spacedim> particle2(position[1],
                                                     reference_position[1],
                                                     1);

        typename Triangulation<dim, spacedim>::active_cell_iterator cell1(
          &tria_pft, 2, 0);
        typename Triangulation<dim, spacedim>::active_cell_iterator cell2(
          &tria_pft, 2, 0);

        particle_handler.insert_particle(particle1, cell1);
        particle_handler.insert_particle(particle2, cell2);

        for (const auto &particle : particle_handler)
          deallog << "Before sort particle id " << particle.get_id()
                  << " is in cell " << particle.get_surrounding_cell()
                  << " on process "
                  << Utilities::MPI::this_mpi_process(
                       tria_pft.get_mpi_communicator())
                  << std::flush << std::endl;
      }



    particle_handler.sort_particles_into_subdomains_and_cells();

    for (const auto &particle : particle_handler)
      deallog << "After sort particle id " << particle.get_id()
              << " is in cell " << particle.get_surrounding_cell()
              << " on process "
              << Utilities::MPI::this_mpi_process(
                   tria_pft.get_mpi_communicator())
              << std::flush << std::endl;

    // Move all points up by 0.5. This will change cell for particle 1 and will
    // move particle 2 out of the domain. Note that we need to change the
    // coordinate dim-1 despite having a spacedim point.
    Point<spacedim> shift;
    shift[dim - 1] = 0.5;
    for (auto &particle : particle_handler)
      particle.set_location(particle.get_location() + shift);

    particle_handler.sort_particles_into_subdomains_and_cells();
    for (const auto &particle : particle_handler)
      deallog << "After shift particle id " << particle.get_id()
              << " is in cell " << particle.get_surrounding_cell()
              << " on process "
              << Utilities::MPI::this_mpi_process(
                   tria_pft.get_mpi_communicator())
              << std::flush << std::endl;
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
