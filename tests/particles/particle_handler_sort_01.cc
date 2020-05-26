// ---------------------------------------------------------------------
//
// Copyright (C) 2017 - 2019 by the deal.II authors
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



// A test that uses the particle handler signals to get notified about
// lost particles. Due to the setup of the test it can easily be used
// to benchmark the performance of the
// ParticleHandler::sort_particles_into_subdomains_and_cells() function
// for different numbers of particles and application scenarios.

#include <deal.II/base/timer.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/fe/mapping_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>

#include <deal.II/particles/generators.h>
#include <deal.II/particles/particle_handler.h>

#include "../tests.h"

template <int dim, int spacedim>
void
lost_particle_notification(
  const typename Particles::ParticleIterator<dim, spacedim> &        particle,
  const typename Triangulation<dim, spacedim>::active_cell_iterator &cell)
{
  deallog << "Particle <" << particle->get_id()
          << "> lost. Current position: " << particle->get_location()
          << std::endl;
}

template <int dim, int spacedim>
void
test()
{
  {
    const unsigned int n_particles     = 1000;
    const double       max_coord_shift = 0.1;
    const bool         measure_time    = false;

    parallel::distributed::Triangulation<dim, spacedim> tr(MPI_COMM_WORLD);
    TimerOutput timer(std::cout, TimerOutput::summary, TimerOutput::wall_times);

    if (measure_time == true)
      timer.enter_subsection("Generate grid");

    GridGenerator::hyper_cube(tr);
    tr.refine_global(2);
    MappingQ<dim, spacedim> mapping(1);

    if (measure_time == true)
      timer.leave_subsection("Generate grid");

    if (measure_time == true)
      timer.enter_subsection("Generate particles");

    Particles::ParticleHandler<dim, spacedim> particle_handler(tr, mapping);
    particle_handler.signals.particle_lost.connect(
      [&](const typename Particles::ParticleIterator<dim, spacedim> &particle,
          const typename Triangulation<dim, spacedim>::active_cell_iterator
            &cell) { lost_particle_notification(particle, cell); });

    // Generate initial particle distribution
    Functions::ConstantFunction<spacedim> particle_density(1.0);
    Particles::Generators::probabilistic_locations(
      tr, particle_density, false, n_particles, particle_handler, mapping);

    if (measure_time == true)
      timer.leave_subsection("Generate particles");

    if (measure_time == true)
      timer.enter_subsection("Move particles");

    // Move particles by random distance up to max_distance
    for (auto &particle : particle_handler)
      {
        Point<spacedim> shift;

        for (unsigned int i = 0; i < dim; ++i)
          shift[i] = random_value(-max_coord_shift, max_coord_shift);

        particle.set_location(particle.get_location() + shift);
      }

    if (measure_time == true)
      timer.leave_subsection("Move particles");

    if (measure_time == true)
      timer.enter_subsection("Sort particles");

    // Measure sort particle function
    particle_handler.sort_particles_into_subdomains_and_cells();

    if (measure_time == true)
      timer.leave_subsection("Sort particles");
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
