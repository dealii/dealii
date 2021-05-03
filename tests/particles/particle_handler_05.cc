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



// Check that particles are correctly transferred during grid
// refinement/coarsening. Note that this also tests some special cases. Since
// many particles lie exactly on cell boundaries they are sorted into the first
// child cell that approximately contains them (i.e the order in which the
// children are checked plays a role when deciding which child the particles is
// sorted into.

#include <deal.II/distributed/tria.h>

#include <deal.II/fe/mapping_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>

#include <deal.II/particles/particle_handler.h>

#include "../tests.h"

template <int dim, int spacedim>
void
create_regular_particle_distribution(
  Particles::ParticleHandler<dim, spacedim> &                particle_handler,
  const parallel::distributed::Triangulation<dim, spacedim> &tr,
  const unsigned int particles_per_direction = 3)
{
  for (unsigned int i = 0; i < particles_per_direction; ++i)
    for (unsigned int j = 0; j < particles_per_direction; ++j)
      {
        Point<spacedim> position;
        Point<dim>      reference_position;
        unsigned int    id = i * particles_per_direction + j;

        position[0] = static_cast<double>(i) /
                      static_cast<double>(particles_per_direction - 1);
        position[1] = static_cast<double>(j) /
                      static_cast<double>(particles_per_direction - 1);

        if (dim > 2)
          for (unsigned int k = 0; k < particles_per_direction; ++k)
            {
              position[2] = static_cast<double>(j) /
                            static_cast<double>(particles_per_direction - 1);
              id = i * particles_per_direction * particles_per_direction +
                   j * particles_per_direction + k;
              Particles::Particle<dim, spacedim> particle(position,
                                                          reference_position,
                                                          id);

              typename parallel::distributed::Triangulation<dim, spacedim>::
                active_cell_iterator cell =
                  GridTools::find_active_cell_around_point(
                    tr, particle.get_location());

              particle_handler.insert_particle(particle, cell);
            }
        else
          {
            Particles::Particle<dim, spacedim> particle(position,
                                                        reference_position,
                                                        id);

            typename parallel::distributed::Triangulation<dim, spacedim>::
              active_cell_iterator cell =
                GridTools::find_active_cell_around_point(
                  tr, particle.get_location());

            particle_handler.insert_particle(particle, cell);
          }
      }
}

template <int dim, int spacedim>
void
test()
{
  {
    parallel::distributed::Triangulation<dim, spacedim> tr(MPI_COMM_WORLD);

    GridGenerator::hyper_cube(tr);
    MappingQ<dim, spacedim> mapping(1);

    Particles::ParticleHandler<dim, spacedim> particle_handler(tr, mapping);

    // TODO: Move this into the Particle handler class. Unfortunately, there are
    // some interactions with the SolutionTransfer class that prevent us from
    // doing this at the moment. When doing this, check that transferring a
    // solution and particles during the same refinement is possible (in
    // particular that the order of serialization/deserialization is preserved).
    tr.signals.pre_distributed_refinement.connect(std::bind(
      &Particles::ParticleHandler<dim,
                                  spacedim>::register_store_callback_function,
      &particle_handler));

    tr.signals.post_distributed_refinement.connect(std::bind(
      &Particles::ParticleHandler<dim,
                                  spacedim>::register_load_callback_function,
      &particle_handler,
      false));

    create_regular_particle_distribution(particle_handler, tr);

    for (const auto &particle : particle_handler)
      deallog << "Before refinement particle id " << particle.get_id()
              << " is in cell " << particle.get_surrounding_cell(tr)
              << std::endl;

    // Check that all particles are moved to children
    tr.refine_global(1);

    for (const auto &particle : particle_handler)
      deallog << "After refinement particle id " << particle.get_id()
              << " is in cell " << particle.get_surrounding_cell(tr)
              << std::endl;

    // Reverse the refinement and check again
    for (auto cell = tr.begin_active(); cell != tr.end(); ++cell)
      cell->set_coarsen_flag();

    tr.execute_coarsening_and_refinement();

    for (const auto &particle : particle_handler)
      deallog << "After coarsening particle id " << particle.get_id()
              << " is in cell " << particle.get_surrounding_cell(tr)
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
