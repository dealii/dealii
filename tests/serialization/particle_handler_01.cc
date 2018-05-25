// ---------------------------------------------------------------------
//
// Copyright (C) 2010 - 2017 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------


// check and illustrate the serialization process for ParticleHandler

#include <deal.II/distributed/tria.h>

#include <deal.II/fe/mapping_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>

#include <deal.II/particles/particle_handler.h>

#include "serialization.h"

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
              Particles::Particle<dim, spacedim> particle(
                position, reference_position, id);

              typename parallel::distributed::Triangulation<dim, spacedim>::
                active_cell_iterator cell =
                  GridTools::find_active_cell_around_point(
                    tr, particle.get_location());

              particle_handler.insert_particle(particle, cell);
            }
        else
          {
            Particles::Particle<dim, spacedim> particle(
              position, reference_position, id);

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
  parallel::distributed::Triangulation<dim, spacedim> tr(MPI_COMM_WORLD);
  GridGenerator::hyper_cube(tr);
  tr.refine_global(2);
  MappingQ<dim, spacedim> mapping(1);

  Particles::ParticleHandler<dim, spacedim> particle_handler(tr, mapping);

  create_regular_particle_distribution(particle_handler, tr);
  particle_handler.sort_particles_into_subdomains_and_cells();

  for (auto particle = particle_handler.begin();
       particle != particle_handler.end();
       ++particle)
    deallog << "Before serialization particle id " << particle->get_id()
            << " is in cell " << particle->get_surrounding_cell(tr)
            << std::endl;

  // TODO: Move this into the Particle handler class. Unfortunately, there are
  // some interactions with the SolutionTransfer class that prevent us from
  // doing this at the moment. When doing this, check that transferring a
  // solution and particles during the same refinement is possible (in
  // particular that the order of serialization/deserialization is preserved).
  tr.signals.pre_distributed_save.connect(std::bind(
    &Particles::ParticleHandler<dim,
                                spacedim>::register_store_callback_function,
    &particle_handler,
    true));

  // save data to archive
  std::ostringstream oss;
  {
    boost::archive::text_oarchive oa(oss, boost::archive::no_header);
    oa << particle_handler;
    tr.save("checkpoint");

    // archive and stream closed when
    // destructors are called
  }
  deallog << oss.str() << std::endl;

  // Now remove all information in tr and particle_handler,
  // this is like creating new objects after a restart
  tr.clear();
  GridGenerator::hyper_cube(tr);

  particle_handler.clear();
  particle_handler.initialize(tr, mapping);

  // This should not produce any output
  for (auto particle = particle_handler.begin();
       particle != particle_handler.end();
       ++particle)
    deallog << "In between particle id " << particle->get_id() << " is in cell "
            << particle->get_surrounding_cell(tr) << std::endl;


  // TODO: Move this into the Particle handler class. Unfortunately, there are
  // some interactions with the SolutionTransfer class that prevent us from
  // doing this at the moment. When doing this, check that transferring a
  // solution and particles during the same refinement is possible (in
  // particular that the order of serialization/deserialization is preserved).
  tr.signals.post_distributed_load.connect(std::bind(
    &Particles::ParticleHandler<dim, spacedim>::register_load_callback_function,
    &particle_handler,
    true));

  // verify correctness of the serialization. Note that the deserialization of
  // the particle handler has to happen before the triangulation (otherwise it
  // does not know if something was stored in the user data of the
  // triangulation).
  {
    std::istringstream            iss(oss.str());
    boost::archive::text_iarchive ia(iss, boost::archive::no_header);

    ia >> particle_handler;
    tr.load("checkpoint");
  }

  for (auto particle = particle_handler.begin();
       particle != particle_handler.end();
       ++particle)
    deallog << "After serialization particle id " << particle->get_id()
            << " is in cell " << particle->get_surrounding_cell(tr)
            << std::endl;

  deallog << "OK" << std::endl << std::endl;
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
