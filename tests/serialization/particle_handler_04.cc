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


// check and illustrate the serialization process for ParticleHandler
// for fully distributed triangulations

#include <deal.II/distributed/fully_distributed_tria.h>

#include <deal.II/fe/mapping_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>

#include <deal.II/particles/particle_handler.h>

#include "serialization.h"

template <int dim, int spacedim>
void
create_regular_particle_distribution(
  Particles::ParticleHandler<dim, spacedim> &particle_handler,
  const parallel::fullydistributed::Triangulation<dim, spacedim> &tr,
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

              typename parallel::fullydistributed::
                Triangulation<dim, spacedim>::active_cell_iterator cell =
                  GridTools::find_active_cell_around_point(
                    tr, particle.get_location());

              particle_handler.insert_particle(particle, cell);
            }
        else
          {
            Particles::Particle<dim, spacedim> particle(position,
                                                        reference_position,
                                                        id);

            typename parallel::fullydistributed::Triangulation<dim, spacedim>::
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
  // Generate fulllydistributed triangulation from serial triangulation
  Triangulation<dim, spacedim> basetria;
  GridGenerator::hyper_cube(basetria);
  basetria.refine_global(2);

  auto construction_data =
    TriangulationDescription::Utilities::create_description_from_triangulation(
      basetria, MPI_COMM_WORLD);

  parallel::fullydistributed::Triangulation<dim, spacedim> tr(MPI_COMM_WORLD);
  tr.create_triangulation(construction_data);

  MappingQ<dim, spacedim> mapping(1);

  Particles::ParticleHandler<dim, spacedim> particle_handler(tr, mapping);

  create_regular_particle_distribution(particle_handler, tr);
  particle_handler.sort_particles_into_subdomains_and_cells();

  for (auto particle = particle_handler.begin();
       particle != particle_handler.end();
       ++particle)
    deallog << "Before serialization particle id " << particle->get_id()
            << " is in cell " << particle->get_surrounding_cell() << std::endl;

  // save data to archive
  std::ostringstream oss;
  {
    boost::archive::text_oarchive oa(oss, boost::archive::no_header);

    oa << particle_handler;
    particle_handler.prepare_for_serialization();
    tr.save("checkpoint");

    // archive and stream closed when
    // destructors are called
  }
  deallog << oss.str() << std::endl;

  // Now remove all information in tr and particle_handler,
  // this is like creating new objects after a restart
  tr.clear();

  particle_handler.clear();
  particle_handler.initialize(tr, mapping);

  // This should not produce any output
  for (auto particle = particle_handler.begin();
       particle != particle_handler.end();
       ++particle)
    deallog << "In between particle id " << particle->get_id() << " is in cell "
            << particle->get_surrounding_cell() << std::endl;

  // verify correctness of the serialization. Note that the deserialization of
  // the particle handler has to happen before the triangulation (otherwise it
  // does not know if something was stored in the user data of the
  // triangulation).
  {
    std::istringstream            iss(oss.str());
    boost::archive::text_iarchive ia(iss, boost::archive::no_header);

    ia >> particle_handler;
    tr.load("checkpoint");
    particle_handler.deserialize();
  }

  for (auto particle = particle_handler.begin();
       particle != particle_handler.end();
       ++particle)
    deallog << "After serialization particle id " << particle->get_id()
            << " is in cell " << particle->get_surrounding_cell() << std::endl;

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
