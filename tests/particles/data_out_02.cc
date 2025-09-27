// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2017 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// Create a particle handler then generate an output
// Tests that scalar particle properties are output correctly.

#include <deal.II/base/data_out_base.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/fe/mapping_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>

#include <deal.II/particles/data_out.h>
#include <deal.II/particles/particle_handler.h>

#include "../tests.h"

template <int dim, int spacedim>
void
create_regular_particle_distribution(
  Particles::ParticleHandler<dim, spacedim> &particle_handler,
  const Triangulation<dim, spacedim>        &tr,
  const unsigned int                         particles_per_direction = 3)
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

              typename Triangulation<dim, spacedim>::active_cell_iterator cell =
                GridTools::find_active_cell_around_point(
                  tr, particle.get_location());

              Particles::ParticleIterator<dim, spacedim> pit =
                particle_handler.insert_particle(particle, cell);

              for (unsigned int i = 0; i < spacedim; ++i)
                pit->get_properties()[i] = pit->get_location()[i];
            }
        else
          {
            Particles::Particle<dim, spacedim> particle(position,
                                                        reference_position,
                                                        id);

            typename Triangulation<dim, spacedim>::active_cell_iterator cell =
              GridTools::find_active_cell_around_point(tr,
                                                       particle.get_location());

            Particles::ParticleIterator<dim, spacedim> pit =
              particle_handler.insert_particle(particle, cell);

            for (unsigned int i = 0; i < spacedim; ++i)
              pit->get_properties()[i] = pit->get_location()[i];
          }
      }
}

template <int dim, int spacedim>
void
test()
{
  {
    Triangulation<dim, spacedim> tr;

    GridGenerator::hyper_cube(tr);
    tr.refine_global(2);
    MappingQ<dim, spacedim> mapping(1);

    // Create a particle handler using two manually created particles
    Particles::ParticleHandler<dim, spacedim> particle_handler(tr,
                                                               mapping,
                                                               spacedim);

    create_regular_particle_distribution(particle_handler, tr, 2);

    std::vector<std::string> data_names;
    std::vector<DataComponentInterpretation::DataComponentInterpretation>
      data_interpretations;

    for (unsigned int i = 0; i < spacedim; ++i)
      {
        data_names.emplace_back("property_coord_" + std::to_string(i));
        data_interpretations.emplace_back(
          DataComponentInterpretation::component_is_scalar);
      }

    Particles::DataOut<dim, spacedim> particle_output;
    particle_output.build_patches(particle_handler,
                                  data_names,
                                  data_interpretations);
    particle_output.write_gnuplot(deallog.get_file_stream());
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
  deallog.push("2d/3d");
  test<2, 3>();
  deallog.pop();
  deallog.push("3d/3d");
  test<3, 3>();
  deallog.pop();
}
