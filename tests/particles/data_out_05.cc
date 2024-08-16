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



// Create a particle handler then generate an output
// Tests that output can be written even if no particles are on the local
// process. See https://github.com/dealii/dealii/issues/10443.

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
test()
{
  {
    parallel::distributed::Triangulation<dim, spacedim> tr(MPI_COMM_WORLD);

    GridGenerator::hyper_cube(tr);
    tr.refine_global(2);
    MappingQ<dim, spacedim> mapping(1);

    // Create a particle handler using one manually created particle
    Particles::ParticleHandler<dim, spacedim> particle_handler(tr,
                                                               mapping,
                                                               spacedim);
    std::vector<Point<spacedim>>              position(1);
    std::vector<Point<dim>>                   reference_position(1);
    std::vector<double>                       properties(spacedim);

    for (unsigned int i = 0; i < dim; ++i)
      {
        position[0][i] = 0.125;
      }

    for (unsigned int i = 0; i < spacedim; ++i)
      {
        properties[i] = 0.125;
      }

    Particles::Particle<dim, spacedim> particle1(position[0],
                                                 reference_position[0],
                                                 0);

    typename Triangulation<dim, spacedim>::active_cell_iterator cell1(&tr,
                                                                      2,
                                                                      0);

    if (Utilities::MPI::this_mpi_process(tr.get_mpi_communicator()) == 0)
      {
        auto particle_it = particle_handler.insert_particle(particle1, cell1);
        particle_it->set_properties(properties);
      }

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
