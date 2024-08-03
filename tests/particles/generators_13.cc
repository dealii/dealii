// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2019 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// like generators_05, but check that empty triangulations do not crash
// (by using more cores than cells).

#include <deal.II/base/function_lib.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_nothing.h>
#include <deal.II/fe/mapping_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>

#include <deal.II/particles/generators.h>
#include <deal.II/particles/particle_handler.h>

#include "../tests.h"

template <int dim, int spacedim>
void
test(Function<spacedim> &probability_density_function)
{
  {
    parallel::distributed::Triangulation<dim, spacedim> tr(MPI_COMM_WORLD);

    GridGenerator::hyper_cube(tr);
    MappingQ<dim, spacedim> mapping(1);

    Particles::ParticleHandler<dim, spacedim> particle_handler(tr, mapping);

    Particles::Generators::probabilistic_locations(
      tr, probability_density_function, true, 16, particle_handler, mapping);

    deallog << "Rank: " << Utilities::MPI::this_mpi_process(MPI_COMM_WORLD)
            << ". Locally owned active cells: "
            << tr.n_locally_owned_active_cells() << ". Local particles: "
            << particle_handler.n_locally_owned_particles() << std::endl;

    for (const auto &cell : tr.active_cell_iterators())
      {
        if (cell->is_locally_owned())
          {
            deallog << "Cell " << cell << " has "
                    << particle_handler.n_particles_in_cell(cell)
                    << " particles." << std::endl;
          }
      }

    for (const auto &particle : particle_handler)
      {
        deallog << "Particle index " << particle.get_id() << " is in cell "
                << particle.get_surrounding_cell() << std::endl;
        deallog << "Particle location: " << particle.get_location()
                << std::endl;
      }
  }

  deallog << "OK" << std::endl;
}



int
main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

  MPILogInitAll init;

  {
    deallog.push("2d/2d");
    Functions::ConstantFunction<2> probability_density_function(1.0);
    test<2, 2>(probability_density_function);
    deallog.pop();
  }
  {
    deallog.push("2d/3d");
    Functions::ConstantFunction<3> probability_density_function(1.0);
    test<2, 3>(probability_density_function);
  }
  {
    deallog.pop();
    deallog.push("3d/3d");
    Functions::ConstantFunction<3> probability_density_function(1.0);
    test<3, 3>(probability_density_function);
    deallog.pop();
  }
}
