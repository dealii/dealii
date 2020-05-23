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



// check the particle generation using the
// Particles::Generator::probabilistic_locations function
// using different probability distributions. In particular
// check that for a non-random cell selection each cell
// has the expected number of particles according to its
// part of the probability density distribution.

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
    tr.refine_global(1);

    MappingQ<dim, spacedim> mapping(1);

    Particles::ParticleHandler<dim, spacedim> particle_handler(tr, mapping);

    Particles::Generators::probabilistic_locations(
      tr, probability_density_function, false, 16, particle_handler, mapping);

    if (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0)
      {
        deallog << "Locally owned active cells: "
                << tr.n_locally_owned_active_cells() << std::endl;

        deallog << "Global particles: " << particle_handler.n_global_particles()
                << std::endl;

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
                    << particle.get_surrounding_cell(tr) << std::endl;
            deallog << "Particle location: " << particle.get_location()
                    << std::endl;
          }
      }
  }

  deallog << "OK" << std::endl;
}



int
main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

  initlog();

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

  {
    deallog.push("2d/2d");
    Functions::Q1WedgeFunction<2> probability_density_function;
    test<2, 2>(probability_density_function);
    deallog.pop();
  }
  {
    deallog.push("2d/3d");
    Functions::Q1WedgeFunction<3> probability_density_function;
    test<2, 3>(probability_density_function);
  }
  {
    deallog.pop();
    deallog.push("3d/3d");
    Functions::Q1WedgeFunction<3> probability_density_function;
    test<3, 3>(probability_density_function);
    deallog.pop();
  }
}
