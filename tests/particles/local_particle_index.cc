// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2021 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// Check ParticleAccessor::local_particle_index

#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/mapping_q.h>

#include <deal.II/grid/filtered_iterator.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>

#include <deal.II/particles/generators.h>
#include <deal.II/particles/particle_handler.h>

#include <fstream>
#include <iomanip>

#include "../tests.h"



template <int dim>
void
test()
{
  // create similar mesh as step-68
  parallel::distributed::Triangulation<dim> background_triangulation(
    MPI_COMM_WORLD);

  GridGenerator::hyper_cube(background_triangulation, 0, 1);
  background_triangulation.refine_global(6 - dim);

  const MappingQ<dim> mapping(1);

  Particles::ParticleHandler<dim> particle_handler;
  particle_handler.initialize(background_triangulation, mapping, 1);

  Point<dim> center;
  for (unsigned int d = 0; d < dim; ++d)
    center[d] = 0.5;

  const double outer_radius = 0.15;
  const double inner_radius = 0.01;

  parallel::distributed::Triangulation<dim> particle_triangulation(
    MPI_COMM_WORLD);

  GridGenerator::hyper_shell(
    particle_triangulation, center, inner_radius, outer_radius, 6);
  particle_triangulation.refine_global(5 - dim);

  const auto my_bounding_box = GridTools::compute_mesh_predicate_bounding_box(
    background_triangulation, IteratorFilters::LocallyOwnedCell());
  const auto global_bounding_boxes =
    Utilities::MPI::all_gather(MPI_COMM_WORLD, my_bounding_box);

  std::vector<std::vector<double>> properties(
    particle_triangulation.n_locally_owned_active_cells(),
    std::vector<double>(1, 0.));
  Particles::Generators::quadrature_points(particle_triangulation,
                                           QMidpoint<dim>(),
                                           global_bounding_boxes,
                                           particle_handler,
                                           mapping,
                                           properties);

  for (unsigned int i = 0; i < 2; ++i)
    {
      // Write out some information about the local particles
      deallog << "Number of global particles: "
              << particle_handler.n_global_particles() << std::endl;
      deallog << "Number of locally owned particles: "
              << particle_handler.n_locally_owned_particles() << std::endl;
      deallog << "Maximal local particle index: "
              << particle_handler.get_max_local_particle_index() << std::endl;

      std::vector<types::particle_index> my_particle_numbers(
        particle_handler.get_max_local_particle_index(),
        numbers::invalid_dof_index);

      // set the id numbers for the particles
      unsigned int count = 0;
      for (auto particle = particle_handler.begin();
           particle != particle_handler.end();
           ++particle, ++count)
        my_particle_numbers[particle->get_local_index()] = particle->get_id();

      for (auto particle = particle_handler.begin_ghost();
           particle != particle_handler.end_ghost();
           ++particle, ++count)
        my_particle_numbers[particle->get_local_index()] = particle->get_id();

      deallog << "Set information of " << count << " particles" << std::endl;

      unsigned int n_unset = 0;
      for (const types::particle_index id : my_particle_numbers)
        if (id == numbers::invalid_dof_index)
          ++n_unset;

      deallog << "Number of particles not set: " << n_unset << std::endl;

      particle_handler.exchange_ghost_particles();
    }
  deallog << "OK" << std::endl;
}



int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi(argc, argv, 1);
  MPILogInitAll                    log;

  test<2>();
  test<3>();
}
