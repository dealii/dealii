// ---------------------------------------------------------------------
//
// Copyright (C) 2019 - 2020 by the deal.II authors
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

// Test that an interpolation sparsity can be constructed for a single
// particle per processor for every valid dimension pair that exist

#include <deal.II/base/point.h>
#include <deal.II/base/tensor.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>

#include <deal.II/grid/grid_generator.h>

#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/sparsity_pattern.h>
#include <deal.II/lac/sparsity_tools.h>

#include <deal.II/particles/generators.h>
#include <deal.II/particles/particle_handler.h>
#include <deal.II/particles/utilities.h>

#include "../tests.h"

using namespace dealii;

template <int dim, int spacedim>
void
test()
{
  const unsigned int my_mpi_id =
    Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);

  parallel::distributed::Triangulation<dim, spacedim> space_tria(
    MPI_COMM_WORLD);
  MappingQ1<dim, spacedim> mapping;
  GridGenerator::hyper_cube(space_tria, -1, 1);
  space_tria.refine_global(2);

  // Start by creating a particle handler that contains one particle per
  // processor
  Particles::ParticleHandler<dim, spacedim> particle_handler(space_tria,
                                                             mapping);

  const unsigned int estimated_n_particles_per_processor = 1;

  Particles::Generators::probabilistic_locations<dim, spacedim>(
    space_tria,
    Functions::ConstantFunction<spacedim, double>(1.),
    true,
    Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD) *
      estimated_n_particles_per_processor,
    particle_handler);


  FE_Q<dim, spacedim> space_fe(1);

  ComponentMask space_mask(space_fe.n_components(), true);

  const auto n_comps = space_mask.n_selected_components();

  DoFHandler<dim, spacedim> space_dh(space_tria);
  space_dh.distribute_dofs(space_fe);

  IndexSet locally_owned_dofs = space_dh.locally_owned_dofs();
  IndexSet locally_relevant_dofs;
  DoFTools::extract_locally_relevant_dofs(space_dh, locally_relevant_dofs);

  if (my_mpi_id == 0)
    deallog << "dim: " << dim << ", spacedim: " << spacedim << std::endl
            << "Total number of particles: "
            << particle_handler.n_global_particles() << std::endl

            << "Space FE: " << space_fe.get_name() << std::endl
            << "Space dofs: " << space_dh.n_dofs() << std::endl;

  DynamicSparsityPattern dsp(particle_handler.n_global_particles() * n_comps,
                             space_dh.n_dofs());

  AffineConstraints<double> constraints;

  // Build the interpolation sparsity
  Particles::Utilities::create_interpolation_sparsity_pattern(
    space_dh, particle_handler, dsp, constraints, space_mask);

  const auto n_local_particles_dofs =
    particle_handler.n_locally_owned_particles() * n_comps;

  auto particle_sizes =
    Utilities::MPI::all_gather(MPI_COMM_WORLD, n_local_particles_dofs);

  const auto my_start = std::accumulate(particle_sizes.begin(),
                                        particle_sizes.begin() + my_mpi_id,
                                        0u);

  IndexSet local_particle_index_set(particle_handler.n_global_particles() *
                                    n_comps);

  local_particle_index_set.add_range(my_start,
                                     my_start + n_local_particles_dofs);

  auto global_particles_index_set =
    Utilities::MPI::all_gather(MPI_COMM_WORLD, n_local_particles_dofs);

  SparsityTools::distribute_sparsity_pattern(dsp,
                                             global_particles_index_set,
                                             MPI_COMM_WORLD,
                                             local_particle_index_set);
  // Actually log the sparsity
  {
    SparsityPattern sparsity_pattern;
    sparsity_pattern.copy_from(dsp);
    sparsity_pattern.print_gnuplot(deallog.get_file_stream());
  }
}

int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

  MPILogInitAll init;

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
