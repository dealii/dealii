// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2020 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

// Test that an interpolation sparsity can be constructed for a single
// particle per processor for every valid dimension pair that exist
// when there are more than one components on the DoFHandler
// that is associated with the triangulation, and we select only one
// component

#include <deal.II/base/point.h>
#include <deal.II/base/tensor.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/mapping_q1.h>

#include <deal.II/grid/grid_generator.h>

#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/sparsity_pattern.h>
#include <deal.II/lac/sparsity_tools.h>

#include <deal.II/particles/generators.h>
#include <deal.II/particles/particle_handler.h>
#include <deal.II/particles/utilities.h>

#include "../tests.h"


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


  FESystem<dim, spacedim> space_fe(FE_Q<dim, spacedim>(1), 3);

  ComponentMask space_mask({false, false, true});

  const auto n_comps = space_mask.n_selected_components();

  DoFHandler<dim, spacedim> space_dh(space_tria);
  space_dh.distribute_dofs(space_fe);

  const IndexSet &locally_owned_dofs = space_dh.locally_owned_dofs();
  const IndexSet  locally_relevant_dofs =
    DoFTools::extract_locally_relevant_dofs(space_dh);

  if (my_mpi_id == 0)
    deallog << "dim: " << dim << ", spacedim: " << spacedim << std::endl
            << "Total number of particles: "
            << particle_handler.n_global_particles() << std::endl

            << "Space FE: " << space_fe.get_name() << std::endl
            << "Space dofs: " << space_dh.n_dofs() << std::endl;

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

  DynamicSparsityPattern dsp(particle_handler.n_global_particles() * n_comps,
                             space_dh.n_dofs(),
                             local_particle_index_set);

  AffineConstraints<double> constraints;

  // Build the interpolation sparsity
  Particles::Utilities::create_interpolation_sparsity_pattern(
    space_dh, particle_handler, dsp, constraints, space_mask);

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
