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

// Test that a Monomial function interpolated to a field
// can be interpolated at the particle position correctly

#include <deal.II/base/function_lib.h>
#include <deal.II/base/point.h>
#include <deal.II/base/tensor.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>

#include <deal.II/grid/grid_generator.h>

#include <deal.II/lac/trilinos_vector.h>

#include <deal.II/numerics/vector_tools.h>

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

  // Create a triangulation with a non-trivial number of degree of freedom
  parallel::distributed::Triangulation<dim, spacedim> tria(MPI_COMM_WORLD);
  MappingQ1<dim, spacedim>                            mapping;
  GridGenerator::hyper_cube(tria, -1, 1);
  tria.refine_global(2);

  // Generate particles at the gauss points where the field will be interpolated
  Particles::ParticleHandler<dim, spacedim> particle_handler(tria, mapping);

  Particles::Generators::regular_reference_locations<dim, spacedim>(
    tria, QGauss<dim>(2).get_points(), particle_handler);


  FE_Q<dim, spacedim> fe(1);

  ComponentMask mask(fe.n_components(), true);
  const auto    n_comps = mask.n_selected_components();

  DoFHandler<dim, spacedim> space_dh(tria);
  space_dh.distribute_dofs(fe);

  IndexSet locally_owned_dofs = space_dh.locally_owned_dofs();
  IndexSet locally_relevant_dofs;
  DoFTools::extract_locally_relevant_dofs(space_dh, locally_relevant_dofs);

  if (my_mpi_id == 0)
    deallog << "dim: " << dim << ", spacedim: " << spacedim << std::endl
            << "Total number of particles: "
            << particle_handler.n_global_particles() << std::endl;

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

  // Create dofs and particle vectors
  TrilinosWrappers::MPI::Vector field_owned(locally_owned_dofs, MPI_COMM_WORLD);
  TrilinosWrappers::MPI::Vector field_relevant(locally_relevant_dofs,
                                               MPI_COMM_WORLD);
  TrilinosWrappers::MPI::Vector interpolation_on_particles(
    local_particle_index_set, MPI_COMM_WORLD);

  // Create interpolation function
  Tensor<1, spacedim> exponents;
  for (unsigned int i = 0; i < spacedim; ++i)
    exponents[i] = 1;
  Functions::Monomial<spacedim> linear(exponents, mask.size());

  // Interpolate function to vector than interpolate field to particles
  VectorTools::interpolate(space_dh, linear, field_owned);
  field_relevant = field_owned;
  Particles::Utilities::interpolate_field_on_particles(
    space_dh, particle_handler, field_relevant, interpolation_on_particles);

  Vector<double> values(mask.size());

  for (auto particle : particle_handler)
    {
      const auto &location = particle.get_location();
      const auto &id       = particle.get_id();
      linear.vector_value(location, values);
      for (unsigned int i = 0, j = 0; i < values.size(); ++i)
        if (mask[i])
          if (std::abs(values[i] - interpolation_on_particles(id * n_comps +
                                                              (j++))) > 1e-10)
            deallog << "NOT OK" << std::endl;
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
