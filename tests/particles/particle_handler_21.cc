/* ------------------------------------------------------------------------
 *
 * SPDX-License-Identifier: LGPL-2.1-or-later
 * Copyright (C) 2022 - 2024 by the deal.II authors
 *
 * This file is part of the deal.II library.
 *
 * Part of the source code is dual licensed under Apache-2.0 WITH
 * LLVM-exception OR LGPL-2.1-or-later. Detailed license information
 * governing the source code and code contributions can be found in
 * LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
 *
 * ------------------------------------------------------------------------

 * This test case is an extremely simplified version of step-68.
 * A ball made of 48 particles is generated. This ball is moved down slightly.
 * The particles remain in the simulation domain and they should not disappear.
 * At the moment of the creation of this test, a bug in the particle_handler
 would
 * make one particle disappear and the number would decrease from 48 to 47 at
 the
 * second time step.
 */

// Include files

#include <deal.II/base/bounding_box.h>
#include <deal.II/base/conditional_ostream.h>
#include <deal.II/base/discrete_time.h>
#include <deal.II/base/mpi.h>
#include <deal.II/base/parameter_acceptor.h>
#include <deal.II/base/timer.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/mapping_q1.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_tools.h>

#include <deal.II/lac/la_parallel_vector.h>
#include <deal.II/lac/vector.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/vector_tools.h>

#include <deal.II/particles/data_out.h>
#include <deal.II/particles/generators.h>
#include <deal.II/particles/particle_handler.h>

#include "../tests.h"


template <int dim>
class VelocityField : public Function<dim>
{
public:
  VelocityField()
    : Function<dim>(dim)
  {}

  virtual void
  vector_value(const Point<dim> &point, Vector<double> &values) const override;
};

template <int dim>
void
VelocityField<dim>::vector_value(const Point<dim> & /*point*/,
                                 Vector<double> &values) const
{
  values[0] = 0;
  values[1] = -1;
  if (dim > 2)
    values[2] = 0;
}

template <int dim>
class ParticleTracking
{
public:
  ParticleTracking();
  void
  run();

private:
  void
  generate_particles();

  void
  euler_step_analytical(const double dt);

  MPI_Comm                                  mpi_communicator;
  parallel::distributed::Triangulation<dim> background_triangulation;
  Particles::ParticleHandler<dim>           particle_handler;

  MappingQ1<dim>     mapping;
  VelocityField<dim> velocity;
};

template <int dim>
ParticleTracking<dim>::ParticleTracking()
  : mpi_communicator(MPI_COMM_WORLD)
  , background_triangulation(mpi_communicator)
{}

// This function generates the tracer particles and the background
// triangulation on which these particles evolve.

template <int dim>
void
ParticleTracking<dim>::generate_particles()
{
  // We create an hyper ball triangulation which we globally refine. The bug
  // that this test tries to reproduce only occurred in curved geometry.
  Point<dim> center_of_triangulation;
  center_of_triangulation[0] = 0.;
  center_of_triangulation[1] = 0.;
  if (dim == 3)
    center_of_triangulation[2] = 0.;

  GridGenerator::hyper_ball(background_triangulation,
                            center_of_triangulation,
                            1);
  background_triangulation.refine_global(3);

  // This initializes the background triangulation where the particles are
  // living and the number of properties of the particles.
  particle_handler.initialize(background_triangulation, mapping, 1 + dim);

  // We create a particle triangulation which is solely used to generate
  // the points which will be used to insert the particles. This
  // triangulation is a hyper shell which is offset from the
  // center of the simulation domain. We generate enough particle to reproduce
  // the bug.

  Point<dim> center_of_particles;
  center_of_particles[0] = 0.0;
  center_of_particles[1] = -.735;
  if (dim == 3)
    center_of_particles[2] = 0.0;

  const double outer_radius = 0.2;
  const double inner_radius = 0.01;

  parallel::distributed::Triangulation<dim> particle_triangulation(
    mpi_communicator);

  GridGenerator::hyper_shell(
    particle_triangulation, center_of_particles, inner_radius, outer_radius, 6);
  particle_triangulation.refine_global(1);

  // We generate the necessary bounding boxes for the particles generator.
  // These bounding boxes are required to quickly identify in which
  // process's subdomain the inserted particle lies, and which cell owns it.
  const auto my_bounding_box = GridTools::compute_mesh_predicate_bounding_box(
    background_triangulation, IteratorFilters::LocallyOwnedCell());
  const auto global_bounding_boxes =
    Utilities::MPI::all_gather(mpi_communicator, my_bounding_box);

  // We generate an empty vector of properties. We will attribute the
  // properties to the particles once they are generated.
  std::vector<std::vector<double>> properties(
    particle_triangulation.n_locally_owned_active_cells(),
    std::vector<double>(dim + 1, 0.));

  // We generate the particles at the position of a single
  // point quadrature. Consequently, one particle will be generated
  // at the centroid of each cell.
  Particles::Generators::quadrature_points(particle_triangulation,
                                           QMidpoint<dim>(),
                                           global_bounding_boxes,
                                           particle_handler,
                                           mapping,
                                           properties);

  deallog << "Number of particles inserted: "
          << particle_handler.n_global_particles() << std::endl;
}

// We integrate the particle trajectories using a first order explicit Euler
// scheme.
template <int dim>
void
ParticleTracking<dim>::euler_step_analytical(const double dt)
{
  const unsigned int this_mpi_rank =
    Utilities::MPI::this_mpi_process(mpi_communicator);
  Vector<double> particle_velocity(dim);

  // Looping over all particles in the domain using a
  // particle iterator
  for (auto &particle : particle_handler)
    {
      // We calculate the velocity of the particles using their current
      // location.
      Point<dim> particle_location = particle.get_location();
      velocity.vector_value(particle_location, particle_velocity);

      // This updates the position of the particles and sets the old position
      // equal to the new position of the particle.
      for (int d = 0; d < dim; ++d)
        particle_location[d] += particle_velocity[d] * dt;

      particle.set_location(particle_location);

      // We store the processor id (a scalar) and the particle velocity (a
      // vector) in the particle properties.
      ArrayView<double> properties = particle.get_properties();
      for (int d = 0; d < dim; ++d)
        properties[d] = particle_velocity[d];
      properties[dim] = this_mpi_rank;
    }
}

template <int dim>
void
ParticleTracking<dim>::run()
{
  DiscreteTime discrete_time(0, 0.015, 0.005);

  generate_particles();

  // The particles are advected by looping over time.
  while (!discrete_time.is_at_end())
    {
      discrete_time.advance_time();
      velocity.set_time(discrete_time.get_previous_time());

      euler_step_analytical(discrete_time.get_previous_step_size());

      unsigned int n_part_before_sort = particle_handler.n_global_particles();

      particle_handler.sort_particles_into_subdomains_and_cells();
      unsigned int n_part_after_sort = particle_handler.n_global_particles();

      deallog << "Number of particles before sort : " << n_part_before_sort
              << " After : " << n_part_after_sort << std::endl;
    }
}

int
main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

  MPILogInitAll all;

  ParticleTracking<3> particle_advection_problem;
  particle_advection_problem.run();
}
