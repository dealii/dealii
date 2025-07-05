/* ------------------------------------------------------------------------
 *
 * SPDX-License-Identifier: LGPL-2.1-or-later
 * Copyright (C) 2023 - 2024 by the deal.II authors
 *
 * This file is part of the deal.II library.
 *
 * Part of the source code is dual licensed under Apache-2.0 WITH
 * LLVM-exception OR LGPL-2.1-or-later. Detailed license information
 * governing the source code and code contributions can be found in
 * LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
 *
 * ------------------------------------------------------------------------

 * This test is based on a modified version of step-68. Modifications were
 * made to ensure that particles were more evenly distributed
 * between the cells to ensure that the computational load could be
 * balanced in some way. This is a more realistic usage of particles.
 * The performance test measures
 * the interpolation of a finite element field to particle location,
 * the displacement of particles and their localizaton within subdomains
 * and cells.
 */


#include <deal.II/base/bounding_box.h>
#include <deal.II/base/conditional_ostream.h>
#include <deal.II/base/discrete_time.h>
#include <deal.II/base/function_lib.h>
#include <deal.II/base/mpi.h>
#include <deal.II/base/parameter_acceptor.h>
#include <deal.II/base/timer.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/mapping_q1.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>

#include <deal.II/lac/la_parallel_vector.h>
#include <deal.II/lac/vector.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/vector_tools.h>

#include <deal.II/particles/data_out.h>
#include <deal.II/particles/generators.h>
#include <deal.II/particles/particle_handler.h>

#include <cmath>
#include <iostream>

#define ENABLE_MPI

#include "performance_test_driver.h"

constexpr bool debug = false;

namespace Step68
{
  using namespace dealii;


  struct ParticleTrackingParameters
  {
    ParticleTrackingParameters()
    {
      if (get_testing_environment() == TestingEnvironment::medium)
        {
          fluid_refinement              = 8;
          particle_insertion_refinement = 10;
        }
      else if (get_testing_environment() == TestingEnvironment::heavy)
        {
          fluid_refinement              = 9;
          particle_insertion_refinement = 11;
        }
    }

    unsigned int velocity_degree       = 1;
    double       time_step             = 0.0001;
    double       final_time            = 0.0050;
    unsigned int output_frequency      = 40;
    unsigned int repartition_frequency = 5;

    unsigned int fluid_refinement              = 7;
    unsigned int particle_insertion_refinement = 9;
  };



  template <int dim>
  class ParticleTracking
  {
  public:
    ParticleTracking(const ParticleTrackingParameters &par,
                     const bool                        interpolated_velocity);

    Measurement
    run();

  private:
    void
    generate_particles();

    void
    setup_background_dofs();

    void
    interpolate_function_to_field();

    void
    euler_step_interpolated(const double dt);

    unsigned int
    cell_weight(
      const typename parallel::distributed::Triangulation<dim>::cell_iterator
        &cell,
      const typename parallel::distributed::Triangulation<dim>::CellStatus
        status) const;

    void
    output_particles();
    void
    output_background();


    const ParticleTrackingParameters &par;

    MPI_Comm                                  mpi_communicator;
    parallel::distributed::Triangulation<dim> background_triangulation;
    Particles::ParticleHandler<dim>           particle_handler;

    DoFHandler<dim>                            fluid_dh;
    FESystem<dim>                              fluid_fe;
    MappingQ1<dim>                             mapping;
    LinearAlgebra::distributed::Vector<double> velocity_field;

    Functions::RayleighKotheVortex<dim> velocity;

    ConditionalOStream pcout;

    bool interpolated_velocity;
  };



  template <int dim>
  ParticleTracking<dim>::ParticleTracking(const ParticleTrackingParameters &par,
                                          const bool interpolated_velocity)
    : par(par)
    , mpi_communicator(MPI_COMM_WORLD)
    , background_triangulation(mpi_communicator)
    , fluid_dh(background_triangulation)
    , fluid_fe(FE_Q<dim>(par.velocity_degree), dim)
    , velocity(4.0)
    , pcout(std::cout,
            debug && Utilities::MPI::this_mpi_process(mpi_communicator) == 0)
    , interpolated_velocity(interpolated_velocity)

  {}



  template <int dim>
  unsigned int
  ParticleTracking<dim>::cell_weight(
    const typename parallel::distributed::Triangulation<dim>::cell_iterator
                                                                        &cell,
    const typename parallel::distributed::Triangulation<dim>::CellStatus status)
    const
  {
    const unsigned int base_weight = 1;

    const unsigned int particle_weight = 1;

    unsigned int n_particles_in_cell = 0;
    switch (status)
      {
        case parallel::distributed::Triangulation<dim>::CELL_PERSIST:
        case parallel::distributed::Triangulation<dim>::CELL_REFINE:
          n_particles_in_cell = particle_handler.n_particles_in_cell(cell);
          break;

        case parallel::distributed::Triangulation<dim>::CELL_INVALID:
          break;

        case parallel::distributed::Triangulation<dim>::CELL_COARSEN:
          for (const auto &child : cell->child_iterators())
            n_particles_in_cell += particle_handler.n_particles_in_cell(child);
          break;

        default:
          Assert(false, ExcInternalError());
          break;
      }

    return base_weight + particle_weight * n_particles_in_cell;
  }



  template <int dim>
  void
  ParticleTracking<dim>::generate_particles()
  {
    GridGenerator::hyper_cube(background_triangulation, 0, 1);
    background_triangulation.refine_global(par.fluid_refinement);

    background_triangulation.signals.weight.connect(
      [&](
        const typename parallel::distributed::Triangulation<dim>::cell_iterator
          &cell,
        const typename parallel::distributed::Triangulation<dim>::CellStatus
          status) -> unsigned int { return this->cell_weight(cell, status); });

    particle_handler.initialize(background_triangulation, mapping, 1 + dim);

    Point<dim> center;
    center[0] = 0.0;
    center[1] = 0.0;
    if (dim == 3)
      center[2] = 0.0;

    parallel::distributed::Triangulation<dim> particle_triangulation(
      MPI_COMM_WORLD);

    GridGenerator::hyper_cube(particle_triangulation, 0, 1);
    particle_triangulation.refine_global(par.particle_insertion_refinement);

    const auto my_bounding_box = GridTools::compute_mesh_predicate_bounding_box(
      background_triangulation, IteratorFilters::LocallyOwnedCell());
    const auto global_bounding_boxes =
      Utilities::MPI::all_gather(MPI_COMM_WORLD, my_bounding_box);

    std::vector<std::vector<double>> properties(
      particle_triangulation.n_locally_owned_active_cells(),
      std::vector<double>(dim + 1, 0.));

    Particles::Generators::quadrature_points(particle_triangulation,
                                             QMidpoint<dim>(),
                                             global_bounding_boxes,
                                             particle_handler,
                                             mapping,
                                             properties);

    pcout << "Number of particles inserted: "
          << particle_handler.n_global_particles() << std::endl;
  }



  template <int dim>
  void
  ParticleTracking<dim>::setup_background_dofs()
  {
    fluid_dh.distribute_dofs(fluid_fe);
    const IndexSet &locally_owned_dofs = fluid_dh.locally_owned_dofs();
    const IndexSet  locally_relevant_dofs =
      DoFTools::extract_locally_relevant_dofs(fluid_dh);

    velocity_field.reinit(locally_owned_dofs,
                          locally_relevant_dofs,
                          mpi_communicator);
  }



  template <int dim>
  void
  ParticleTracking<dim>::interpolate_function_to_field()
  {
    velocity_field.zero_out_ghost_values();
    VectorTools::interpolate(mapping, fluid_dh, velocity, velocity_field);
    velocity_field.update_ghost_values();
  }

  template <int dim>
  void
  ParticleTracking<dim>::euler_step_interpolated(const double dt)
  {
    Vector<double> local_dof_values(fluid_fe.dofs_per_cell);

    FEPointEvaluation<dim, dim> evaluator(mapping, fluid_fe, update_values);
    std::vector<Point<dim>>     particle_positions;

    const double this_mpi_process =
      static_cast<double>(Utilities::MPI::this_mpi_process(mpi_communicator));

    auto particle = particle_handler.begin();
    while (particle != particle_handler.end())
      {
        const auto cell = particle->get_surrounding_cell();
        const auto dh_cell =
          typename DoFHandler<dim>::cell_iterator(*cell, &fluid_dh);

        dh_cell->get_dof_values(velocity_field, local_dof_values);

        const auto pic = particle_handler.particles_in_cell(cell);
        Assert(pic.begin() == particle, ExcInternalError());
        particle_positions.clear();
        for (auto &p : pic)
          particle_positions.push_back(p.get_reference_location());

        evaluator.reinit(cell, particle_positions);
        evaluator.evaluate(make_array_view(local_dof_values),
                           EvaluationFlags::values);

        for (unsigned int particle_index = 0; particle != pic.end();
             ++particle, ++particle_index)
          {
            Point<dim>           &particle_location = particle->get_location();
            const Tensor<1, dim> &particle_velocity =
              evaluator.get_value(particle_index);
            particle_location += particle_velocity * dt;

            ArrayView<double> properties = particle->get_properties();
            for (int d = 0; d < dim; ++d)
              properties[d] = particle_velocity[d];

            properties[dim] = this_mpi_process;
          }
      }
  }



  template <int dim>
  void
  ParticleTracking<dim>::output_particles()
  {
    Particles::DataOut<dim, dim> particle_output;

    std::vector<std::string> solution_names(dim, "velocity");
    solution_names.emplace_back("process_id");

    std::vector<DataComponentInterpretation::DataComponentInterpretation>
      data_component_interpretation(
        dim, DataComponentInterpretation::component_is_part_of_vector);
    data_component_interpretation.push_back(
      DataComponentInterpretation::component_is_scalar);

    particle_output.build_patches(particle_handler,
                                  solution_names,
                                  data_component_interpretation);
  }



  template <int dim>
  void
  ParticleTracking<dim>::output_background()
  {
    std::vector<std::string> solution_names(dim, "velocity");
    std::vector<DataComponentInterpretation::DataComponentInterpretation>
      data_component_interpretation(
        dim, DataComponentInterpretation::component_is_part_of_vector);

    DataOut<dim> data_out;

    data_out.attach_dof_handler(fluid_dh);
    data_out.add_data_vector(velocity_field,
                             solution_names,
                             DataOut<dim>::type_dof_data,
                             data_component_interpretation);
    Vector<float> subdomain(background_triangulation.n_active_cells());
    for (unsigned int i = 0; i < subdomain.size(); ++i)
      subdomain(i) = background_triangulation.locally_owned_subdomain();
    data_out.add_data_vector(subdomain, "subdomain");

    data_out.build_patches(mapping);
  }



  template <int dim>
  Measurement
  ParticleTracking<dim>::run()
  {
    std::map<std::string, dealii::Timer> timer;

    DiscreteTime discrete_time(0, par.final_time, par.time_step);

    timer["generate_particles"].start();
    generate_particles();
    timer["generate_particles"].stop();

    pcout << "Repartitioning triangulation after particle generation"
          << std::endl;

    particle_handler.prepare_for_coarsening_and_refinement();
    background_triangulation.repartition();
    particle_handler.unpack_after_coarsening_and_refinement();

    setup_background_dofs();
    interpolate_function_to_field();
    euler_step_interpolated(0.);


    while (!discrete_time.is_at_end())
      {
        discrete_time.advance_time();
        velocity.set_time(discrete_time.get_previous_time());

        if ((discrete_time.get_step_number() % par.repartition_frequency) == 0)
          {
            timer["load_balance"].start();
            particle_handler.prepare_for_coarsening_and_refinement();
            background_triangulation.repartition();
            particle_handler.unpack_after_coarsening_and_refinement();

            setup_background_dofs();
            timer["load_balance"].stop();
          }


        timer["advect"].start();
        interpolate_function_to_field();
        euler_step_interpolated(discrete_time.get_previous_step_size());
        timer["advect"].stop();


        timer["sort"].start();
        particle_handler.sort_particles_into_subdomains_and_cells();
        timer["sort"].stop();

        if ((discrete_time.get_step_number() % par.output_frequency) == 0)
          {
            timer["output"].start();
            output_particles();
            output_background();
            timer["output"].stop();
          }
      }

    return {timer["generate_particles"].wall_time(),
            timer["load_balance"].wall_time(),
            timer["advect"].wall_time(),
            timer["sort"].wall_time(),
            timer["output"].wall_time()};
  }

} // namespace Step68


std::tuple<Metric, unsigned int, std::vector<std::string>>
describe_measurements()
{
  return {Metric::timing,
          4,
          {"generate_particles", "load_balance", "advect", "sort", "output"}};
}


Measurement
perform_single_measurement()
{
  using namespace Step68;
  using namespace dealii;

  ParticleTrackingParameters par;

  Step68::ParticleTracking<2> particle_tracking(par, true);
  return particle_tracking.run();
}
