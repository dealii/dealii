/* ---------------------------------------------------------------------
 *
 * Copyright (C) 2020 by the deal.II authors
 *
 * This file is part of the deal.II library.
 *
 * The deal.II library is free software; you can use it, redistribute
 * it, and/or modify it under the terms of the GNU Lesser General
 * Public License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 * The full text of the license can be found in the file LICENSE.md at
 * the top level directory of deal.II.
 *
 * ---------------------------------------------------------------------

 *
 * Authors: Bruno Blais, Toni El Geitani Nehme, Rene Gassmoeller, Luca Heltai,
 Wolfgang Banghert 2020
 */

#include <deal.II/base/bounding_box.h>
#include <deal.II/base/conditional_ostream.h>
#include <deal.II/base/discrete_time.h>
#include <deal.II/base/index_set.h>
#include <deal.II/base/mpi.h>
#include <deal.II/base/parameter_acceptor.h>
#include <deal.II/base/std_cxx14/memory.h>
#include <deal.II/base/timer.h>

#include <deal.II/distributed/cell_weights.h>
#include <deal.II/distributed/grid_refinement.h>
#include <deal.II/distributed/solution_transfer.h>
#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/mapping_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/generic_linear_algebra.h>
#include <deal.II/lac/petsc_vector.h>
#include <deal.II/lac/sparsity_tools.h>
#include <deal.II/lac/trilinos_vector.h>
#include <deal.II/lac/vector.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/vector_tools.h>

#include <deal.II/particles/data_out.h>
#include <deal.II/particles/generators.h>
#include <deal.II/particles/particle_handler.h>
#include <deal.II/particles/utilities.h>

#define FORCE_USE_OF_TRILINOS

namespace LA
{
#if defined(DEAL_II_WITH_PETSC) && !defined(DEAL_II_PETSC_WITH_COMPLEX) && \
  !(defined(DEAL_II_WITH_TRILINOS) && defined(FORCE_USE_OF_TRILINOS))
  using namespace dealii::LinearAlgebraPETSc;
#  define USE_PETSC_LA
#elif defined(DEAL_II_WITH_TRILINOS)
  using namespace dealii::LinearAlgebraTrilinos;
#else
#  error DEAL_II_WITH_PETSC or DEAL_II_WITH_TRILINOS required
#endif
} // namespace LA

#include <cmath>
#include <fstream>
#include <iostream>
#include <memory>

namespace Stepx
{
  using namespace dealii;

  class ParticleTrackingParameters : public ParameterAcceptor
  {
  public:
    ParticleTrackingParameters();

    // This class consists largely of member variables that
    // describe the details of the particle tracking simulation and its
    // discretization. The following parameters are about where output should
    // land, the spatial discretization of the velocity (the default is $Q_1$),
    // the time step, finally,  the output frequency (how many time steps should
    // elapse before we generate graphical output again):
    std::string output_directory = "./";

    unsigned int velocity_degree  = 1;
    double       time_step        = 0.002;
    double       final_time       = 4.0;
    unsigned int output_frequency = 10;
    unsigned int repartition_frequency = 5;

    // We allow every grid to be refined independently. In this tutorial, no
    // physics is resolved on the fluid grid, and its velocity is calculated
    // analytically.
    unsigned int fluid_refinement              = 4;
    unsigned int particle_insertion_refinement = 3;
  };



  // There remains the task of declaring what run-time parameters we can accept
  // in input files.
  ParticleTrackingParameters::ParticleTrackingParameters()
    : ParameterAcceptor("Particle Tracking Problem/")
  {
    add_parameter(
      "Velocity degree", velocity_degree, "", this->prm, Patterns::Integer(1));

    add_parameter("Output frequency",
                  output_frequency,
                  "Iteration frequency at which output results are written");

                      add_parameter("Repartition frequency",
                  repartition_frequency,
                  "Iteration frequency at which the mesh is load balanced");

    add_parameter("Output directory", output_directory);

    add_parameter("Time step", time_step);

    add_parameter("Final time", final_time, "End time of the simulation");

    add_parameter("Fluid refinement",
                  fluid_refinement,
                  "Refinement level of the fluid domain");

    add_parameter(
      "Particle insertion refinement",
      particle_insertion_refinement,
      "Refinement of the volumetric mesh used to insert the particles");
  }

  // Creating the function for the velocity profile.
  template <int dim>
  class SingleVortex : public Function<dim>
  {
  public:
    SingleVortex()
      : Function<dim>(dim)
    {}
    virtual void
    vector_value(const Point<dim> &point,
                 Vector<double> &  values) const override;
  };

  template <int dim>
  void
  SingleVortex<dim>::vector_value(const Point<dim> &point,
                                  Vector<double> &  values) const
  {
    const double T = 4;
    const double t = this->get_time();

    const double px = numbers::PI * point(0);
    const double py = numbers::PI * point(1);
    const double pt = numbers::PI / T * t;

    values[0] = -2 * cos(pt) * pow(sin(px), 2) * sin(py) * cos(py);
    values[1] = 2 * cos(pt) * pow(sin(py), 2) * sin(px) * cos(px);
    if (dim == 3)
      {
        values[2] = 0;
      }
  }

  // Solver
  template <int dim>
  class ParticleTracking
  {
  public:
    ParticleTracking(const ParticleTrackingParameters &par,
                     const bool                        interpolated_velocity);
    void
    run_analytical_velocity();

    unsigned int
    cell_weight(const typename parallel::distributed::Triangulation<dim>::cell_iterator &cell,
                const typename parallel::distributed::Triangulation<dim>::CellStatus status);

  private:
    void
    particles_generation();
    void
    setup_background_dofs();
    void
    interpolate_function_to_field();
    void
    euler_interpolated(double dt);
    void
    euler_analytical(double dt);
    void
    field_euler(double t, double dt, double T);
    void
    output_particles(unsigned int it);
    void
    output_background(unsigned int it);

    const ParticleTrackingParameters &par;


    MPI_Comm                                  mpi_communicator;
    parallel::distributed::Triangulation<dim> background_triangulation;
    Particles::ParticleHandler<dim>           particle_handler;


    DoFHandler<dim> fluid_dh;
    FESystem<dim>   fluid_fe;
    MappingQ<dim>   mapping;
    LA::MPI::Vector field_owned;
    LA::MPI::Vector field_relevant;

    SingleVortex<dim> velocity;

    ConditionalOStream pcout;

    bool interpolated_velocity = false;
  };

  template <int dim>
  ParticleTracking<dim>::ParticleTracking(const ParticleTrackingParameters &par,
                                          const bool interpolated_velocity)
    : par(par)
    , mpi_communicator(MPI_COMM_WORLD)
    , background_triangulation(MPI_COMM_WORLD)
    , fluid_dh(background_triangulation)
    , fluid_fe(FE_Q<dim>(par.velocity_degree), dim)
    , mapping(par.velocity_degree)
    , pcout({std::cout, Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0})
    , interpolated_velocity(interpolated_velocity)

  {}

      template <int dim>
    unsigned int
    ParticleTracking<dim>::cell_weight(const typename parallel::distributed::Triangulation<dim>::cell_iterator &cell,
                            const typename parallel::distributed::Triangulation<dim>::CellStatus status)
    {
      if (cell->is_active() && !cell->is_locally_owned())
        return 0;

      // This determines how important particle distribution is compared to cell distribution
      // (1 cell == 1000). We set this number much higher to indicate the particle load is the
      // only one that is important to distribute.
        const unsigned int particle_weight = 10000;

      if (status == parallel::distributed::Triangulation<dim>::CELL_PERSIST
          || status == parallel::distributed::Triangulation<dim>::CELL_REFINE)
        {
          const unsigned int n_particles_in_cell = particle_handler.n_particles_in_cell(cell);
          return n_particles_in_cell * particle_weight;
        }
      else if (status == parallel::distributed::Triangulation<dim>::CELL_COARSEN)
        {
          unsigned int n_particles_in_cell = 0;

          for (unsigned int child_index = 0; child_index < GeometryInfo<dim>::max_children_per_cell; ++child_index)
            n_particles_in_cell += particle_handler.n_particles_in_cell(cell->child(child_index));

          return n_particles_in_cell * particle_weight;
        }

      Assert (false, ExcInternalError());
      return 0;
    }

  // Generation of particles using the grid where particles are generated at the
  // locations of the degrees of freedom.
  template <int dim>
  void
  ParticleTracking<dim>::particles_generation()
  {
    // Create a square triangulation
    GridGenerator::hyper_cube(background_triangulation, 0, 1);
    background_triangulation.refine_global(par.fluid_refinement);

    // In order to consider the particles when repartitioning the triangulation
    // the algorithm needs to know three things:
    // 1. How much weight to assign to each cell (how many particles are in there)
    // 2. How to pack the particles before shipping data around
    // 3. How to unpack the particles after repartitioning
    // Attach the correct functions to the signals inside parallel::distributed::Triangulation,
    // which will be called every time the repartition() function is called.
                    background_triangulation.signals.cell_weight.connect(
          [&] (const typename parallel::distributed::Triangulation<dim>::cell_iterator &cell,
               const typename parallel::distributed::Triangulation<dim>::CellStatus status)
          -> unsigned int
        {
          return this->cell_weight(cell, status);
        });

            background_triangulation.signals.pre_distributed_repartition.connect(std::bind(
      &Particles::ParticleHandler<dim>::register_store_callback_function,
      &particle_handler));

          background_triangulation.signals.post_distributed_repartition.connect(std::bind(
      &Particles::ParticleHandler<dim>::register_load_callback_function,
      &particle_handler,
      false));

    // Establish where the particles are living
    particle_handler.initialize(background_triangulation, mapping);

    // Generate the necessary bounding boxes for the generator of the particles
    const auto my_bounding_box = GridTools::compute_mesh_predicate_bounding_box(
      background_triangulation, IteratorFilters::LocallyOwnedCell());
    const auto global_bounding_boxes =
      Utilities::MPI::all_gather(MPI_COMM_WORLD, my_bounding_box);

    Point<dim> center;
    center[0] = 0.5;
    center[1] = 0.75;
    if (dim == 3)
      {
        center[2] = 0.5;
      }

    const double outer_radius = 0.15;
    const double inner_radius = 0.001;

    // Generation and refinement of the grid where the particles will be
    // created.
    parallel::distributed::Triangulation<dim> particle_triangulation(
      MPI_COMM_WORLD);

    GridGenerator::hyper_shell(
      particle_triangulation, center, inner_radius, outer_radius, 6);
    particle_triangulation.refine_global(par.particle_insertion_refinement);

    DoFHandler<dim> particles_dof_handler(particle_triangulation);
    FE_Q<dim>       particles_fe(1);

    particles_dof_handler.distribute_dofs(particles_fe);

    // Generation of the particles using the Particles::Generators
    Particles::Generators::dof_support_points(particles_dof_handler,
                                              global_bounding_boxes,
                                              particle_handler);

    // Displaying the total number of generated particles in the domain
    pcout << "Number of particles inserted: "
          << particle_handler.n_global_particles() << std::endl;
  }

  // Sets up the background degree of freedom using their interpolation
  // And allocated a vector where you can store the entire solution
  // of the velocity field
  template <int dim>
  void
  ParticleTracking<dim>::setup_background_dofs()
  {
    fluid_dh.distribute_dofs(fluid_fe);
    IndexSet locally_owned_dofs = fluid_dh.locally_owned_dofs();
    IndexSet locally_relevant_dofs;
    DoFTools::extract_locally_relevant_dofs(fluid_dh, locally_relevant_dofs);

    field_owned.reinit(locally_owned_dofs, mpi_communicator);
    field_relevant.reinit(locally_owned_dofs,
                          locally_relevant_dofs,
                          mpi_communicator);

    pcout << "Number of degrees of freedom in background grid: "
          << fluid_dh.n_dofs() << std::endl;
  }

  template <int dim>
  void
  ParticleTracking<dim>::interpolate_function_to_field()
  {
    const MappingQ<dim> mapping(fluid_fe.degree);

    VectorTools::interpolate(mapping, fluid_dh, velocity, field_owned);
    field_relevant = field_owned;
  }

  template <int dim>
  void
  ParticleTracking<dim>::euler_interpolated(double dt)
  {
    std::vector<types::global_dof_index> dof_indices(fluid_fe.dofs_per_cell);

    Vector<double> dof_data_per_cell(fluid_fe.dofs_per_cell);

    Tensor<1, dim> particle_velocity;


    auto particle = particle_handler.begin();
    while (particle != particle_handler.end())
      {
        const auto &cell =
          particle->get_surrounding_cell(background_triangulation);
        const auto &dh_cell =
          typename DoFHandler<dim>::cell_iterator(*cell, &fluid_dh);
        dh_cell->get_dof_indices(dof_indices);

        // Gather the DOF information in a local vector to prevent dynamically
        // re-accessing everything when there are multiple particles in a cell
        for (unsigned int j = 0; j < fluid_fe.dofs_per_cell; ++j)
          {
            dof_data_per_cell[j] = field_relevant(dof_indices[j]);
          }

        const auto pic = particle_handler.particles_in_cell(cell);
        for (; particle != pic.end(); ++particle)
          {
            const auto &reference_location = particle->get_reference_location();
            particle_velocity              = 0.;
            for (unsigned int j = 0; j < fluid_fe.dofs_per_cell; ++j)
              {
                const auto comp_j = fluid_fe.system_to_component_index(j);

                particle_velocity[comp_j.first] +=
                  fluid_fe.shape_value(j, reference_location) *
                  dof_data_per_cell[j];
              }

            Point<dim> particle_location = particle->get_location();
            for (int d = 0; d < dim; ++d)
              particle_location[d] += particle_velocity[d] * dt;
            particle->set_location(particle_location);
          }
      }
  }

  template <int dim>
  void
  ParticleTracking<dim>::euler_analytical(double dt)
  {
    Vector<double> particle_velocity(dim);

    // Looping over all particles in the domain using a particle iterator
    for (auto particle = particle_handler.begin();
         particle != particle_handler.end();
         ++particle)
      {
        // Get the velocity using the current location of particle
        velocity.vector_value(particle->get_location(), particle_velocity);

        Point<dim> particle_location = particle->get_location();
        // Updating the position of the particles and Setting the old position
        // equal to the new position of the particle
        for (int d = 0; d < dim; ++d)
          particle_location[d] += particle_velocity[d] * dt;

        particle->set_location(particle_location);
      }
  }

  //  template <int dim>
  //  void
  //  ParticleTracking<dim>::parallel_weight()
  //  {
  //    parallel::CellWeights<dim> cell_weights(background_dh);

  //  }

  template <int dim>
  void
  ParticleTracking<dim>::output_particles(unsigned int it)
  {
    Particles::DataOut<dim, dim> particle_output;
    particle_output.build_patches(particle_handler);
    std::string output_folder(par.output_directory);
    std::string file_name(interpolated_velocity ? "interpolated-particles" :
                                                  "analytical-particles");

                                                  pcout << "Writing particle output file: " << file_name << "-" << it << std::endl;

    particle_output.write_vtu_with_pvtu_record(
      output_folder, file_name, it, mpi_communicator, 6);
  }

  template <int dim>
  void
  ParticleTracking<dim>::output_background(unsigned int it)
  {
    std::vector<std::string> solution_names(dim, "velocity");
    std::vector<DataComponentInterpretation::DataComponentInterpretation>
      data_component_interpretation(
        dim, DataComponentInterpretation::component_is_part_of_vector);

    DataOut<dim> data_out;

    // Attach the solution data to data_out object
    data_out.attach_dof_handler(fluid_dh);
    data_out.add_data_vector(field_relevant,
                             solution_names,
                             DataOut<dim>::type_dof_data,
                             data_component_interpretation);
    Vector<float> subdomain(background_triangulation.n_active_cells());
    for (unsigned int i = 0; i < subdomain.size(); ++i)
      subdomain(i) = background_triangulation.locally_owned_subdomain();
    data_out.add_data_vector(subdomain, "subdomain");

    MappingQ<dim> mapping(fluid_fe.degree);

    data_out.build_patches(mapping);

    std::string output_folder(par.output_directory);
    std::string file_name("background");

    pcout << "Writing background field file: " << file_name << "-" << it << std::endl;

    data_out.write_vtu_with_pvtu_record(
      output_folder, file_name, it, mpi_communicator, 6);
  }

  template <int dim>
  void
  ParticleTracking<dim>::run_analytical_velocity()
  {
    DiscreteTime discrete_time(0, par.final_time, par.time_step);

    particles_generation();

    pcout << "Repartitioning triangulation after particle generation" << std::endl;
          background_triangulation.repartition();

    setup_background_dofs();
    interpolate_function_to_field();

    output_particles(discrete_time.get_step_number());
    output_background(discrete_time.get_step_number());

    // Looping over time in order to move the particles
    while (!discrete_time.is_at_end())
      {
        discrete_time.advance_time();
        velocity.set_time(discrete_time.get_previous_time());

        if ((discrete_time.get_step_number() % par.repartition_frequency) == 0)
          {
                        pcout << "Repartitioning triangulation after particle advection" << std::endl;
        background_triangulation.repartition();
            setup_background_dofs();
          }

        interpolate_function_to_field();

        if (interpolated_velocity)
          euler_interpolated(discrete_time.get_previous_step_size());
        else
          euler_analytical(discrete_time.get_previous_step_size());

        particle_handler.sort_particles_into_subdomains_and_cells();

        if ((discrete_time.get_step_number() % par.output_frequency) == 0)
          {
            output_particles(discrete_time.get_step_number());
            output_background(discrete_time.get_step_number());
          }
      }
  }

} // namespace Stepx

// @sect3{The main() function}

// The remainder of the code, the `main()` function, is standard.
int
main(int argc, char *argv[])
{
  using namespace Stepx;
  using namespace dealii;
  deallog.depth_console(1);

  try
    {
      Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

      std::string prm_file;
      if (argc > 1)
        prm_file = argv[1];
      else
        prm_file = "parameters.prm";

      ParticleTrackingParameters par;
      ParameterAcceptor::initialize(prm_file);
      {
        Stepx::ParticleTracking<2> particle_tracking(par, false);
        particle_tracking.run_analytical_velocity();
      }
      {
        Stepx::ParticleTracking<2> particle_tracking(par, true);
        particle_tracking.run_analytical_velocity();
      }
    }
  catch (std::exception &exc)
    {
      std::cerr << std::endl
                << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Exception on processing: " << std::endl
                << exc.what() << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;

      return 1;
    }
  catch (...)
    {
      std::cerr << std::endl
                << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Unknown exception!" << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;
      return 1;
    }

  return 0;
}
