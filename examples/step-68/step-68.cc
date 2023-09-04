/* ---------------------------------------------------------------------
 *
 * Copyright (C) 2020 - 2023 by the deal.II authors
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
 * Authors: Bruno Blais, Toni El Geitani Nehme, Rene Gassmoeller, Peter Munch
 */

// @sect3{Include files}

#include <deal.II/base/bounding_box.h>
#include <deal.II/base/conditional_ostream.h>
#include <deal.II/base/discrete_time.h>
#include <deal.II/base/function_lib.h>
#include <deal.II/base/mpi.h>
#include <deal.II/base/parameter_acceptor.h>
#include <deal.II/base/timer.h>

#include <deal.II/distributed/solution_transfer.h>
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

// From the following include file we import the ParticleHandler class
// that allows you to manage
// a collection of particles (objects of type Particles::Particle), representing
// a collection of points with some attached properties (e.g., an id) floating
// on a parallel::distributed::Triangulation. The methods and classes in the
// namespace Particles allows one to easily implement Particle-In-Cell methods
// and particle tracing on distributed triangulations:
#include <deal.II/particles/particle_handler.h>

// We import the particles generator
// which allow us to insert the particles. In the present step, the particle
// are globally inserted using a non-matching hyper-shell triangulation:
#include <deal.II/particles/generators.h>

// Since the particles do not form a triangulation, they have their
// own specific DataOut class which will enable us to write them
// to commonly used parallel vtu format (or any number of other file formats):
#include <deal.II/particles/data_out.h>

#include <cmath>
#include <iostream>



namespace Step68
{
  using namespace dealii;

  // @sect3{Run-time parameter handling}

  // Similarly to what is done in step-60, we set up a class that holds
  // all the parameters of our problem and derive it from the ParameterAcceptor
  // class to simplify the management and creation of parameter files.
  //
  // The ParameterAcceptor paradigm requires all parameters to be writable by
  // the ParameterAcceptor methods. In order to avoid bugs that would be very
  // difficult to track down (such as writing things like `if (time = 0)`
  // instead of `if(time == 0)`), we declare all the parameters in an external
  // class, which is initialized before the actual `ParticleTracking` class, and
  // pass it to the main class as a `const` reference.
  //
  // The constructor of the class is responsible for the connection between the
  // members of this class and the corresponding entries in the
  // ParameterHandler. Thanks to the use of the
  // ParameterHandler::add_parameter() method, this connection is trivial, but
  // requires all members of this class to be writable.
  class ParticleTrackingParameters : public ParameterAcceptor
  {
  public:
    ParticleTrackingParameters();

    // This class consists largely of member variables that
    // describe the details of the particle tracking simulation and its
    // discretization. The following parameters are about where output should
    // written to, the spatial discretization of the velocity (the default is
    // $Q_1$), the time step and the output interval (how many time
    // steps should elapse before we generate graphical output again):
    std::string output_directory = "./";

    unsigned int velocity_degree      = 1;
    double       time_step            = 0.002;
    double       final_time           = 4.0;
    unsigned int output_interval      = 10;
    unsigned int repartition_interval = 5;

    // We allow every grid to be refined independently. In this tutorial, no
    // physics is resolved on the fluid grid, and its velocity is calculated
    // analytically.
    unsigned int fluid_refinement              = 4;
    unsigned int particle_insertion_refinement = 3;
  };



  // There remains the task of declaring what run-time parameters we can accept
  // in input files. Since we have a very limited number of parameters, all
  // parameters are declared in the same section.
  ParticleTrackingParameters::ParticleTrackingParameters()
    : ParameterAcceptor("Particle Tracking Problem/")
  {
    add_parameter(
      "Velocity degree", velocity_degree, "", prm, Patterns::Integer(1));

    add_parameter("Output interval",
                  output_interval,
                  "Iteration interval between which output results are written",
                  prm,
                  Patterns::Integer(1));

    add_parameter("Repartition interval",
                  repartition_interval,
                  "Iteration interval at which the mesh is load balanced",
                  prm,
                  Patterns::Integer(1));

    add_parameter("Output directory", output_directory);

    add_parameter("Time step", time_step, "", prm, Patterns::Double());

    add_parameter("Final time",
                  final_time,
                  "End time of the simulation",
                  prm,
                  Patterns::Double());

    add_parameter("Fluid refinement",
                  fluid_refinement,
                  "Refinement level of the fluid domain",
                  prm,
                  Patterns::Integer(0));

    add_parameter(
      "Particle insertion refinement",
      particle_insertion_refinement,
      "Refinement of the volumetric mesh used to insert the particles",
      prm,
      Patterns::Integer(0));
  }



  // @sect3{The <code>ParticleTracking</code> class declaration}

  // We are now ready to introduce the main class of our tutorial program.
  template <int dim>
  class ParticleTracking
  {
  public:
    ParticleTracking(const ParticleTrackingParameters &par,
                     const bool                        interpolated_velocity);
    void run();

  private:
    // This function is responsible for the initial
    // generation of the particles on top of the background grid.
    void generate_particles();

    // When the velocity profile is interpolated to the position of the
    // particles, it must first be stored using degrees of freedom.
    // Consequently, as is the case for other parallel case (e.g. step-40) we
    // initialize the degrees of freedom on the background grid.
    void setup_background_dofs();

    // In one of the test cases, the function is mapped to the background grid
    // and a finite element interpolation is used to calculate the velocity
    // at the particle location. This function calculates the value of the
    // function at the support point of the triangulation.
    void interpolate_function_to_field();

    // The next two functions are responsible for carrying out step of explicit
    // Euler time integration for the cases where the velocity field is
    // interpolated at the positions of the particles or calculated
    // analytically, respectively.
    void euler_step_interpolated(const double dt);
    void euler_step_analytical(const double dt);

    // The `cell_weight()` function indicates to the triangulation how much
    // computational work is expected to happen on this cell, and consequently
    // how the domain needs to be partitioned so that every MPI rank receives a
    // roughly equal amount of work (potentially not an equal number of cells).
    // While the function is called from the outside, it is connected to the
    // corresponding signal from inside this class, therefore it can be
    // `private`.
    unsigned int cell_weight(
      const typename parallel::distributed::Triangulation<dim>::cell_iterator
                      &cell,
      const CellStatus status) const;

    // The following two functions are responsible for outputting the simulation
    // results for the particles and for the velocity profile on the background
    // mesh, respectively.
    void output_particles(const unsigned int it);
    void output_background(const unsigned int it);

    // The private members of this class are similar to other parallel deal.II
    // examples. The parameters are stored as a `const` member. It is important
    // to note that we keep the `Vortex` class as a member since its time
    // must be modified as the simulation proceeds.

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



  // @sect3{The <code>PatricleTracking</code> class implementation}

  // @sect4{Constructor}

  // The constructors and destructors are rather trivial. They are very similar
  // to what is done in step-40. We set the processors we want to work on
  // to all machines available (`MPI_COMM_WORLD`) and
  // initialize the <code>pcout</code> variable to only allow processor zero
  // to output anything to the standard output.

  template <int dim>
  ParticleTracking<dim>::ParticleTracking(const ParticleTrackingParameters &par,
                                          const bool interpolated_velocity)
    : par(par)
    , mpi_communicator(MPI_COMM_WORLD)
    , background_triangulation(mpi_communicator)
    , fluid_dh(background_triangulation)
    , fluid_fe(FE_Q<dim>(par.velocity_degree) ^ dim)
    , velocity(4.0)
    , pcout(std::cout, Utilities::MPI::this_mpi_process(mpi_communicator) == 0)
    , interpolated_velocity(interpolated_velocity)

  {}



  // @sect4{Cell weight}

  // This function is the key component that allow us to dynamically balance the
  // computational load for this example. The function attributes a weight to
  // every cell that represents the computational work on this cell. Here the
  // majority of work is expected to happen on the particles, therefore the
  // return value of this function (representing "work for this cell") is
  // calculated based on the number of particles in the current cell.
  // The function is
  // connected to the `weight` signal inside the triangulation, and will be
  // called once per cell, whenever the triangulation repartitions the domain
  // between ranks (the connection is created inside the
  // generate_particles() function of this class).
  template <int dim>
  unsigned int ParticleTracking<dim>::cell_weight(
    const typename parallel::distributed::Triangulation<dim>::cell_iterator
                    &cell,
    const CellStatus status) const
  {
    // First, we introduce a base weight that will be assigned to every cell.
    const unsigned int base_weight = 1;

    // The following variable then determines how important particle work is
    // compared to cell work. We set the weight per particle much higher to
    // indicate that the particle load is the only one that is important to
    // distribute the cells in this example. The optimal value of this number
    // depends on the application and can range from 0 (cheap particle
    // operations, expensive cell operations) to much larger than the base
    // weight of 1 (expensive particle operations, cheap cell operations, like
    // presumed in this example).
    const unsigned int particle_weight = 10;

    // This example does not use adaptive refinement, therefore every cell
    // should have the status `CellStatus::cell_will_persist`. However this
    // function can also be used to distribute load during refinement, therefore
    // we consider refined or coarsened cells as well.
    unsigned int n_particles_in_cell = 0;
    switch (status)
      {
        case CellStatus::cell_will_persist:
        case CellStatus::cell_will_be_refined:
          n_particles_in_cell = particle_handler.n_particles_in_cell(cell);
          break;

        case CellStatus::cell_invalid:
          break;

        case CellStatus::children_will_be_coarsened:
          for (const auto &child : cell->child_iterators())
            n_particles_in_cell += particle_handler.n_particles_in_cell(child);
          break;

        default:
          Assert(false, ExcInternalError());
          break;
      }

    return base_weight + particle_weight * n_particles_in_cell;
  }



  // @sect4{Particles generation}

  // This function generates the tracer particles and the background
  // triangulation on which these particles evolve.
  template <int dim>
  void ParticleTracking<dim>::generate_particles()
  {
    // We create a hyper cube triangulation which we globally refine. This
    // triangulation covers the full trajectory of the particles.
    GridGenerator::hyper_cube(background_triangulation, 0, 1);
    background_triangulation.refine_global(par.fluid_refinement);

    // In order to consider the particles when repartitioning the triangulation
    // the algorithm needs to know how much weight to assign to each cell
    // (how many particles are in there).
    //
    // We attach a weight function to the signal inside
    // parallel::distributed::Triangulation. This signal will be called every
    // time the repartition() function is called. This connection only needs to
    // be created once, so we might as well have set it up in the constructor
    // of this class, but for the purpose of this example we want to group the
    // particle related instructions.
    background_triangulation.signals.weight.connect(
      [&](const typename parallel::distributed::Triangulation<
            dim>::cell_iterator &cell,
          const CellStatus       status) -> unsigned int {
        return this->cell_weight(cell, status);
      });

    // This initializes the background triangulation where the particles are
    // living and the number of properties of the particles.
    particle_handler.initialize(background_triangulation, mapping, 1 + dim);

    // We create a particle triangulation which is solely used to generate
    // the points which will be used to insert the particles. This
    // triangulation is a hyper shell which is offset from the
    // center of the simulation domain. This will be used to generate a
    // disk filled with particles which will allow an easy monitoring
    // of the motion due to the vortex.
    Point<dim> center;
    center[0] = 0.5;
    center[1] = 0.75;
    if (dim == 3)
      center[2] = 0.5;

    const double outer_radius = 0.15;
    const double inner_radius = 0.01;

    parallel::distributed::Triangulation<dim> particle_triangulation(
      MPI_COMM_WORLD);

    GridGenerator::hyper_shell(
      particle_triangulation, center, inner_radius, outer_radius, 6);
    particle_triangulation.refine_global(par.particle_insertion_refinement);

    // We generate the necessary bounding boxes for the particles generator.
    // These bounding boxes are required to quickly identify in which
    // process's subdomain the inserted particle lies, and which cell owns it.
    const auto my_bounding_box = GridTools::compute_mesh_predicate_bounding_box(
      background_triangulation, IteratorFilters::LocallyOwnedCell());
    const auto global_bounding_boxes =
      Utilities::MPI::all_gather(MPI_COMM_WORLD, my_bounding_box);

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

    pcout << "Number of particles inserted: "
          << particle_handler.n_global_particles() << std::endl;
  }



  // @sect4{Background DOFs and interpolation}

  // This function sets up the background degrees of freedom used for the
  // velocity interpolation and allocates the field vector where the entire
  // solution of the velocity field is stored.
  template <int dim>
  void ParticleTracking<dim>::setup_background_dofs()
  {
    fluid_dh.distribute_dofs(fluid_fe);
    const IndexSet locally_owned_dofs = fluid_dh.locally_owned_dofs();
    const IndexSet locally_relevant_dofs =
      DoFTools::extract_locally_relevant_dofs(fluid_dh);

    velocity_field.reinit(locally_owned_dofs,
                          locally_relevant_dofs,
                          mpi_communicator);
  }



  // This function takes care of interpolating the
  // vortex velocity field to the field vector. This is achieved rather easily
  // by using the VectorTools::interpolate() function.
  template <int dim>
  void ParticleTracking<dim>::interpolate_function_to_field()
  {
    velocity_field.zero_out_ghost_values();
    VectorTools::interpolate(mapping, fluid_dh, velocity, velocity_field);
    velocity_field.update_ghost_values();
  }



  // @sect4{Time integration of the trajectories}

  // We integrate the particle trajectories
  // using an analytically defined velocity field. This demonstrates a
  // relatively trivial usage of the particles.
  template <int dim>
  void ParticleTracking<dim>::euler_step_analytical(const double dt)
  {
    const unsigned int this_mpi_rank =
      Utilities::MPI::this_mpi_process(mpi_communicator);
    Vector<double> particle_velocity(dim);

    // Looping over all particles in the domain using a particle iterator
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
        // vector) in the particle properties. In this example, this is done
        // purely for visualization purposes.
        ArrayView<double> properties = particle.get_properties();
        for (int d = 0; d < dim; ++d)
          properties[d] = particle_velocity[d];
        properties[dim] = this_mpi_rank;
      }
  }



  // In contrast to the previous function in this function we
  // integrate the particle trajectories by interpolating the value of
  // the velocity field at the degrees of freedom to the position of
  // the particles. This is achieved using the FEPointEvaluation object.
  template <int dim>
  void ParticleTracking<dim>::euler_step_interpolated(const double dt)
  {
    Vector<double> local_dof_values(fluid_fe.dofs_per_cell);

    FEPointEvaluation<dim, dim> evaluator(mapping, fluid_fe, update_values);
    std::vector<Point<dim>>     particle_positions;

    // We loop over all the local particles. Although this could be achieved
    // directly by looping over all the cells, this would force us
    // to loop over numerous cells which do not contain particles.
    // Rather, we loop over all the particles, but, we get the reference
    // of the cell in which the particle lies and then loop over all particles
    // within that cell. This enables us to gather the values of the velocity
    // out of the `velocity_field` vector once and use them for all particles
    // that lie within the cell.
    auto particle = particle_handler.begin();
    while (particle != particle_handler.end())
      {
        const auto cell = particle->get_surrounding_cell();
        const auto dh_cell =
          typename DoFHandler<dim>::cell_iterator(*cell, &fluid_dh);

        dh_cell->get_dof_values(velocity_field, local_dof_values);

        // Next, compute the velocity at the particle locations by evaluating
        // the finite element solution at the position of the particles.
        // This is achieved using FEPointEvaluation object.
        const auto pic = particle_handler.particles_in_cell(cell);
        Assert(pic.begin() == particle, ExcInternalError());
        particle_positions.clear();
        for (auto &p : pic)
          particle_positions.push_back(p.get_reference_location());

        evaluator.reinit(cell, particle_positions);
        evaluator.evaluate(make_array_view(local_dof_values),
                           EvaluationFlags::values);

        // We move the particles using the interpolated velocity field
        for (unsigned int particle_index = 0; particle != pic.end();
             ++particle, ++particle_index)
          {
            Point<dim>            particle_location = particle->get_location();
            const Tensor<1, dim> &particle_velocity =
              evaluator.get_value(particle_index);
            particle_location += particle_velocity * dt;
            particle->set_location(particle_location);

            // Again, we store the particle velocity and the processor id in the
            // particle properties for visualization purposes.
            ArrayView<double> properties = particle->get_properties();
            for (int d = 0; d < dim; ++d)
              properties[d] = particle_velocity[d];

            properties[dim] =
              Utilities::MPI::this_mpi_process(mpi_communicator);
          }
      }
  }



  // @sect4{Data output}

  // The next two functions take care of writing both the particles
  // and the background mesh to vtu with a pvtu record. This ensures
  // that the simulation results can be visualized when the simulation is
  // launched in parallel.
  template <int dim>
  void ParticleTracking<dim>::output_particles(const unsigned int it)
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
    const std::string output_folder(par.output_directory);
    const std::string file_name(interpolated_velocity ?
                                  "interpolated-particles" :
                                  "analytical-particles");

    pcout << "Writing particle output file: " << file_name << '-' << it
          << std::endl;

    particle_output.write_vtu_with_pvtu_record(
      output_folder, file_name, it, mpi_communicator, 6);
  }



  template <int dim>
  void ParticleTracking<dim>::output_background(const unsigned int it)
  {
    std::vector<std::string> solution_names(dim, "velocity");
    std::vector<DataComponentInterpretation::DataComponentInterpretation>
      data_component_interpretation(
        dim, DataComponentInterpretation::component_is_part_of_vector);

    DataOut<dim> data_out;

    // Attach the solution data to data_out object
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

    const std::string output_folder(par.output_directory);
    const std::string file_name("background");

    pcout << "Writing background field file: " << file_name << '-' << it
          << std::endl;

    data_out.write_vtu_with_pvtu_record(
      output_folder, file_name, it, mpi_communicator, 6);
  }



  // @sect4{Running the simulation}
  // This function orchestrates the entire simulation. It is very similar
  // to the other time dependent tutorial programs -- take step-21 or step-26 as
  // an example. Note that we use the DiscreteTime class to monitor the time,
  // the time-step and the step-number. This function is relatively
  // straightforward.

  template <int dim>
  void ParticleTracking<dim>::run()
  {
    DiscreteTime discrete_time(0, par.final_time, par.time_step);

    generate_particles();

    pcout << "Repartitioning triangulation after particle generation"
          << std::endl;

    particle_handler.prepare_for_coarsening_and_refinement();
    background_triangulation.repartition();
    particle_handler.unpack_after_coarsening_and_refinement();

    // We set the initial property of the particles by doing an
    // explicit Euler iteration with a time-step of 0 both in the case
    // of the analytical and the interpolated approach.
    if (interpolated_velocity)
      {
        setup_background_dofs();
        interpolate_function_to_field();
        euler_step_interpolated(0.);
      }
    else
      euler_step_analytical(0.);

    output_particles(discrete_time.get_step_number());
    if (interpolated_velocity)
      output_background(discrete_time.get_step_number());

    // The particles are advected by looping over time.
    while (!discrete_time.is_at_end())
      {
        discrete_time.advance_time();
        velocity.set_time(discrete_time.get_previous_time());

        if ((discrete_time.get_step_number() % par.repartition_interval) == 0)
          {
            particle_handler.prepare_for_coarsening_and_refinement();
            background_triangulation.repartition();
            particle_handler.unpack_after_coarsening_and_refinement();

            if (interpolated_velocity)
              setup_background_dofs();
          }

        if (interpolated_velocity)
          {
            interpolate_function_to_field();
            euler_step_interpolated(discrete_time.get_previous_step_size());
          }
        else
          euler_step_analytical(discrete_time.get_previous_step_size());

        // After the particles have been moved, it is necessary to identify
        // in which cell they now reside. This is achieved by calling
        // <code>sort_particles_into_subdomains_and_cells</code>
        particle_handler.sort_particles_into_subdomains_and_cells();

        if ((discrete_time.get_step_number() % par.output_interval) == 0)
          {
            output_particles(discrete_time.get_step_number());
            if (interpolated_velocity)
              output_background(discrete_time.get_step_number());
          }
      }
  }

} // namespace Step68



// @sect3{The main() function}

// The remainder of the code, the `main()` function, is standard.
// We note that we run the particle tracking with the analytical velocity
// and the interpolated velocity and produce both results
int main(int argc, char *argv[])
{
  using namespace Step68;
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
        Step68::ParticleTracking<2> particle_tracking(par, false);
        particle_tracking.run();
      }
      {
        Step68::ParticleTracking<2> particle_tracking(par, true);
        particle_tracking.run();
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
