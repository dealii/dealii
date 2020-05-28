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
 * Authors: Bruno Blais, Toni El Geitani Nehme, Rene Gassmoeller, Peter Munch
 */


// @sect3{Include files}

// The majority of the include files are generic

#include <deal.II/base/bounding_box.h>
#include <deal.II/base/conditional_ostream.h>
#include <deal.II/base/discrete_time.h>
#include <deal.II/base/mpi.h>
#include <deal.II/base/parameter_acceptor.h>
#include <deal.II/base/timer.h>

#include <deal.II/distributed/cell_weights.h>
#include <deal.II/distributed/solution_transfer.h>
#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/mapping_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/lac/generic_linear_algebra.h>
#include <deal.II/lac/petsc_vector.h>
#include <deal.II/lac/trilinos_vector.h>
#include <deal.II/lac/vector.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/vector_tools.h>

// From the following include file we import the ParticleHandler class
// that allows you to manage
// a collection of particles (objects of type Particles::Particle), representing
// a collection of points with some attached properties (e.g., an id) floating
// on a parallel::distributed::Triangulation. The methods and classes in the
// namespace Particles allows one to easily implement Particle-In-Cell methods
// and particle tracing on distributed triangulations
#include <deal.II/particles/particle_handler.h>

// We import the particles generator
// which allow us to insert the particles. In the present step, the particle
// are globally inserted using a non-matching hyper-shell triangulation
#include <deal.II/particles/generators.h>

// Since the particles do not form a triangulation, they have their
// own specific data out class which will enable us to write them
// to commonly used parallel vtu format
#include <deal.II/particles/data_out.h>


// This step uses parallel vector to interpolate the velocity field
// at the position of the particles. This step supports the use of both
// Trilinos and PETSC distributed vectors
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
  // difficult to track down (such as writing things like `time = 0` instead of
  // `time == 0`), we declare all the parameters in an external class, which is
  // initialized before the actual `ParticleTracking` class, and pass it to
  // the main class as a `const` reference.
  //
  // The constructor of the class is responsible for the connection between the
  // members of this class and the corresponding entries in the
  // ParameterHandler. Thanks to the use of the
  // ParameterHandler::add_parameter() method, this connection is trivial, but
  // requires all members of this class to be writeable.
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

    unsigned int velocity_degree       = 1;
    double       time_step             = 0.002;
    double       final_time            = 4.0;
    unsigned int output_frequency      = 10;
    unsigned int repartition_frequency = 5;

    // We allow every grid to be refined independently. In this tutorial, no
    // physics is resolved on the fluid grid, and its velocity is calculated
    // analytically.
    unsigned int fluid_refinement              = 4;
    unsigned int particle_insertion_refinement = 3;
  };



  // There remains the task of declaring what run-time parameters we can accept
  // in input files. Since we have a very limited number of parameters, all
  // parameters are declared in the same category.
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


  // @sect3{Velocity profile}

  // The velocity profile is provided as a Function object. We provide the
  // velocity profile. In the present step, this function is hard-coded within
  // the example. However, it could have been easily made using a ParsedFunction
  template <int dim>
  class Vortex : public Function<dim>
  {
  public:
    Vortex()
      : Function<dim>(dim)
    {}
    virtual void vector_value(const Point<dim> &point,
                              Vector<double> &  values) const override;
  };

  template <int dim>
  void Vortex<dim>::vector_value(const Point<dim> &point,
                                 Vector<double> &  values) const
  {
    const double T = 4;
    // Since the velocity profile is time dependant, the present time in the
    // simulation must be gathered from the Function object.
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

  // @sect3{The <code>PatricleTracking</code> class declaration}

  // We are now ready to introduce the main class of our tutorial program.
  // Contrarily to some other steps, there is an additional function that is
  // left public other than the constructor and the `run()` method, which is the
  // cell_weight() function. This function is connected to the triangulation and
  // must be callable from outside of the scope of this class. Everything else
  // is left `private`, and accessed through the run method itself.
  template <int dim>
  class ParticleTracking
  {
  public:
    ParticleTracking(const ParticleTrackingParameters &par,
                     const bool                        interpolated_velocity);
    void run();

    // The cell_weight() function indicates to the triangulation how much
    // computational work is expected to happen on this cell, and consequently
    // how the domain needs to be partitioned so that every MPI rank receives a
    // roughly equal amount of work (potentially not an equal number of cells).
    unsigned int cell_weight(
      const typename parallel::distributed::Triangulation<dim>::cell_iterator
        &cell,
      const typename parallel::distributed::Triangulation<dim>::CellStatus
        status);

  private:
    // The particles_generation function is responsible for the initial
    // generation of the particles on top of the background grid
    void particles_generation();

    // When the velocity profile is interpolated to the position of the
    // particles, it must first be stored using degrees of freedom.
    // Consequently, as is the case for other parallel case (e.g. step-40) we
    // initialize the degrees of freedom on the background grid
    void setup_background_dofs();

    void interpolate_function_to_field();

    // The next two functions are responsible for carrying out explicit Euler
    // time integration for the cases where the velocity field is interpolated
    // at the positions of the particles or calculated analytically,
    // respectively
    void euler_interpolated(double dt);
    void euler_analytical(double dt);

    // The following two functions are responsible for outputting the simulation
    // results for the particles and for the velocity profile on the background
    // mesh, respectively.
    void output_particles(unsigned int it);
    void output_background(unsigned int it);

    // The private members of this class are similar to other parallel deal.II
    // examples. The parameters are stored as a const member. It is important
    // to note that we keep the Vortex class as a member since its time
    // must be modified as the simulation proceeds.

    const ParticleTrackingParameters &par;

    MPI_Comm                                  mpi_communicator;
    parallel::distributed::Triangulation<dim> background_triangulation;
    Particles::ParticleHandler<dim>           particle_handler;

    DoFHandler<dim> fluid_dh;
    FESystem<dim>   fluid_fe;
    MappingQ<dim>   mapping;
    LA::MPI::Vector field_owned;
    LA::MPI::Vector field_relevant;

    Vortex<dim> velocity;

    ConditionalOStream pcout;

    bool interpolated_velocity;
  };

  // @sect3{The <code>PatricleTracking</code> class implementation}

  // @sect4{Constructor}

  // Constructors and destructors are rather trivial. They are very similar
  // to what is done in step-40. we set the set of processors we want to work on
  // to all machines available (MPI_COMM_WORLD) and
  // initialize the <code>pcout</code> variable to only allow processor zero
  // to output anything to the standard output.

  template <int dim>
  ParticleTracking<dim>::ParticleTracking(const ParticleTrackingParameters &par,
                                          const bool interpolated_velocity)
    : par(par)
    , mpi_communicator(MPI_COMM_WORLD)
    , background_triangulation(mpi_communicator)
    , fluid_dh(background_triangulation)
    , fluid_fe(FE_Q<dim>(par.velocity_degree), dim)
    , mapping(par.velocity_degree)
    , pcout(
        {std::cout, Utilities::MPI::this_mpi_process(mpi_communicator) == 0})
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
  // connected to the cell_weight() signal inside the triangulation, and will be
  // called once per cell, whenever the triangulation repartitions the domain
  // between ranks (the connection is created inside the
  // particles_generation() function of this class).
  template <int dim>
  unsigned int ParticleTracking<dim>::cell_weight(
    const typename parallel::distributed::Triangulation<dim>::cell_iterator
      &                                                                  cell,
    const typename parallel::distributed::Triangulation<dim>::CellStatus status)
  {
    if (cell->is_active() && !cell->is_locally_owned())
      return 0;

    // This determines how important particle work is compared to cell
    // work (by default every cell has a weight of 1000). 
    // We set the weight per particle much higher to indicate that
    // the particle load is the only one that is important to distribute
    // in this example. The optimal value of this number depends on the
    // application and can range from 0 (cheap particle operations,
    // expensive cell operations) to much larger than 1000 (expensive
    // particle operations, cheap cell operations, like in this example).
    const unsigned int particle_weight = 10000;

    // This example does not use adaptive refinement, therefore every cell
    // should have the status CELL_PERSIST. However this function can also
    // be used to distribute load during refinement, therefore we consider
    // refined or coarsened cells as well.
    if (status == parallel::distributed::Triangulation<dim>::CELL_PERSIST ||
        status == parallel::distributed::Triangulation<dim>::CELL_REFINE)
      {
        const unsigned int n_particles_in_cell =
          particle_handler.n_particles_in_cell(cell);
        return n_particles_in_cell * particle_weight;
      }
    else if (status == parallel::distributed::Triangulation<dim>::CELL_COARSEN)
      {
        unsigned int n_particles_in_cell = 0;

        for (unsigned int child_index = 0;
             child_index < GeometryInfo<dim>::max_children_per_cell;
             ++child_index)
          n_particles_in_cell +=
            particle_handler.n_particles_in_cell(cell->child(child_index));

        return n_particles_in_cell * particle_weight;
      }

    Assert(false, ExcInternalError());
    return 0;
  }

  // @sect4{Particles generation}

  // This function generates the tracer particles and the background
  // triangulation on which these particles evolve.
  template <int dim>
  void ParticleTracking<dim>::particles_generation()
  {
    // We create an hyper_cube triangulation which we globally define. This
    // triangulation englobes the full trajectory of the particles.
    GridGenerator::hyper_cube(background_triangulation, 0, 1);
    background_triangulation.refine_global(par.fluid_refinement);

    // In order to consider the particles when repartitioning the triangulation
    // the algorithm needs to know three things:
    // 1. How much weight to assign to each cell (how many particles are in
    // there)
    // 2. How to pack the particles before shipping data around
    // 3. How to unpack the particles after repartitioning
    // Attach the correct functions to the signals inside
    // parallel::distributed::Triangulation, which will be called every time the
    // repartition() function is called.
    // These connections only need to be created once, so we might as well
    // have set them up in the constructor of this class, but for the purpose
    // of this example we want to group the particle related instructions.
    background_triangulation.signals.cell_weight.connect(
      [&](
        const typename parallel::distributed::Triangulation<dim>::cell_iterator
          &cell,
        const typename parallel::distributed::Triangulation<dim>::CellStatus
          status) -> unsigned int { return this->cell_weight(cell, status); });

    background_triangulation.signals.pre_distributed_repartition.connect(
      std::bind(
        &Particles::ParticleHandler<dim>::register_store_callback_function,
        &particle_handler));

    background_triangulation.signals.post_distributed_repartition.connect(
      std::bind(
        &Particles::ParticleHandler<dim>::register_load_callback_function,
        &particle_handler,
        false));

    // Establish the background triangulation where the particles are living
    particle_handler.initialize(background_triangulation, mapping);

    // We create a particle triangulation which is solely used to generate
    // the points which will be used to insert the particles. This
    // triangulation is an hyper_shell which is off-set from the
    // center of the simulation domain.

    Point<dim> center;
    center[0] = 0.5;
    center[1] = 0.75;
    if (dim == 3)
      {
        center[2] = 0.5;
      }

    const double outer_radius = 0.15;
    const double inner_radius = 0.01;

    parallel::distributed::Triangulation<dim> particle_triangulation(
      MPI_COMM_WORLD);

    GridGenerator::hyper_shell(
      particle_triangulation, center, inner_radius, outer_radius, 6);
    particle_triangulation.refine_global(par.particle_insertion_refinement);

    DoFHandler<dim> particles_dof_handler(particle_triangulation);
    FE_Q<dim>       particles_fe(1);

    particles_dof_handler.distribute_dofs(particles_fe);

    // We generate the necessary bounding boxes for the particles generator
    // These bounding boxes are required to quickly identify in which
    // processors and which cell the inserted particle lies.
    const auto my_bounding_box = GridTools::compute_mesh_predicate_bounding_box(
      background_triangulation, IteratorFilters::LocallyOwnedCell());
    const auto global_bounding_boxes =
      Utilities::MPI::all_gather(MPI_COMM_WORLD, my_bounding_box);

    // We generate the particles at the position of the degree of
    // freedom of the dummy particle triangulation
    Particles::Generators::dof_support_points(particles_dof_handler,
                                              global_bounding_boxes,
                                              particle_handler);

    // Displaying the total number of generated particles in the domain
    pcout << "Number of particles inserted: "
          << particle_handler.n_global_particles() << std::endl;
  }

  // @sect4{Background DOFs and interpolation}

  // Sets up the background degree of freedom used for the velocity
  // interpolation And allocate the field vector where the entire
  // solution of the velocity field is stored
  template <int dim>
  void ParticleTracking<dim>::setup_background_dofs()
  {
    fluid_dh.distribute_dofs(fluid_fe);
    IndexSet locally_owned_dofs = fluid_dh.locally_owned_dofs();
    IndexSet locally_relevant_dofs;
    DoFTools::extract_locally_relevant_dofs(fluid_dh, locally_relevant_dofs);

    field_owned.reinit(locally_owned_dofs, mpi_communicator);
    field_relevant.reinit(locally_owned_dofs,
                          locally_relevant_dofs,
                          mpi_communicator);
  }

  // Interpolates the Vortex velocity field to the field vector
  template <int dim>
  void ParticleTracking<dim>::interpolate_function_to_field()
  {
    const MappingQ<dim> mapping(fluid_fe.degree);

    VectorTools::interpolate(mapping, fluid_dh, velocity, field_owned);
    field_relevant = field_owned;
  }

  // @sect4{Time integration of the trajectories}

  // We integrate the particle trajectories
  // using an analytically defined velocity field. This is a relatively trivial
  // usage of the particles.
  template <int dim>
  void ParticleTracking<dim>::euler_analytical(double dt)
  {
    Vector<double> particle_velocity(dim);

    // Looping over all particles in the domain using a particle iterator
    for (auto particle = particle_handler.begin();
         particle != particle_handler.end();
         ++particle)
      {
        // Get the velocity using the current location of particle
        Point<dim> particle_location = particle->get_location();
        velocity.vector_value(particle_location, particle_velocity);

        // Updating the position of the particles and Setting the old position
        // equal to the new position of the particle
        for (int d = 0; d < dim; ++d)
          particle_location[d] += particle_velocity[d] * dt;

        particle->set_location(particle_location);
      }
  }


  // We integrate the particle trajectories by interpolating the value of the
  // velocity field at the degrees of freedom to the position of the particles.
  template <int dim>
  void ParticleTracking<dim>::euler_interpolated(double dt)
  {
    std::vector<types::global_dof_index> dof_indices(fluid_fe.dofs_per_cell);
    Vector<double> dof_data_per_cell(fluid_fe.dofs_per_cell);

    // We loop over all the local particles. Although this could be achieved
    // directly by looping over all the cells, this would force us
    // to loop over numerous cells which do not contain particles.
    // Consequently, we loop over all the particles, but, we get the reference
    // of the cell in which the particle lies and then loop over all particles
    // within that cell. This enables us to gather the values of the velocity
    // out of the field_relevant vector once and use them for all particles
    // that lie within the cell. Once we are done with all particles on one
    // cell, we advance the `particle` iterator to the particle past the end of
    // the ones on the current cell (this is the last line of the `while` loop's
    // body).
    auto particle = particle_handler.begin();
    while (particle != particle_handler.end())
      {
        const auto &cell =
          particle->get_surrounding_cell(background_triangulation);
        const auto &dh_cell =
          typename DoFHandler<dim>::cell_iterator(*cell, &fluid_dh);
        dh_cell->get_dof_indices(dof_indices);

        // Gather the velocity information in a local vector to prevent
        // dynamically re-accessing everything when there are multiple particles
        // in a cell
        for (unsigned int j = 0; j < fluid_fe.dofs_per_cell; ++j)
          {
            dof_data_per_cell[j] = field_relevant(dof_indices[j]);
          }

        // Compute the velocity at the particle locations by evaluating
        // the finite element solution at the position of the particles.
        // This is essentially an optimized version of the particle advection
        // functionality in step-19, but instead of creating quadrature
        // objects and FEValues objects for each cell, we do the
        // evaluation by hand, which is somewhat more efficient and only
        // matters for this tutorial, because the particle work is the
        // dominant cost of the whole program.
        const auto pic = particle_handler.particles_in_cell(cell);
        for (; particle != pic.end(); ++particle)
          {
            const auto &reference_location = particle->get_reference_location();
            Tensor<1, dim> particle_velocity;
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


  // @sect4{Data output}

  // These two functions take care of writing both the particles
  // and the background mesh to vtu with a pvtu record

  template <int dim>
  void ParticleTracking<dim>::output_particles(unsigned int it)
  {
    Particles::DataOut<dim, dim> particle_output;
    particle_output.build_patches(particle_handler);
    std::string output_folder(par.output_directory);
    std::string file_name(interpolated_velocity ? "interpolated-particles" :
                                                  "analytical-particles");

    pcout << "Writing particle output file: " << file_name << "-" << it
          << std::endl;

    particle_output.write_vtu_with_pvtu_record(
      output_folder, file_name, it, mpi_communicator, 6);
  }

  template <int dim>
  void ParticleTracking<dim>::output_background(unsigned int it)
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

    pcout << "Writing background field file: " << file_name << "-" << it
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

    particles_generation();

    pcout << "Repartitioning triangulation after particle generation"
          << std::endl;
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
