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
 * Authors: Luca Heltai, Bruno Blais, Rene Gassmoeller, 2020
 */

// @sect3{Include files}
// Most of these have been introduced elsewhere, we'll comment only on the new
// ones. The switches close to the top that allow selecting between PETSc
// and Trilinos linear algebra capabilities are similar to the ones in
// step-40 and step-50.

#include <deal.II/base/function.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/timer.h>

#include <deal.II/lac/block_linear_operator.h>
#include <deal.II/lac/generic_linear_algebra.h>
#include <deal.II/lac/linear_operator.h>
#include <deal.II/lac/linear_operator_tools.h>

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

#include <deal.II/base/conditional_ostream.h>
#include <deal.II/base/index_set.h>
#include <deal.II/base/parameter_acceptor.h>
#include <deal.II/base/parsed_function.h>
#include <deal.II/base/utilities.h>

#include <deal.II/distributed/grid_refinement.h>
#include <deal.II/distributed/solution_transfer.h>
#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_renumbering.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_nothing.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_fe_field.h>
#include <deal.II/fe/mapping_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/manifold_lib.h>

#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/petsc_precondition.h>
#include <deal.II/lac/petsc_solver.h>
#include <deal.II/lac/petsc_sparse_matrix.h>
#include <deal.II/lac/petsc_vector.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/solver_gmres.h>
#include <deal.II/lac/solver_minres.h>
#include <deal.II/lac/sparsity_tools.h>
#include <deal.II/lac/vector.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/error_estimator.h>
#include <deal.II/numerics/vector_tools.h>

// These are the only new include files with regard to step-60. In this
// tutorial, the non-matching coupling between the solid and the fluid is
// computed using an intermediate data structure that keeps track of how the
// locations of quadrature points of the solid evolve within the fluid mesh.
// This data structure needs to keep track of the position of the quadrature
// points on each cell describing the solid domain, of the quadrature weights,
// and possibly of the normal vector to each point, if the solid domain is of
// co-dimension one.
//
// Deal.II offers these facilities in the Particles namespace, through the
// ParticleHandler class. ParticleHandler is a class that allows you to manage
// a collection of particles (objects of type Particles::Particle), representing
// a collection of points with some attached properties (e.g., an id) floating
// on a parallel::distributed::Triangulation. The methods and classes in the
// namespace Particles allows one to easily implement Particle-In-Cell methods
// and particle tracing on distributed triangulations.
//
// We "abuse" this data structure to store information about the location of
// solid quadrature points embedded in the surrounding fluid grid, including
// integration weights, and possibly surface normals. The reason why we use this
// additional data structure is related to the fact that the solid and the fluid
// grids might be non-overlapping, and if we were using two separate
// triangulation objects, would be distributed independently among parallel
// processes.
//
// In order to couple the two problems, we rely on the ParticleHandler class,
// storing in each particle the position of a solid quadrature point (which is
// in general not aligned to any of the fluid quadrature points), its weight,
// and any other information that may be required to couple the two problems.
// These locations are then propagated along with the (prescribed) velocity
// of the solid impeller.
//
// Ownership of the solid quadrature points is initially inherited from the MPI
// partitioning on the solid mesh itself. The Particles so generated are later
// distributed to the fluid mesh using the methods of the ParticleHandler class.
// This allows transparent exchange of information between MPI processes about
// the overlapping pattern between fluid cells and solid quadrature points.
#include <deal.II/particles/data_out.h>
#include <deal.II/particles/generators.h>
#include <deal.II/particles/particle_handler.h>
#include <deal.II/particles/utilities.h>

// When generating the grids, we allow reading it from a file, and if deal.II
// has been built with OpenCASCADE support, we also allow reading CAD files and
// use them as manifold descriptors for the grid (see step-54 for a detailed
// description of the various Manifold descriptors that are available in the
// OpenCASCADE namespace)
#include <deal.II/opencascade/manifold_lib.h>
#include <deal.II/opencascade/utilities.h>
#ifdef DEAL_II_WITH_OPENCASCADE
#  include <TopoDS.hxx>
#endif

#include <cmath>
#include <fstream>
#include <iostream>
#include <memory>

namespace Step70
{
  using namespace dealii;

  // @sect3{Run-time parameter handling}

  // Similarly to what we have done in step-60, we set up a class that holds
  // all the parameters of our problem and derive it from the ParameterAcceptor
  // class to simplify the management and creation of parameter files.
  //
  // The ParameterAcceptor paradigm requires all parameters to be writable by
  // the ParameterAcceptor methods. In order to avoid bugs that would be very
  // difficult to track down (such as writing things like `time = 0` instead of
  // `time == 0`), we declare all the parameters in an external class, which is
  // initialized before the actual `StokesImmersedProblem` class, and pass it to
  // the main class as a `const` reference.
  //
  // The constructor of the class is responsible for the connection between the
  // members of this class and the corresponding entries in the
  // ParameterHandler. Thanks to the use of the
  // ParameterHandler::add_parameter() method, this connection is trivial, but
  // requires all members of this class to be writeable.
  template <int dim, int spacedim = dim>
  class StokesImmersedProblemParameters : public ParameterAcceptor
  {
  public:
    StokesImmersedProblemParameters();

    // however, since this class will be passed as a `const` reference to the
    // StokesImmersedProblem class, we have to make sure we can still set the
    // time correctly in the objects derived by the Function class defined
    // herein. In order to do so, we declare both the
    // `StokesImmersedProblemParameters::rhs` and
    // `StokesImmersedProblemParameters::angular_velocity` members to be
    // `mutable`, and define the following little helper method that sets their
    // time to the correct value.
    void set_time(const double &time) const
    {
      rhs.set_time(time);
      angular_velocity.set_time(time);
    }

    // The remainder of the class consists largely of member variables that
    // describe the details of the simulation and its discretization. The
    // following parameters are about where output should land, the spatial and
    // temporal discretization (the default is the $Q_2\times Q_1$ Taylor-Hood
    // discretization which uses a polynomial degree of 2 for the velocity), and
    // how many time steps should elapse before we generate graphical output
    // again:
    std::string output_directory = ".";

    unsigned int velocity_degree = 2;

    unsigned int number_of_time_steps = 501;
    double       final_time           = 1.0;

    unsigned int output_frequency = 1;

    // We allow every grid to be refined independently. In this tutorial, no
    // physics is resolved on the solid grid, and its velocity is given as a
    // datum. However it is relatively straightforward to incorporate some
    // elasticity model in this tutorial, and transform it into a fully fledged
    // FSI solver.
    unsigned int initial_fluid_refinement      = 5;
    unsigned int initial_solid_refinement      = 5;
    unsigned int particle_insertion_refinement = 3;

    // To provide a rough description of the fluid domain, we use the method
    // extract_rtree_level() applied to the tree of bounding boxes of each
    // locally owned cell of the fluid triangulation. The higher the level of
    // the tree, the larger the number of extracted bounding boxes, and the more
    // accurate is the description of the fluid domain.
    // However, a large number of bounding boxes also implies a large
    // communication cost, since the collection of bounding boxes is gathered by
    // all processes.
    unsigned int fluid_rtree_extraction_level = 1;

    // The only two numerical parameters used in the equations are the viscosity
    // of the fluid, and the penalty term $\beta$ used in the Nitsche
    // formulation:
    double viscosity    = 1.0;
    double penalty_term = 100;

    // By default, we create a hyper_cube without colorization, and we use
    // homogeneous Dirichlet boundary conditions. In this set we store the
    // boundary ids to use when setting the boundary conditions:
    std::list<types::boundary_id> homogeneous_dirichlet_ids{0};

    // We illustrate here another way to create a Triangulation from a parameter
    // file, using the method GridGenerator::generate_from_name_and_arguments(),
    // that takes the name of a function in the GridGenerator namespace, and its
    // arguments as a single string representing the arguments as a tuple.
    //
    // The mechanism with which the arguments are parsed from and to a string is
    // explained in detail in the Patterns::Tools::Convert class, which is
    // used to translate from strings to most of the basic STL types (vectors,
    // maps, tuples) and basic deal.II types (Point, Tensor, BoundingBox, etc.).
    //
    // In general objects that can be represented by rank 1 uniform elements
    // (i.e., std::vector<double>, Point<dim>, std::set<int>, etc.) are comma
    // separated. Additional ranks take a semicolon, allowing you to parse
    // strings into objects of type `std::vector<std::vector<double>>`, or,
    // for example, `std::vector<Point<dim>>`, as `0.0, 0.1; 0.1, 0.2`. This
    // string could be interpreted as a vector of two Point objects, or a vector
    // of vector of doubles.
    //
    // When the entries are not uniform, as in the tuple case, we use a colon
    // to separate the various entries. For example, a string like `5: 0.1, 0.2`
    // could be used to parse an object of type `std::pair<int, Point<2>>` or a
    // `std::tuple<int, std::vector<double>>`.
    //
    // In our case most of the arguments are Point objects (representing
    // centers, corners, subdivision elements, etc.), integer values (number of
    // subdivisions), double values (radius, lengths, etc.), or boolean options
    // (such as the `colorize` option that many GridGenerator functions take).
    //
    // In the example below, we set reasonable default values, but these can be
    // changed at run time by selecting any other supported function of the
    // GridGenerator namespace. If the GridGenerator function fails, this
    // program will interpret the name of the grid as a vtk grid filename, and
    // the arguments as a map from manifold_id to the CAD files describing the
    // geometry of the domain. Every CAD file will be analyzed and a Manifold of
    // the OpenCASCADE namespace will be generated according to the content of
    // the CAD file itself.
    //
    // To be as generic as possible, we do this for each of the generated grids:
    // the fluid grid, the solid grid, but also the tracer particles which are
    // also generated using a triangulation.
    std::string name_of_fluid_grid       = "hyper_cube";
    std::string arguments_for_fluid_grid = "-1: 1: false";
    std::string name_of_solid_grid       = "hyper_rectangle";
    std::string arguments_for_solid_grid = spacedim == 2 ?
                                             "-.5, -.1: .5, .1: false" :
                                             "-.5, -.1, -.1: .5, .1, .1: false";
    std::string name_of_particle_grid = "hyper_ball";
    std::string arguments_for_particle_grid =
      spacedim == 2 ? "0.3, 0.3: 0.1: false" : "0.3, 0.3, 0.3 : 0.1: false";

    // Similarly, we allow for different local refinement strategies. In
    // particular, we limit the maximum number of refinement levels, in order
    // to control the minimum size of the fluid grid, and guarantee that it is
    // compatible with the solid grid. The minimum number of refinement levels
    // is also controlled to ensured sufficient accuracy in the
    // bulk of the flow. Additionally, we perform local refinement
    // based on standard error estimators on the fluid velocity field.
    //
    // We permit the user to choose between the
    // two most common refinement strategies, namely `fixed_number` or
    // `fixed_fraction`, that refer to the methods
    // GridRefinement::refine_and_coarsen_fixed_fraction() and
    // GridRefinement::refine_and_coarsen_fixed_number().
    //
    // Refinement may be done every few time steps, instead of continuously, and
    // we control this value by the `refinement_frequency` parameter:
    int          max_level_refinement = 8;
    int          min_level_refinement = 5;
    std::string  refinement_strategy  = "fixed_fraction";
    double       coarsening_fraction  = 0.3;
    double       refinement_fraction  = 0.3;
    unsigned int max_cells            = 20000;
    int          refinement_frequency = 5;

    // Finally, the following two function objects are used to control the
    // source term of Stokes flow and the angular velocity at which we move the
    // solid body. In a more realistic simulation, the solid velocity or its
    // deformation would come from the solution of an auxiliary problem on the
    // solid domain. In this example step we leave this part aside, and simply
    // impose a fixed rotational velocity field along the z-axis on the immersed
    // solid, governed by a function that can be specified in the parameter
    // file:
    mutable ParameterAcceptorProxy<Functions::ParsedFunction<spacedim>> rhs;
    mutable ParameterAcceptorProxy<Functions::ParsedFunction<spacedim>>
      angular_velocity;
  };



  // There remains the task of declaring what run-time parameters we can accept
  // in input files. We split the parameters in various categories, by putting
  // them in different sections of the ParameterHandler class. We begin by
  // declaring all the global parameters used by StokesImmersedProblem
  // in the global scope:
  template <int dim, int spacedim>
  StokesImmersedProblemParameters<dim,
                                  spacedim>::StokesImmersedProblemParameters()
    : ParameterAcceptor("Stokes Immersed Problem/")
    , rhs("Right hand side", spacedim + 1)
    , angular_velocity("Angular velocity")
  {
    add_parameter(
      "Velocity degree", velocity_degree, "", this->prm, Patterns::Integer(1));

    add_parameter("Number of time steps", number_of_time_steps);
    add_parameter("Output frequency", output_frequency);

    add_parameter("Output directory", output_directory);

    add_parameter("Final time", final_time);

    add_parameter("Viscosity", viscosity);

    add_parameter("Nitsche penalty term", penalty_term);

    add_parameter("Initial fluid refinement",
                  initial_fluid_refinement,
                  "Initial mesh refinement used for the fluid domain Omega");

    add_parameter("Initial solid refinement",
                  initial_solid_refinement,
                  "Initial mesh refinement used for the solid domain Gamma");

    add_parameter("Fluid bounding boxes extraction level",
                  fluid_rtree_extraction_level,
                  "Extraction level of the rtree used to construct global "
                  "bounding boxes");

    add_parameter(
      "Particle insertion refinement",
      particle_insertion_refinement,
      "Refinement of the volumetric mesh used to insert the particles");

    add_parameter(
      "Homogeneous Dirichlet boundary ids",
      homogeneous_dirichlet_ids,
      "Boundary Ids over which homogeneous Dirichlet boundary conditions are applied");

    // Next section is dedicated to the parameters used to create the
    // various grids. We will need three different triangulations: `Fluid
    // grid` is used to define the fluid domain, `Solid grid` defines the
    // solid domain, and `Particle grid` is used to distribute some tracer
    // particles, that are advected with the velocity and only used as
    // passive tracers.
    enter_subsection("Grid generation");
    {
      add_parameter("Fluid grid generator", name_of_fluid_grid);
      add_parameter("Fluid grid generator arguments", arguments_for_fluid_grid);

      add_parameter("Solid grid generator", name_of_solid_grid);
      add_parameter("Solid grid generator arguments", arguments_for_solid_grid);

      add_parameter("Particle grid generator", name_of_particle_grid);
      add_parameter("Particle grid generator arguments",
                    arguments_for_particle_grid);
    }
    leave_subsection();



    enter_subsection("Refinement and remeshing");
    {
      add_parameter("Refinement step frequency", refinement_frequency);
      add_parameter("Refinement maximal level", max_level_refinement);
      add_parameter("Refinement minimal level", min_level_refinement);
      add_parameter("Refinement strategy",
                    refinement_strategy,
                    "",
                    this->prm,
                    Patterns::Selection("fixed_fraction|fixed_number"));
      add_parameter("Refinement coarsening fraction", coarsening_fraction);
      add_parameter("Refinement fraction", refinement_fraction);
      add_parameter("Maximum number of cells", max_cells);
    }
    leave_subsection();

    // The final task is to correct the default dimension for the right hand
    // side function and define a meaningful default angular velocity instead of
    // zero.
    rhs.declare_parameters_call_back.connect([&]() {
      Functions::ParsedFunction<spacedim>::declare_parameters(this->prm,
                                                              spacedim + 1);
    });
    angular_velocity.declare_parameters_call_back.connect([&]() {
      this->prm.set("Function expression",
                    "t < .500001 ? 6.283185 : -6.283185");
    });
  }


  // Once the angular velocity is provided as a Function object, we reconstruct
  // the pointwise solid velocity through the following class which derives
  // from the Function class. It provides the value of the velocity of
  // the solid body at a given position by assuming that the body rotates
  // around the origin (or the $z$ axis in 3d) with a given angular velocity.
  template <int spacedim>
  class SolidVelocity : public Function<spacedim>
  {
  public:
    static_assert(spacedim > 1,
                  "Cannot instantiate SolidVelocity for spacedim == 1");

    SolidVelocity(const Functions::ParsedFunction<spacedim> &angular_velocity)
      : angular_velocity(angular_velocity)
    {}

    virtual double value(const Point<spacedim> &p,
                         unsigned int           component = 0) const override
    {
      Tensor<1, spacedim> velocity;

      // We assume that the angular velocity is directed along the z-axis, i.e.,
      // we model the actual angular velocity as if it was a two-dimensional
      // rotation, irrespective of the actual value of `spacedim`.
      const double omega = angular_velocity.value(p);
      velocity[0]        = -omega * p[1];
      velocity[1]        = omega * p[0];

      return velocity[component];
    }

  private:
    const Functions::ParsedFunction<spacedim> &angular_velocity;
  };


  // Similarly, we assume that the solid position can be computed explicitly at
  // each time step, exploiting the knowledge of the angular velocity. We
  // compute the exact position of the solid particle assuming that the solid is
  // rotated by an amount equal to the time step multiplied by the angular
  // velocity computed at the point `p`:
  template <int spacedim>
  class SolidPosition : public Function<spacedim>
  {
  public:
    static_assert(spacedim > 1,
                  "Cannot instantiate SolidPosition for spacedim == 1");

    SolidPosition(const Functions::ParsedFunction<spacedim> &angular_velocity,
                  const double                               time_step)
      : Function<spacedim>(spacedim)
      , angular_velocity(angular_velocity)
      , time_step(time_step)
    {}

    virtual double value(const Point<spacedim> &p,
                         unsigned int           component = 0) const override
    {
      Point<spacedim> new_position = p;

      double dtheta = angular_velocity.value(p) * time_step;

      new_position[0] = std::cos(dtheta) * p[0] - std::sin(dtheta) * p[1];
      new_position[1] = std::sin(dtheta) * p[0] + std::cos(dtheta) * p[1];

      return new_position[component];
    }

    void set_time_step(const double new_time_step)
    {
      time_step = new_time_step;
    }

  private:
    const Functions::ParsedFunction<spacedim> &angular_velocity;
    double                                     time_step;
  };


  // @sect3{The StokesImmersedProblem class declaration}

  // We are now ready to introduce the main class of our tutorial program. As
  // usual, other than the constructor, we leave a single public entry point:
  // the `run()` method. Everything else is left `private`, and accessed through
  // the run method itself.
  template <int dim, int spacedim = dim>
  class StokesImmersedProblem
  {
  public:
    StokesImmersedProblem(
      const StokesImmersedProblemParameters<dim, spacedim> &par);

    void run();

    // The next section contains the `private` members of the class.
    // The first method is similar to what is present in previous example.
    // However it not only takes care of generating the grid for the fluid, but
    // also the grid for the solid. The second computes the largest time step
    // that guarantees that each particle moves of at most one cell. This is
    // important to ensure that the Particles::ParticleHandler can find which
    // cell a particle ends up in, as it can only look from one cell to its
    // immediate neighbors (because, in a parallel setting, every MPI process
    // only knows about the cells it owns as well as their immediate neighbors).
  private:
    void make_grid();

    double compute_time_step() const;

    // The next two functions initialize the
    // Particles::ParticleHandler objects used in this class. We have two such
    // objects: One represents passive tracers, used to plot the trajectories
    // of fluid particles, while the the other represents material particles
    // of the solid, which are placed at quadrature points of the solid grid.
    void setup_tracer_particles();
    void setup_solid_particles();

    // The remainder of the set up is split in two parts: The first of the
    // following two functions creates all objects that are needed once per
    // simulation, whereas the other sets up all objects that need to be
    // reinitialized at every refinement step.
    void initial_setup();
    void setup_dofs();

    // The assembly routine is very similar to other Stokes assembly routines,
    // with the exception of the Nitsche restriction part, which exploits one of
    // the particle handlers to integrate on a non-matching part of the fluid
    // domain, corresponding to the position of the solid. We split these two
    // parts into two separate functions.
    void assemble_stokes_system();
    void assemble_nitsche_restriction();

    // The remaining functions solve the linear system (which looks almost
    // identical to the one in step-60) and then postprocess the solution: The
    // refine_and_transfer() method is called only every `refinement_frequency`
    // steps to adapt the mesh and also make sure that all the fields that were
    // computed on the time step before refinement are transferred correctly to
    // the new grid. This includes vector fields, as well as particle
    // information. Similarly, we call the two output methods only every
    // `output_frequency` steps.
    void solve();

    void refine_and_transfer();

    void output_results(const unsigned int cycle, const double time) const;
    void output_particles(const Particles::ParticleHandler<spacedim> &particles,
                          std::string                                 fprefix,
                          const unsigned int                          iter,
                          const double time) const;

    // Let us then move on to the member functions of the class. The first
    // deals with run-time parameters that are read from a parameter file.
    // As noted before, we make sure we cannot modify this object from within
    // this class, by making it a `const` reference.
    const StokesImmersedProblemParameters<dim, spacedim> &par;

    // Then there is also the MPI communicator object that we will use to
    // let processes send information across the network if the program runs
    // in parallel, along with the `pcout` object and timer information
    // that has also been employed by step-40, for example:
    MPI_Comm mpi_communicator;

    ConditionalOStream pcout;

    mutable TimerOutput computing_timer;

    // Next is one of the main novelties with regard to step-60. Here we
    // assume that both the solid and the fluid are fully distributed
    // triangulations. This allows the problem to scale to a very large number
    // of degrees of freedom, at the cost of communicating all the overlapping
    // regions between non matching triangulations. This is especially tricky,
    // since we make no assumptions on the relative position or distribution of
    // the various subdomains of the two triangulations. In particular, we
    // assume that every process owns only a part of the `solid_tria`, and only
    // a part of the `fluid_tria`, not necessarily in the same physical region,
    // and not necessarily overlapping.
    //
    // We could in principle try to create the initial subdivisions in such a
    // way that each process's subdomains overlap between the solid and the
    // fluid regions. However, this overlap would be destroyed during the
    // simulation, and we would have to redistribute the DoFs again and again.
    // The approach we follow in this tutorial is more flexible, and not much
    // more expensive. We make two all-to-all communications at the beginning of
    // the simulation to exchange information about an (approximate) information
    // of the geometrical occupancy of each processor (done through a collection
    // of bounding boxes).
    //
    // This information is used by the Particles::ParticleHandler class
    // to exchange (using a some-to-some communication pattern) all particles,
    // so that every process knows about the particles that live on the
    // region occupied by the fluid subdomain that it owns.
    //
    // In order to couple the overlapping regions, we exploit the facilities
    // implemented in the ParticleHandler class.
    parallel::distributed::Triangulation<spacedim>      fluid_tria;
    parallel::distributed::Triangulation<dim, spacedim> solid_tria;

    // Next come descriptions of the finite elements in use, along with
    // appropriate quadrature formulas and the corresponding DoFHandler objects.
    // For the current implementation, only `fluid_fe` is really necessary. For
    // completeness, and to allow easy extension, we also keep the `solid_fe`
    // around, which is however initialized to a FE_Nothing finite element
    // space, i.e., one that has no degrees of freedom.
    //
    // We declare both finite element spaces as `std::unique_ptr` objects rather
    // than regular member variables, to allow their generation after
    // `StokesImmersedProblemParameters` has been initialized. In particular,
    // they will be initialized in the `initial_setup()` method.
    std::unique_ptr<FiniteElement<spacedim>>      fluid_fe;
    std::unique_ptr<FiniteElement<dim, spacedim>> solid_fe;

    std::unique_ptr<Quadrature<spacedim>> fluid_quadrature_formula;
    std::unique_ptr<Quadrature<dim>>      solid_quadrature_formula;

    DoFHandler<spacedim>      fluid_dh;
    DoFHandler<dim, spacedim> solid_dh;

    std::unique_ptr<MappingFEField<dim, spacedim>> solid_mapping;

    // Similarly to how things are done in step-22, we use a block system to
    // treat the Stokes part of the problem, and follow very closely what was
    // done there.
    std::vector<IndexSet> fluid_owned_dofs;
    std::vector<IndexSet> solid_owned_dofs;

    std::vector<IndexSet> fluid_relevant_dofs;
    std::vector<IndexSet> solid_relevant_dofs;

    // Using this partitioning of degrees of freedom, we can then define all of
    // the objects necessary to describe the linear systems in question:
    AffineConstraints<double> constraints;

    LA::MPI::BlockSparseMatrix system_matrix;
    LA::MPI::BlockSparseMatrix preconditioner_matrix;

    LA::MPI::BlockVector solution;
    LA::MPI::BlockVector locally_relevant_solution;
    LA::MPI::BlockVector system_rhs;

    // Let us move to the particles side of this program. There are two
    // Particles::ParticleHandler objects used to couple the solid with the
    // fluid, and to describe the passive tracers. These, in many ways, play a
    // role similar to the DoFHandler class used in the discretization, i.e.,
    // they provide for an enumeration of particles and allow querying
    // information about each particle.
    Particles::ParticleHandler<spacedim> tracer_particle_handler;
    Particles::ParticleHandler<spacedim> solid_particle_handler;

    // For every tracer particle, we need to compute the velocity field in its
    // current position, and update its position using a discrete time stepping
    // scheme. We do this using distributed linear algebra objects that store
    // the coordinates of each particle's location or velocity. That is, these
    // vectors have `tracer_particle_handler.n_global_particles() * spacedim`
    // entries that we will store in a way so that parts of the vector are
    // partitioned across all processes. (Implicitly, we here make the
    // assumption that the `spacedim` coordinates of each particle are stored in
    // consecutive entries of the vector.) Thus, we need to determine who the
    // owner of each vector entry is. We set this owner to be equal to the
    // process that generated that particle at time $t=0$. This information is
    // stored for every process in the
    // `locally_owned_tracer_particle_coordinates` IndexSet.
    //
    // Once the particles have been distributed around to match the process that
    // owns the region where the particle lives, we will need read access from
    // that process to the corresponding velocity field. We achieve this by
    // filling a read only velocity vector field that contains the relevant
    // information in ghost entries. This is achieved using the
    // `locally_relevant_tracer_particle_coordinates` IndexSet, that keeps track
    // of how things change during the simulation, i.e., it keeps track of where
    // particles that the current process owns have ended up being, and who owns
    // the particles that ended up in my subdomain.
    //
    // While this is not the most efficient strategy, we keep it this way to
    // illustrate how things would work in a real fluid-structure
    // interaction (FSI) problem. If a particle is linked to a specific solid
    // degree of freedom, we are not free to choose who owns it, and we have to
    // communicate this information around. We illustrate this here, and show
    // that the communication pattern is point-to-point, and negligible in terms
    // of total cost of the algorithm.
    //
    // The vectors defined based on these subdivisions are then used to store
    // the particles velocities (read-only, with ghost entries) and their
    // displacement (read/write, no ghost entries).
    IndexSet locally_owned_tracer_particle_coordinates;
    IndexSet locally_relevant_tracer_particle_coordinates;

    LA::MPI::Vector tracer_particle_velocities;
    LA::MPI::Vector relevant_tracer_particle_displacements;

    // One of the key points of this tutorial program is the coupling between
    // two independent parallel::distributed::Triangulation objects, one of
    // which may be moving and deforming (with possibly large deformations) with
    // respect to the other. When both the fluid and the solid triangulations
    // are of type parallel::distributed::Triangulation, every process has
    // access only to its fraction of locally owned cells of each of the two
    // triangulations. As mentioned above, in general, the locally owned domains
    // are not overlapping.
    //
    // In order to allow for the efficient exchange of information between
    // non-overlapping parallel::distributed::Triangulation objects, some
    // algorithms of the library require the user to provide a rough description
    // of the area occupied by the locally owned part of the triangulation, in
    // the form of a collection of axis-aligned bounding boxes for each process,
    // that provide a full covering of the locally owned part of the domain.
    // This kind of information can then be used in situations where one needs
    // to send information to the owner of the cell surrounding a known
    // location, without knowing who that owner may in fact be. But, if one
    // knows a collection of bounding boxes for the geometric area or volume
    // each process owns, then we can determine a subset of all processes that
    // might possibly own the cell in which that location lies: namely, all of
    // those processes whose bounding boxes contain that point. Instead of
    // sending the information associated to that location to all processes, one
    // can then get away with only sending it to a small subset of the processes
    // with point-to-point communication primitives. (You will notice that this
    // also allows for the typical time-vs-memory trade-off: The more data we
    // are willing to store about each process's owned area -- in the form of
    // more refined bounding box information -- the less communication we have
    // to perform.)
    //
    // We construct this information by gathering a vector (of length
    // Utilities::MPI::n_mpi_processes()) of vectors of BoundingBox objects.
    // We fill this vector using the extract_rtree_level() function, and allow
    // the user to select what level of the tree to extract. The "level"
    // corresponds to how coarse/fine the overlap of the area with bounding
    // boxes should be.
    //
    // As an example, this is what would be extracted by the
    // extract_rtree_level() function applied to a two dimensional hyper ball,
    // distributed over three processes. Each image shows in green the bounding
    // boxes associated to the locally owned cells of the triangulation on each
    // process, and in violet the bounding boxes extracted from the rtree:
    //
    // @image html rtree-process-0.png
    // @image html rtree-process-1.png
    // @image html rtree-process-2.png
    //
    // We store these boxes in a global member variable, which is updated at
    // every refinement step:
    std::vector<std::vector<BoundingBox<spacedim>>> global_fluid_bounding_boxes;
  };



  // @sect3{The StokesImmersedProblem class implementation}

  // @sect4{Object construction and mesh initialization functions}

  // In the constructor, we create the mpi_communicator as well as
  // the triangulations and dof_handler for both the fluid and the solid.
  // Using the mpi_communicator, both the ConditionalOStream and TimerOutput
  // object are constructed.
  template <int dim, int spacedim>
  StokesImmersedProblem<dim, spacedim>::StokesImmersedProblem(
    const StokesImmersedProblemParameters<dim, spacedim> &par)
    : par(par)
    , mpi_communicator(MPI_COMM_WORLD)
    , pcout(std::cout,
            (Utilities::MPI::this_mpi_process(mpi_communicator) == 0))
    , computing_timer(mpi_communicator,
                      pcout,
                      TimerOutput::summary,
                      TimerOutput::wall_times)
    , fluid_tria(mpi_communicator,
                 typename Triangulation<spacedim>::MeshSmoothing(
                   Triangulation<spacedim>::smoothing_on_refinement |
                   Triangulation<spacedim>::smoothing_on_coarsening))
    , solid_tria(mpi_communicator,
                 typename Triangulation<dim, spacedim>::MeshSmoothing(
                   Triangulation<dim, spacedim>::smoothing_on_refinement |
                   Triangulation<dim, spacedim>::smoothing_on_coarsening))
    , fluid_dh(fluid_tria)
    , solid_dh(solid_tria)
  {}


  // In order to generate the grid, we first try to use the functions in the
  // deal.II GridGenerator namespace, by leveraging the
  // GridGenerator::generate_from_name_and_argument(). If this function fails,
  // then we use the following method, where the name is interpreted as a
  // filename, and the arguments are interpreted as a map from manifold ids to
  // CAD files, and are converted to Manifold descriptors using the OpenCASCADE
  // namespace facilities. At the top, we read the file into a triangulation:
  template <int dim, int spacedim>
  void read_grid_and_cad_files(const std::string &grid_file_name,
                               const std::string &ids_and_cad_file_names,
                               Triangulation<dim, spacedim> &tria)
  {
    GridIn<dim, spacedim> grid_in;
    grid_in.attach_triangulation(tria);
    grid_in.read(grid_file_name);

    // If we got to this point, then the Triangulation has been read, and we are
    // ready to attach to it the correct manifold descriptions. We perform the
    // next lines of code only if deal.II has been built with OpenCASCADE
    // support. For each entry in the map, we try to open the corresponding CAD
    // file, we analyze it, and according to its content, opt for either a
    // OpenCASCADE::ArcLengthProjectionLineManifold (if the CAD file contains a
    // single `TopoDS_Edge` or a single `TopoDS_Wire`) or a
    // OpenCASCADE::NURBSPatchManifold, if the file contains a single face.
    // Notice that if the CAD files do not contain single wires, edges, or
    // faces, an assertion will be throw in the generation of the Manifold.
    //
    // We use the Patterns::Tools::Convert class to do the conversion from the
    // string to a map between manifold ids and file names for us:
#ifdef DEAL_II_WITH_OPENCASCADE
    using map_type  = std::map<types::manifold_id, std::string>;
    using Converter = Patterns::Tools::Convert<map_type>;

    for (const auto &pair : Converter::to_value(ids_and_cad_file_names))
      {
        const auto &manifold_id   = pair.first;
        const auto &cad_file_name = pair.second;

        const auto extension = boost::algorithm::to_lower_copy(
          cad_file_name.substr(cad_file_name.find_last_of('.') + 1));

        TopoDS_Shape shape;
        if (extension == "iges" || extension == "igs")
          shape = OpenCASCADE::read_IGES(cad_file_name);
        else if (extension == "step" || extension == "stp")
          shape = OpenCASCADE::read_STEP(cad_file_name);
        else
          AssertThrow(false,
                      ExcNotImplemented("We found an extension that we "
                                        "do not recognize as a CAD file "
                                        "extension. Bailing out."));

        // Now we check how many faces are contained in the `Shape`. OpenCASCADE
        // is intrinsically 3D, so if this number is zero, we interpret this as
        // a line manifold, otherwise as a
        // OpenCASCADE::NormalToMeshProjectionManifold in `spacedim` = 3, or
        // OpenCASCADE::NURBSPatchManifold in `spacedim` = 2.
        const auto n_elements = OpenCASCADE::count_elements(shape);
        if ((std::get<0>(n_elements) == 0))
          tria.set_manifold(
            manifold_id,
            OpenCASCADE::ArclengthProjectionLineManifold<dim, spacedim>(shape));
        else if (spacedim == 3)
          {
            // We use this trick, because
            // OpenCASCADE::NormalToMeshProjectionManifold is only implemented
            // for spacedim = 3. The check above makes sure that things actually
            // work correctly.
            const auto t = reinterpret_cast<Triangulation<dim, 3> *>(&tria);
            t->set_manifold(manifold_id,
                            OpenCASCADE::NormalToMeshProjectionManifold<dim, 3>(
                              shape));
          }
        else
          // We also allow surface descriptions in two dimensional spaces based
          // on single NURBS patches. For this to work, the CAD file must
          // contain a single `TopoDS_Face`.
          tria.set_manifold(manifold_id,
                            OpenCASCADE::NURBSPatchManifold<dim, spacedim>(
                              TopoDS::Face(shape)));
      }
#else
    (void)ids_and_cad_file_names;
    AssertThrow(false, ExcNotImplemented("Generation of the grid failed."));
#endif
  }



  // Now let's put things together, and make all the necessary grids. As
  // mentioned above, we first try to generate the grid internally, and if we
  // fail (i.e., if we end up in the `catch` clause), then we proceed with the
  // above function.
  //
  // We repeat this pattern for both the fluid and the solid mesh.
  template <int dim, int spacedim>
  void StokesImmersedProblem<dim, spacedim>::make_grid()
  {
    try
      {
        GridGenerator::generate_from_name_and_arguments(
          fluid_tria, par.name_of_fluid_grid, par.arguments_for_fluid_grid);
      }
    catch (...)
      {
        pcout << "Generating from name and argument failed." << std::endl
              << "Trying to read from file name." << std::endl;
        read_grid_and_cad_files(par.name_of_fluid_grid,
                                par.arguments_for_fluid_grid,
                                fluid_tria);
      }
    fluid_tria.refine_global(par.initial_fluid_refinement);

    try
      {
        GridGenerator::generate_from_name_and_arguments(
          solid_tria, par.name_of_solid_grid, par.arguments_for_solid_grid);
      }
    catch (...)
      {
        read_grid_and_cad_files(par.name_of_solid_grid,
                                par.arguments_for_solid_grid,
                                solid_tria);
      }

    solid_tria.refine_global(par.initial_solid_refinement);
  }

  // @sect4{Particle initialization functions}

  // Once the solid and fluid grids have been created, we start filling the
  // Particles::ParticleHandler objects. The first one we take care of is the
  // one we use to keep track of passive tracers in the fluid. These are
  // simply transported along, and in some sense their locations are
  // unimportant: We just want to use them to see where flow is being
  // transported. We could use any way we choose to determine where they are
  // initially located. A convenient one is to create the initial locations as
  // the vertices of a mesh in a shape of our choice -- a choice determined by
  // one of the run-time parameters in the parameter file.
  //
  // In this implementation, we create tracers using the support points of a
  // FE_Q finite element space defined on a temporary grid, which is then
  // discarded. Of this grid, we only keep around the Particles::Particle
  // objects (stored in a Particles::ParticleHandler class) associated to the
  // support points.
  //
  // The Particles::ParticleHandler class offers the possibility to insert a set
  // of particles that live physically in the part of the domain owned by the
  // active process. However, in this case this function would not suffice. The
  // particles generated as the locally owned support points of an FE_Q object
  // on an arbitrary grid (non-matching with regard to the fluid grid) have no
  // reasons to lie in the same physical region of the locally owned subdomain
  // of the fluid grid. In fact this will almost never be the case, especially
  // since we want to keep track of what is happening to the particles
  // themselves.
  //
  // In particle-in-cell methods (PIC), it is often customary to assign
  // ownership of the particles to the process where the particles lie. In this
  // tutorial we illustrate a different approach, which is useful if one wants
  // to keep track of information related to the particles (for example, if a
  // particle is associated to a given degree of freedom, which is owned by a
  // specific process and not necessarily the same process that owns the fluid
  // cell where the particle happens to be at any given time).
  // In the approach used here, ownership of the particles is assigned once at
  // the beginning, and one-to-one communication happens whenever the original
  // owner needs information from the process that owns the cell where the
  // particle lives. We make sure that we set ownership of the particles using
  // the initial particle distribution, and keep the same ownership throughout
  // the execution of the program.
  //
  // With this overview out of the way, let us see what the function does. At
  // the top, we create a temporary triangulation and DoFHandler object from
  // which we will take the node locations for initial particle locations:
  template <int dim, int spacedim>
  void StokesImmersedProblem<dim, spacedim>::setup_tracer_particles()
  {
    parallel::distributed::Triangulation<spacedim> particle_insert_tria(
      mpi_communicator);
    GridGenerator::generate_from_name_and_arguments(
      particle_insert_tria,
      par.name_of_particle_grid,
      par.arguments_for_particle_grid);
    particle_insert_tria.refine_global(par.particle_insertion_refinement);

    FE_Q<spacedim>       particles_fe(1);
    DoFHandler<spacedim> particles_dof_handler(particle_insert_tria);
    particles_dof_handler.distribute_dofs(particles_fe);

    // This is where things start to get complicated. Since we may run
    // this program in a parallel environment, every parallel process will now
    // have created these temporary triangulations and DoFHandlers. But, in
    // fully distributed triangulations, the active process only knows about the
    // locally owned cells, and has no idea of how other processes have
    // distributed their own cells. This is true for both the temporary
    // triangulation created above as well as the fluid triangulation into which
    // we want to embed the particles below. On the other hand, these locally
    // known portions of the two triangulations will, in general, not overlap.
    // That is, the locations of the particles we will create from the node
    // locations of the temporary mesh are arbitrary, and may fall within a
    // region of the fluid triangulation that the current process doesn't have
    // access to (i.e., a region of the fluid domain where cells are
    // artificial). In order to understand who to send those particles to, we
    // need to have a (rough) idea of how the fluid grid is distributed among
    // processors.
    //
    // We construct this information by first building an index tree of boxes
    // bounding the locally owned cells, and then extracting one of the first
    // levels of the tree:
    std::vector<BoundingBox<spacedim>> all_boxes;
    all_boxes.reserve(fluid_tria.n_locally_owned_active_cells());
    for (const auto &cell : fluid_tria.active_cell_iterators())
      if (cell->is_locally_owned())
        all_boxes.emplace_back(cell->bounding_box());

    const auto tree = pack_rtree(all_boxes);
    const auto local_boxes =
      extract_rtree_level(tree, par.fluid_rtree_extraction_level);

    // Each process now has a collection of bounding boxes that completely
    // enclose all locally owned processes (but that may overlap the bounding
    // boxes of other processes). We then exchange this information between all
    // participating processes so that every process knows the bounding boxes of
    // all other processes.
    //
    // Equipped with this knowledge, we can then initialize the
    // `tracer_particle_handler` to the fluid mesh and generate the particles
    // from the support points of the (temporary) tracer particles
    // triangulation. This function call uses the `global_bounding_boxes` object
    // we just constructed to figure out where to send the particles whose
    // locations were derived from the locally owned part of the
    // `particles_dof_handler`. At the end of this call, every particle will
    // have been distributed to the correct process (i.e., the process that owns
    // the fluid cell where the particle lives). We also output their number to
    // the screen at this point.
    global_fluid_bounding_boxes =
      Utilities::MPI::all_gather(mpi_communicator, local_boxes);

    tracer_particle_handler.initialize(fluid_tria,
                                       StaticMappingQ1<spacedim>::mapping);

    Particles::Generators::dof_support_points(particles_dof_handler,
                                              global_fluid_bounding_boxes,
                                              tracer_particle_handler);

    pcout << "Tracer particles: "
          << tracer_particle_handler.n_global_particles() << std::endl;

    // Each particle so created has a unique ID. At some point in the
    // algorithm below, we will need vectors containing position and velocity
    // information for each particle. This vector will have size `n_particles *
    // spacedim`, and we will have to store the elements of this vector in a way
    // so that each parallel process "owns" those elements that correspond to
    // coordinates of the particles it owns. In other words, we have to
    // partition the index space between zero and `n_particles * spacedim` among
    // all processes. We can do this by querying the `tracer_particle_handler`
    // for the IDs of its locally relevant particles, and construct the indices
    // that would be needed to store in a (parallel distributed) vector of the
    // position and velocity of all particles where we implicitly assume that we
    // store the coordinates of each location or velocity in `spacedim`
    // successive vector elements (this is what the IndexSet::tensor_priduct()
    // function does).
    locally_owned_tracer_particle_coordinates =
      tracer_particle_handler.locally_owned_particle_ids().tensor_product(
        complete_index_set(spacedim));

    // At the beginning of the simulation, all particles are in their original
    // position. When particles move, they may traverse to a part of the domain
    // which is owned by another process. If this happens, the current process
    // keeps formally "ownership" of the particles, but may need read access
    // from the process where the particle has landed. We keep this information
    // in another index set, which stores the indices of all particles that are
    // currently on the current process's subdomain, independently if they have
    // always been here or not.
    //
    // Keeping this index set around allows us to leverage linear algebra
    // classes for all communications regarding positions and velocities of the
    // particles. This mimics what would happen in the case where another
    // problem was solved in the solid domain (as in fluid-structure
    // interaction. In this latter case, additional DOFs on the solid domain
    // would be coupled to what is occurring in the fluid domain.
    locally_relevant_tracer_particle_coordinates =
      locally_owned_tracer_particle_coordinates;

    // Finally, we make sure that upon refinement, particles are correctly
    // transferred. When performing local refinement or coarsening, particles
    // will land in another cell. We could in principle redistribute all
    // particles after refining, however this would be overly expensive.
    //
    // The Particles::ParticleHandler class has a way to transfer information
    // from a cell to its children or to its parent upon refinement, without the
    // need to reconstruct the entire data structure. This is done by
    // registering two callback functions to the triangulation. These
    // functions will receive a signal when refinement is about to happen, and
    // when it has just happened, and will take care of transferring all
    // information to the newly refined grid with minimal computational cost.
    fluid_tria.signals.pre_distributed_refinement.connect(
      [&]() { tracer_particle_handler.register_store_callback_function(); });

    fluid_tria.signals.post_distributed_refinement.connect([&]() {
      tracer_particle_handler.register_load_callback_function(false);
    });
  }


  // Similarly to what we have done for passive tracers, we next set up the
  // particles that track the quadrature points of the solid mesh. The main
  // difference here is that we also want to attach a weight value (the "JxW"
  // value of the quadrature point) to each of particle, so that we can compute
  // integrals even without direct access to the original solid grid.
  //
  // This is achieved by leveraging the "properties" concept of the
  // Particles::Particle class. It is possible to store (in a memory
  // efficient way) an arbitrary number of `double` numbers for each of the
  // Particles::Particle objects inside a Particles::ParticleHandler object. We
  // use this possibility to store the JxW values of the quadrature points of
  // the solid grid.
  //
  // In our case, we only need to store one property per particle: the JxW value
  // of the integration on the solid grid. This is passed at construction time
  // to the solid_particle_handler object as the last argument
  template <int dim, int spacedim>
  void StokesImmersedProblem<dim, spacedim>::setup_solid_particles()
  {
    QGauss<dim> quadrature(fluid_fe->degree + 1);

    const unsigned int n_properties = 1;
    solid_particle_handler.initialize(fluid_tria,
                                      StaticMappingQ1<spacedim>::mapping,
                                      n_properties);

    // The number of particles that we generate locally is equal to the total
    // number of locally owned cells times the number of quadrature points used
    // in each cell. We store all these points in a vector, and their
    // corresponding properties in a vector of vectors:
    std::vector<Point<spacedim>> quadrature_points_vec;
    quadrature_points_vec.reserve(quadrature.size() *
                                  solid_tria.n_locally_owned_active_cells());

    std::vector<std::vector<double>> properties;
    properties.reserve(quadrature.size() *
                       solid_tria.n_locally_owned_active_cells());

    FEValues<dim, spacedim> fe_v(*solid_fe,
                                 quadrature,
                                 update_JxW_values | update_quadrature_points);
    for (const auto &cell : solid_dh.active_cell_iterators())
      if (cell->is_locally_owned())
        {
          fe_v.reinit(cell);
          const auto &points = fe_v.get_quadrature_points();
          const auto &JxW    = fe_v.get_JxW_values();

          for (unsigned int q = 0; q < points.size(); ++q)
            {
              quadrature_points_vec.emplace_back(points[q]);
              properties.emplace_back(
                std::vector<double>(n_properties, JxW[q]));
            }
        }

    // We proceed in the same way we did with the tracer particles, reusing the
    // computed bounding boxes. However, we first check that the
    // `global_fluid_bounding_boxes` object has been actually filled. This
    // should certainly be the case here, since this method is called after the
    // one that initializes the tracer particles. However, we want to make sure
    // that if in the future someone decides (for whatever reason) to initialize
    // first the solid particle handler, or to copy just this part of the
    // tutorial, a meaningful exception is thrown when things don't work as
    // expected
    //
    // Since we have already stored the position of the quadrature points,
    // we can use these positions to insert the particles directly using
    // the `solid_particle_handler` instead of having to go through a
    // Particles::Generators function:
    Assert(!global_fluid_bounding_boxes.empty(),
           ExcInternalError(
             "I was expecting the "
             "global_fluid_bounding_boxes to be filled at this stage. "
             "Make sure you fill this vector before trying to use it "
             "here. Bailing out."));

    solid_particle_handler.insert_global_particles(quadrature_points_vec,
                                                   global_fluid_bounding_boxes,
                                                   properties);


    // As in the previous function, we end by making sure that upon refinement,
    // particles are correctly transferred:
    fluid_tria.signals.pre_distributed_refinement.connect(
      [&]() { solid_particle_handler.register_store_callback_function(); });

    fluid_tria.signals.post_distributed_refinement.connect(
      [&]() { solid_particle_handler.register_load_callback_function(false); });

    pcout << "Solid particles: " << solid_particle_handler.n_global_particles()
          << std::endl;
  }



  // @sect4{DoF initialization functions}

  // We set up the finite element space and the quadrature formula to be
  // used throughout the step. For the fluid, we use Taylor-Hood elements (e.g.
  // $Q_k \times Q_{k-1}$). Since we do not solve any equation on the solid
  // domain, an empty finite element space is generated. A natural extension of
  // this program would be to solve a fluid structure interaction problem, which
  // would require that the `solid_fe` use more useful FiniteElement class.
  //
  // Like for many other functions, we store the time necessary to carry out the
  // operations we perform here. The current function puts its timing
  // information into a section with label "Initial setup". Numerous other calls
  // to this timer are made in various functions. They allow to monitor the
  // absolute and relative cost of each individual function to identify
  // bottlenecks.
  template <int dim, int spacedim>
  void StokesImmersedProblem<dim, spacedim>::initial_setup()
  {
    TimerOutput::Scope t(computing_timer, "Initial setup");

    fluid_fe =
      std::make_unique<FESystem<spacedim>>(FE_Q<spacedim>(par.velocity_degree),
                                           spacedim,
                                           FE_Q<spacedim>(par.velocity_degree -
                                                          1),
                                           1);


    solid_fe = std::make_unique<FE_Nothing<dim, spacedim>>();
    solid_dh.distribute_dofs(*solid_fe);

    fluid_quadrature_formula =
      std::make_unique<QGauss<spacedim>>(par.velocity_degree + 1);
    solid_quadrature_formula =
      std::make_unique<QGauss<dim>>(par.velocity_degree + 1);
  }


  // We next construct the distributed block matrices and vectors which are used
  // to solve the linear equations that arise from the problem. This function is
  // adapted from step-55 and we refer to this step for a thorough explanation.
  template <int dim, int spacedim>
  void StokesImmersedProblem<dim, spacedim>::setup_dofs()
  {
    TimerOutput::Scope t(computing_timer, "Setup dofs");

    fluid_dh.distribute_dofs(*fluid_fe);

    std::vector<unsigned int> stokes_sub_blocks(spacedim + 1, 0);
    stokes_sub_blocks[spacedim] = 1;
    DoFRenumbering::component_wise(fluid_dh, stokes_sub_blocks);

    auto dofs_per_block =
      DoFTools::count_dofs_per_fe_block(fluid_dh, stokes_sub_blocks);

    const unsigned int n_u = dofs_per_block[0], n_p = dofs_per_block[1];

    pcout << "   Number of degrees of freedom: " << fluid_dh.n_dofs() << " ("
          << n_u << '+' << n_p << " -- "
          << solid_particle_handler.n_global_particles() << '+'
          << tracer_particle_handler.n_global_particles() << ')' << std::endl;

    fluid_owned_dofs.resize(2);
    fluid_owned_dofs[0] = fluid_dh.locally_owned_dofs().get_view(0, n_u);
    fluid_owned_dofs[1] =
      fluid_dh.locally_owned_dofs().get_view(n_u, n_u + n_p);

    IndexSet locally_relevant_dofs;
    DoFTools::extract_locally_relevant_dofs(fluid_dh, locally_relevant_dofs);
    fluid_relevant_dofs.resize(2);
    fluid_relevant_dofs[0] = locally_relevant_dofs.get_view(0, n_u);
    fluid_relevant_dofs[1] = locally_relevant_dofs.get_view(n_u, n_u + n_p);

    {
      constraints.reinit(locally_relevant_dofs);

      FEValuesExtractors::Vector velocities(0);
      DoFTools::make_hanging_node_constraints(fluid_dh, constraints);
      VectorTools::interpolate_boundary_values(
        fluid_dh,
        0,
        Functions::ZeroFunction<spacedim>(spacedim + 1),
        constraints,
        fluid_fe->component_mask(velocities));
      constraints.close();
    }

    auto locally_owned_dofs_per_processor =
      Utilities::MPI::all_gather(mpi_communicator,
                                 fluid_dh.locally_owned_dofs());
    {
      system_matrix.clear();

      Table<2, DoFTools::Coupling> coupling(spacedim + 1, spacedim + 1);
      for (unsigned int c = 0; c < spacedim + 1; ++c)
        for (unsigned int d = 0; d < spacedim + 1; ++d)
          if (c == spacedim && d == spacedim)
            coupling[c][d] = DoFTools::none;
          else if (c == spacedim || d == spacedim || c == d)
            coupling[c][d] = DoFTools::always;
          else
            coupling[c][d] = DoFTools::none;

      BlockDynamicSparsityPattern dsp(dofs_per_block, dofs_per_block);

      DoFTools::make_sparsity_pattern(
        fluid_dh, coupling, dsp, constraints, false);

      SparsityTools::distribute_sparsity_pattern(
        dsp,
        locally_owned_dofs_per_processor,
        mpi_communicator,
        locally_relevant_dofs);

      system_matrix.reinit(fluid_owned_dofs, dsp, mpi_communicator);
    }

    {
      preconditioner_matrix.clear();

      Table<2, DoFTools::Coupling> coupling(spacedim + 1, spacedim + 1);
      for (unsigned int c = 0; c < spacedim + 1; ++c)
        for (unsigned int d = 0; d < spacedim + 1; ++d)
          if (c == spacedim && d == spacedim)
            coupling[c][d] = DoFTools::always;
          else
            coupling[c][d] = DoFTools::none;

      BlockDynamicSparsityPattern dsp(dofs_per_block, dofs_per_block);

      DoFTools::make_sparsity_pattern(
        fluid_dh, coupling, dsp, constraints, false);
      SparsityTools::distribute_sparsity_pattern(
        dsp,
        locally_owned_dofs_per_processor,
        mpi_communicator,
        locally_relevant_dofs);
      preconditioner_matrix.reinit(fluid_owned_dofs, dsp, mpi_communicator);
    }

    locally_relevant_solution.reinit(fluid_owned_dofs,
                                     fluid_relevant_dofs,
                                     mpi_communicator);
    system_rhs.reinit(fluid_owned_dofs, mpi_communicator);
    solution.reinit(fluid_owned_dofs, mpi_communicator);
  }


  // @sect4{Assembly functions}

  // We assemble the system matrix, the preconditioner matrix, and the right
  // hand side. The code is adapted from step-55, which is essentially what
  // step-27 also has, and is pretty standard if you know what the Stokes
  // equations look like.
  template <int dim, int spacedim>
  void StokesImmersedProblem<dim, spacedim>::assemble_stokes_system()
  {
    system_matrix         = 0;
    preconditioner_matrix = 0;
    system_rhs            = 0;

    TimerOutput::Scope t(computing_timer, "Assemble Stokes terms");

    FEValues<spacedim> fe_values(*fluid_fe,
                                 *fluid_quadrature_formula,
                                 update_values | update_gradients |
                                   update_quadrature_points |
                                   update_JxW_values);

    const unsigned int dofs_per_cell = fluid_fe->n_dofs_per_cell();
    const unsigned int n_q_points    = fluid_quadrature_formula->size();

    FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);
    FullMatrix<double> cell_matrix2(dofs_per_cell, dofs_per_cell);
    Vector<double>     cell_rhs(dofs_per_cell);

    std::vector<Vector<double>> rhs_values(n_q_points,
                                           Vector<double>(spacedim + 1));

    std::vector<Tensor<2, spacedim>> grad_phi_u(dofs_per_cell);
    std::vector<double>              div_phi_u(dofs_per_cell);
    std::vector<double>              phi_p(dofs_per_cell);

    std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);
    const FEValuesExtractors::Vector     velocities(0);
    const FEValuesExtractors::Scalar     pressure(spacedim);

    for (const auto &cell : fluid_dh.active_cell_iterators())
      if (cell->is_locally_owned())
        {
          cell_matrix  = 0;
          cell_matrix2 = 0;
          cell_rhs     = 0;

          fe_values.reinit(cell);
          par.rhs.vector_value_list(fe_values.get_quadrature_points(),
                                    rhs_values);
          for (unsigned int q = 0; q < n_q_points; ++q)
            {
              for (unsigned int k = 0; k < dofs_per_cell; ++k)
                {
                  grad_phi_u[k] = fe_values[velocities].gradient(k, q);
                  div_phi_u[k]  = fe_values[velocities].divergence(k, q);
                  phi_p[k]      = fe_values[pressure].value(k, q);
                }

              for (unsigned int i = 0; i < dofs_per_cell; ++i)
                {
                  for (unsigned int j = 0; j < dofs_per_cell; ++j)
                    {
                      cell_matrix(i, j) +=
                        (par.viscosity *
                           scalar_product(grad_phi_u[i], grad_phi_u[j]) -
                         div_phi_u[i] * phi_p[j] - phi_p[i] * div_phi_u[j]) *
                        fe_values.JxW(q);

                      cell_matrix2(i, j) += 1.0 / par.viscosity * phi_p[i] *
                                            phi_p[j] * fe_values.JxW(q);
                    }

                  const unsigned int component_i =
                    fluid_fe->system_to_component_index(i).first;
                  cell_rhs(i) += fe_values.shape_value(i, q) *
                                 rhs_values[q](component_i) * fe_values.JxW(q);
                }
            }


          cell->get_dof_indices(local_dof_indices);
          constraints.distribute_local_to_global(cell_matrix,
                                                 cell_rhs,
                                                 local_dof_indices,
                                                 system_matrix,
                                                 system_rhs);

          constraints.distribute_local_to_global(cell_matrix2,
                                                 local_dof_indices,
                                                 preconditioner_matrix);
        }

    system_matrix.compress(VectorOperation::add);
    preconditioner_matrix.compress(VectorOperation::add);
    system_rhs.compress(VectorOperation::add);
  }


  // The following method is then the one that deals with the penalty terms that
  // result from imposing the velocity on the impeller. It is, in a sense, the
  // heart of the tutorial, but it is relatively straightforward. Here we
  // exploit the `solid_particle_handler` to compute the Nitsche restriction or
  // the penalization in the embedded domain.
  template <int dim, int spacedim>
  void StokesImmersedProblem<dim, spacedim>::assemble_nitsche_restriction()
  {
    TimerOutput::Scope t(computing_timer, "Assemble Nitsche terms");

    const FEValuesExtractors::Vector velocities(0);
    const FEValuesExtractors::Scalar pressure(spacedim);

    SolidVelocity<spacedim> solid_velocity(par.angular_velocity);

    std::vector<types::global_dof_index> fluid_dof_indices(
      fluid_fe->n_dofs_per_cell());

    FullMatrix<double>     local_matrix(fluid_fe->n_dofs_per_cell(),
                                    fluid_fe->n_dofs_per_cell());
    dealii::Vector<double> local_rhs(fluid_fe->n_dofs_per_cell());

    const auto penalty_parameter =
      1.0 / GridTools::minimal_cell_diameter(fluid_tria);

    // We loop over all the local particles. Although this could be achieved
    // directly by looping over all the cells, this would force us
    // to loop over numerous cells which do not contain particles.
    // Consequently, we loop over all the particles, but, we get the reference
    // of the cell in which the particle lies and then loop over all particles
    // within that cell. This enables us to skip the cells which do not contain
    // particles, yet to assemble the local matrix and rhs of each cell to apply
    // the Nitsche restriction. Once we are done with all particles on one cell,
    // we advance the `particle` iterator to the particle past the end of the
    // ones on the current cell (this is the last line of the `while` loop's
    // body).
    auto particle = solid_particle_handler.begin();
    while (particle != solid_particle_handler.end())
      {
        local_matrix = 0;
        local_rhs    = 0;

        // We get an iterator to the cell within which the particle lies from
        // the particle itself. We can then assemble the additional
        // terms in the system matrix and the right hand side as we would
        // normally.
        const auto &cell = particle->get_surrounding_cell(fluid_tria);
        const auto &dh_cell =
          typename DoFHandler<spacedim>::cell_iterator(*cell, &fluid_dh);
        dh_cell->get_dof_indices(fluid_dof_indices);

        // So then let us get the collection of cells that are located on this
        // cell and iterate over them. From each particle we gather the location
        // and the reference location of the particle as well as the additional
        // information that is attached to the particle. In the present case,
        // this information is the "JxW" of the quadrature points which were
        // used to generate the particles.
        //
        // Using this information, we can add the contribution of the quadrature
        // point to the local_matrix and local_rhs. We can evaluate the value of
        // the shape function at the position of each particle easily by using
        // its reference location.
        const auto pic = solid_particle_handler.particles_in_cell(cell);
        Assert(pic.begin() == particle, ExcInternalError());
        for (const auto &p : pic)
          {
            const auto &ref_q  = p.get_reference_location();
            const auto &real_q = p.get_location();
            const auto &JxW    = p.get_properties()[0];

            for (unsigned int i = 0; i < fluid_fe->n_dofs_per_cell(); ++i)
              {
                const auto comp_i =
                  fluid_fe->system_to_component_index(i).first;
                if (comp_i < spacedim)
                  {
                    for (unsigned int j = 0; j < fluid_fe->n_dofs_per_cell();
                         ++j)
                      {
                        const auto comp_j =
                          fluid_fe->system_to_component_index(j).first;
                        if (comp_i == comp_j)
                          local_matrix(i, j) +=
                            penalty_parameter * par.penalty_term *
                            fluid_fe->shape_value(i, ref_q) *
                            fluid_fe->shape_value(j, ref_q) * JxW;
                      }
                    local_rhs(i) += penalty_parameter * par.penalty_term *
                                    solid_velocity.value(real_q, comp_i) *
                                    fluid_fe->shape_value(i, ref_q) * JxW;
                  }
              }
          }

        constraints.distribute_local_to_global(local_matrix,
                                               local_rhs,
                                               fluid_dof_indices,
                                               system_matrix,
                                               system_rhs);
        particle = pic.end();
      }

    system_matrix.compress(VectorOperation::add);
    system_rhs.compress(VectorOperation::add);
  }


  // @sect4{Solving the linear system}

  // This function solves the linear system with FGMRES with a block diagonal
  // preconditioner and an algebraic multigrid (AMG) method for the diagonal
  // blocks. The preconditioner applies a V cycle to the $(0,0)$ (i.e., the
  // velocity-velocity) block and a CG with the mass matrix for the $(1,1)$
  // block (which is our approximation to the Schur complement: the pressure
  // mass matrix assembled above).
  template <int dim, int spacedim>
  void StokesImmersedProblem<dim, spacedim>::solve()
  {
    TimerOutput::Scope t(computing_timer, "Solve");

    LA::MPI::PreconditionAMG prec_A;
    {
      LA::MPI::PreconditionAMG::AdditionalData data;

#ifdef USE_PETSC_LA
      data.symmetric_operator = true;
#endif
      prec_A.initialize(system_matrix.block(0, 0), data);
    }

    LA::MPI::PreconditionAMG prec_S;
    {
      LA::MPI::PreconditionAMG::AdditionalData data;

#ifdef USE_PETSC_LA
      data.symmetric_operator = true;
#endif
      prec_S.initialize(preconditioner_matrix.block(1, 1), data);
    }

    const auto A = linear_operator<LA::MPI::Vector>(system_matrix.block(0, 0));
    const auto amgA = linear_operator(A, prec_A);

    const auto S =
      linear_operator<LA::MPI::Vector>(preconditioner_matrix.block(1, 1));
    const auto amgS = linear_operator(S, prec_S);

    ReductionControl          inner_solver_control(100,
                                          1e-8 * system_rhs.l2_norm(),
                                          1.e-2);
    SolverCG<LA::MPI::Vector> cg(inner_solver_control);

    const auto invS = inverse_operator(S, cg, amgS);

    const auto P = block_diagonal_operator<2, LA::MPI::BlockVector>(
      std::array<
        dealii::LinearOperator<typename LA::MPI::BlockVector::BlockType>,
        2>{{amgA, amgS}});

    SolverControl solver_control(system_matrix.m(),
                                 1e-10 * system_rhs.l2_norm());

    SolverFGMRES<LA::MPI::BlockVector> solver(solver_control);

    constraints.set_zero(solution);

    solver.solve(system_matrix, solution, system_rhs, P);


    pcout << "   Solved in " << solver_control.last_step() << " iterations."
          << std::endl;

    constraints.distribute(solution);

    locally_relevant_solution = solution;
    const double mean_pressure =
      VectorTools::compute_mean_value(fluid_dh,
                                      QGauss<spacedim>(par.velocity_degree + 2),
                                      locally_relevant_solution,
                                      spacedim);
    solution.block(1).add(-mean_pressure);
    locally_relevant_solution.block(1) = solution.block(1);
  }



  // @sect4{Mesh refinement}

  // We deal with mesh refinement in a completely standard way:
  template <int dim, int spacedim>
  void StokesImmersedProblem<dim, spacedim>::refine_and_transfer()
  {
    TimerOutput::Scope               t(computing_timer, "Refine");
    const FEValuesExtractors::Vector velocity(0);

    Vector<float> error_per_cell(fluid_tria.n_active_cells());
    KellyErrorEstimator<spacedim>::estimate(fluid_dh,
                                            QGauss<spacedim - 1>(
                                              par.velocity_degree + 1),
                                            {},
                                            locally_relevant_solution,
                                            error_per_cell,
                                            fluid_fe->component_mask(velocity));

    if (par.refinement_strategy == "fixed_fraction")
      {
        parallel::distributed::GridRefinement::
          refine_and_coarsen_fixed_fraction(fluid_tria,
                                            error_per_cell,
                                            par.refinement_fraction,
                                            par.coarsening_fraction);
      }
    else if (par.refinement_strategy == "fixed_number")
      {
        parallel::distributed::GridRefinement::refine_and_coarsen_fixed_number(
          fluid_tria,
          error_per_cell,
          par.refinement_fraction,
          par.coarsening_fraction,
          par.max_cells);
      }

    for (const auto &cell : fluid_tria.active_cell_iterators())
      {
        if (cell->refine_flag_set() &&
            cell->level() == par.max_level_refinement)
          cell->clear_refine_flag();
        if (cell->coarsen_flag_set() &&
            cell->level() == par.min_level_refinement)
          cell->clear_coarsen_flag();
      }

    parallel::distributed::SolutionTransfer<spacedim, LA::MPI::BlockVector>
      transfer(fluid_dh);
    fluid_tria.prepare_coarsening_and_refinement();
    transfer.prepare_for_coarsening_and_refinement(locally_relevant_solution);
    fluid_tria.execute_coarsening_and_refinement();

    setup_dofs();

    transfer.interpolate(solution);
    constraints.distribute(solution);
    locally_relevant_solution = solution;
  }


  // @sect4{Creating output for visualization}

  // We output the results (velocity and pressure) on the fluid domain
  // using the standard parallel capabilities of deal.II. A single compressed
  // vtu file is written that agglomerates the information of all processors. An
  // additional `.pvd` record is written to associate the physical time to the
  // vtu files.
  template <int dim, int spacedim>
  void
  StokesImmersedProblem<dim, spacedim>::output_results(const unsigned int cycle,
                                                       double time) const
  {
    TimerOutput::Scope t(computing_timer, "Output fluid");

    std::vector<std::string> solution_names(spacedim, "velocity");
    solution_names.emplace_back("pressure");
    std::vector<DataComponentInterpretation::DataComponentInterpretation>
      data_component_interpretation(
        spacedim, DataComponentInterpretation::component_is_part_of_vector);
    data_component_interpretation.push_back(
      DataComponentInterpretation::component_is_scalar);

    DataOut<spacedim> data_out;
    data_out.attach_dof_handler(fluid_dh);
    data_out.add_data_vector(locally_relevant_solution,
                             solution_names,
                             DataOut<spacedim>::type_dof_data,
                             data_component_interpretation);


    Vector<float> subdomain(fluid_tria.n_active_cells());
    for (unsigned int i = 0; i < subdomain.size(); ++i)
      subdomain(i) = fluid_tria.locally_owned_subdomain();
    data_out.add_data_vector(subdomain, "subdomain");

    data_out.build_patches();

    const std::string filename =
      "solution-" + Utilities::int_to_string(cycle) + ".vtu";
    data_out.write_vtu_in_parallel(par.output_directory + "/" + filename,
                                   mpi_communicator);

    static std::vector<std::pair<double, std::string>> times_and_names;
    times_and_names.push_back(std::make_pair(time, filename));
    std::ofstream ofile(par.output_directory + "/" + "solution.pvd");
    DataOutBase::write_pvd_record(ofile, times_and_names);
  }


  // Similarly, we write the particles (either from the solid or the tracers)
  // as a single compressed vtu file through the Particles::DataOut object.
  // This simple object does not write the additional information
  // attached as "properties" to the particles, but only writes their id -- but
  // then, we don't care about the "JxW" values of these particle locations
  // anyway, so no information that we may have wanted to visualize is lost.
  template <int dim, int spacedim>
  void StokesImmersedProblem<dim, spacedim>::output_particles(
    const Particles::ParticleHandler<spacedim> &particles,
    std::string                                 fprefix,
    const unsigned int                          iter,
    const double                                time) const
  {
    Particles::DataOut<spacedim> particles_out;
    particles_out.build_patches(particles);
    const std::string filename =
      (fprefix + "-" + Utilities::int_to_string(iter) + ".vtu");
    particles_out.write_vtu_in_parallel(par.output_directory + "/" + filename,
                                        mpi_communicator);


    static std::map<std::string, std::vector<std::pair<double, std::string>>>
      times_and_names;
    if (times_and_names.find(fprefix) != times_and_names.end())
      times_and_names[fprefix].push_back(std::make_pair(time, filename));
    else
      times_and_names[fprefix] = {std::make_pair(time, filename)};
    std::ofstream ofile(par.output_directory + "/" + fprefix + ".pvd");
    DataOutBase::write_pvd_record(ofile, times_and_names[fprefix]);
  }


  // @sect4{The "run" function}

  // This function now orchestrates the entire simulation. It is very similar
  // to the other time dependent tutorial programs -- take step-21 or step-26 as
  // an example. At the beginning, we output some status information and also
  // save all current parameters to a file in the output directory, for
  // reproducibility.
  template <int dim, int spacedim>
  void StokesImmersedProblem<dim, spacedim>::run()
  {
#ifdef USE_PETSC_LA
    pcout << "Running StokesImmersedProblem<"
          << Utilities::dim_string(dim, spacedim) << "> using PETSc."
          << std::endl;
#else
    pcout << "Running StokesImmersedProblem<"
          << Utilities::dim_string(dim, spacedim) << "> using Trilinos."
          << std::endl;
#endif
    par.prm.print_parameters(par.output_directory + "/" + "used_parameters_" +
                               std::to_string(dim) + std::to_string(spacedim) +
                               ".prm",
                             ParameterHandler::Short);

    // We then start the time loop. We initialize all the elements of the
    // simulation in the first cycle
    const double time_step    = par.final_time / (par.number_of_time_steps - 1);
    double       time         = 0;
    unsigned int output_cycle = 0;

    for (unsigned int cycle = 0; cycle < par.number_of_time_steps;
         ++cycle, time += time_step)
      {
        par.set_time(time);
        pcout << "Cycle " << cycle << ':' << std::endl
              << "Time : " << time << ", time step: " << time_step << std::endl;

        if (cycle == 0)
          {
            make_grid();
            initial_setup();
            setup_dofs();
            setup_tracer_particles();
            setup_solid_particles();
            tracer_particle_velocities.reinit(
              locally_owned_tracer_particle_coordinates, mpi_communicator);
            output_results(output_cycle, time);
            {
              TimerOutput::Scope t(computing_timer, "Output tracer particles");
              output_particles(tracer_particle_handler,
                               "tracer",
                               output_cycle,
                               time);
            }
            {
              TimerOutput::Scope t(computing_timer, "Output solid particles");
              output_particles(solid_particle_handler,
                               "solid",
                               output_cycle,
                               time);
            }
          }
        // After the first time step, we displace the solid body at the
        // beginning of each time step to take into account the fact that is has
        // moved.
        else
          {
            TimerOutput::Scope t(computing_timer,
                                 "Set solid particle position");

            SolidPosition<spacedim> solid_position(par.angular_velocity,
                                                   time_step);
            solid_particle_handler.set_particle_positions(solid_position,
                                                          false);
          }

        // In order to update the state of the system, we first
        // interpolate the fluid velocity at the position of the tracer
        // particles and, with a naive explicit Euler scheme, advect the
        // massless tracer particles.
        {
          TimerOutput::Scope t(computing_timer, "Set tracer particle motion");
          Particles::Utilities::interpolate_field_on_particles(
            fluid_dh,
            tracer_particle_handler,
            locally_relevant_solution,
            tracer_particle_velocities,
            fluid_fe->component_mask(FEValuesExtractors::Vector(0)));

          tracer_particle_velocities *= time_step;

          locally_relevant_tracer_particle_coordinates =
            tracer_particle_handler.locally_owned_particle_ids().tensor_product(
              complete_index_set(spacedim));

          relevant_tracer_particle_displacements.reinit(
            locally_owned_tracer_particle_coordinates,
            locally_relevant_tracer_particle_coordinates,
            mpi_communicator);

          relevant_tracer_particle_displacements = tracer_particle_velocities;

          tracer_particle_handler.set_particle_positions(
            relevant_tracer_particle_displacements);
        }

        // Using these new locations, we can then assemble the Stokes system and
        // solve it.
        assemble_stokes_system();
        assemble_nitsche_restriction();
        solve();

        // With the appropriate frequencies, we then write the information of
        // the solid particles, the tracer particles, and the fluid domain into
        // files for visualization, and end the time step by adapting the mesh.
        if (cycle % par.output_frequency == 0)
          {
            output_results(output_cycle, time);
            {
              TimerOutput::Scope t(computing_timer, "Output tracer particles");
              output_particles(tracer_particle_handler,
                               "tracer",
                               output_cycle,
                               time);
            }
            {
              TimerOutput::Scope t(computing_timer, "Output solid particles");
              output_particles(solid_particle_handler,
                               "solid",
                               output_cycle,
                               time);
            }
            ++output_cycle;
          }
        if (cycle % par.refinement_frequency == 0 &&
            cycle != par.number_of_time_steps - 1)
          refine_and_transfer();
      }
  }

} // namespace Step70


// @sect3{The main() function}

// The remainder of the code, the `main()` function, is standard, with the
// exception of the handling of input parameter files. We allow the user to
// specify an optional parameter file as an argument to the program. If
// nothing is specified, we use the default file "parameters.prm", which is
// created if non existent. The file name is scanned for the the string "23"
// first, and "3" afterwards. If the filename contains the string "23", the
// problem classes are instantiated with template arguments 2 and 3
// respectively. If only the string "3" is found, then both template arguments
// are set to 3, otherwise both are set to 2.
//
// If the program is called without any command line arguments (i.e.,
// `argc==1`), then we just use "parameters.prm" by default.
int main(int argc, char *argv[])
{
  using namespace Step70;
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

      if (prm_file.find("23") != std::string::npos)
        {
          StokesImmersedProblemParameters<2, 3> par;
          ParameterAcceptor::initialize(prm_file);

          StokesImmersedProblem<2, 3> problem(par);
          problem.run();
        }
      else if (prm_file.find("3") != std::string::npos)
        {
          StokesImmersedProblemParameters<3> par;
          ParameterAcceptor::initialize(prm_file);

          StokesImmersedProblem<3> problem(par);
          problem.run();
        }
      else
        {
          StokesImmersedProblemParameters<2> par;
          ParameterAcceptor::initialize(prm_file);

          StokesImmersedProblem<2> problem(par);
          problem.run();
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
