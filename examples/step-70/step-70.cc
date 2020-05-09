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
 * Authors: Luca Heltai, Bruno Blais, 2020
 */

// @sect3{Include files}
// Most of these have been introduced elsewhere, we'll comment only on the new
// ones.

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

#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_renumbering.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/fe_nothing.h>
#include <deal.II/fe/mapping_fe_field.h>
#include <deal.II/fe/mapping_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

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

// These are the only new include files w.r.t. step-60. In this tutorial,
// the non-matching coupling between the solid and the fluid is computed using
// an intermediate data structure that keeps track of how the quadrature points
// of the solid evolve w.r.t. the fluid mesh. This data structure needs to keep
// track of the position of the quadrature points on each cell describing the
// solid domain, of the quadrature weights, and possibly of the normal vector
// to each point, if the solid domain is of co-dimension one.
//
// Deal.II offers these facilities in the Particles namespace, through the
// ParticleHandler class. ParticleHandler is a class that allows you to manage
// a collection of particles (objects of type Particles::Particle), representing
// a collection of points with some attached properties (e.g. an id) floating on
// a parallel::distributed::Triangulation. The methods and classes in the
// namespace Particles allows one to easily implement Particle In Cell methods
// and particle tracing on distributed triangulations.
//
// We "abuse" this data structure to store information about the location of
// solid quadrature points w.r.t. to the surrounding fluid grid, including
// integration weights, and possibly surface normals. The reason why we use this
// additional data structure is related to the fact that the solid and the fluid
// grids might be non-overlapping, and are distributed independently among
// processes.
//
// In order to couple the two problems, we rely on the ParticleHandler class,
// storing in each particle the position of a solid quadrature point (which is
// in general not aligned to any of the fluid quadrature points), its weight,
// and any other information that may be required to couple the two problems.
//
// Ownership of the solid quadrature points is inherited by the MPI partitioning
// on the solid mesh itself. The Particles so generated are later distributed to
// the fluid mesh using the methods of the ParticleHandler class. This allows
// transparent exchange of information between MPI processes about the
// overlapping pattern between fluid cells and solid quadrature points.
#include <deal.II/particles/data_out.h>
#include <deal.II/particles/generators.h>
#include <deal.II/particles/particle_handler.h>

// When generating the grids, we allow reading it from a file, and if deal.II
// has been built with OpenCASCADE support, we also allow reading cad files and
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

  // REMOVE THIS FUNCTION ONCE #9891 is merged.
  template <int dim,
            int spacedim,
            typename InputVectorType,
            typename OutputVectorType>
  void interpolate_field_on_particles(
    const DoFHandler<dim, spacedim> &                field_dh,
    const Particles::ParticleHandler<dim, spacedim> &particle_handler,
    const InputVectorType &                          field_vector,
    OutputVectorType &                               interpolated_field,
    const ComponentMask &                            field_comps)
  {
    if (particle_handler.n_locally_owned_particles() == 0)
      {
        interpolated_field.compress(VectorOperation::add);
        return; // nothing else to do here
      }

    const auto &tria     = field_dh.get_triangulation();
    const auto &fe       = field_dh.get_fe();
    auto        particle = particle_handler.begin();

    // Take care of components
    const ComponentMask comps =
      (field_comps.size() == 0 ? ComponentMask(fe.n_components(), true) :
                                 field_comps);
    AssertDimension(comps.size(), fe.n_components());
    const auto n_comps = comps.n_selected_components();

    AssertDimension(field_vector.size(), field_dh.n_dofs());
    AssertDimension(interpolated_field.size(),
                    particle_handler.get_next_free_particle_index() * n_comps);
    // Add check on locally owned indices

    // Global to local indices
    std::vector<unsigned int> space_gtl(fe.n_components(),
                                        numbers::invalid_unsigned_int);
    for (unsigned int i = 0, j = 0; i < space_gtl.size(); ++i)
      if (comps[i])
        space_gtl[i] = j++;

    std::vector<types::global_dof_index> dof_indices(fe.dofs_per_cell);

    while (particle != particle_handler.end())
      {
        const auto &cell = particle->get_surrounding_cell(tria);
        const auto &dh_cell =
          typename DoFHandler<dim, spacedim>::cell_iterator(*cell, &field_dh);
        dh_cell->get_dof_indices(dof_indices);
        const auto pic = particle_handler.particles_in_cell(cell);
        Assert(pic.begin() == particle, ExcInternalError());
        for (unsigned int i = 0; particle != pic.end(); ++particle, ++i)
          {
            const auto &reference_location = particle->get_reference_location();

            const auto id = particle->get_id();

            for (unsigned int j = 0; j < fe.dofs_per_cell; ++j)
              {
                const auto comp_j =
                  space_gtl[fe.system_to_component_index(j).first];
                if (comp_j != numbers::invalid_unsigned_int)
                  interpolated_field[id * n_comps + comp_j] +=
                    fe.shape_value(j, reference_location) *
                    field_vector(dof_indices[j]);
              }
          }
      }
    interpolated_field.compress(VectorOperation::add);
  }

  // Similiarly to what we have done in step-60, we set up a class that holds
  // all the parameters of our problem and derive it from the ParameterAcceptor
  // class to simplify the management and creation of parameter files.
  //
  // The ParameterAcceptor paradigm requires all parameters to be writeable by
  // the ParameterAcceptor methods. In order to avoid bugs that would be very
  // difficult to trace down (such as writing things like `time = 0` instead of
  // `time == 0`), we declare all the parameters in an external class, which is
  // initialized before the actual StokesImmersedProblem class, and pass it to
  // the main class as a const reference.
  template <int dim, int spacedim = dim>
  class StokesImmersedProblemParameters : public ParameterAcceptor
  {
  public:
    // The constructor is responsible for the connection between the members of
    // this class and the corresponding entries in the ParameterHandler. Thanks
    // to the use of the ParameterHandler::add_parameter() method, this
    // connection is trivial, but requires all members of this class to be
    // writeable
    StokesImmersedProblemParameters();

    // however, since this class will be passed as a const reference to the
    // StokesImmersedProblem class, we have to make sure we can still set the
    // time correctly in the objects derived by the Function class defined
    // here. In order to do so, we declare both the
    // StokesImmersedProblemParameters::rhs and
    // StokesImmersedProblemParameters::angular_velocity members to be mutable,
    // and define this little helper method that sets their time to the correct
    // value.
    void set_time(const double &time) const
    {
      rhs.set_time(time);
      angular_velocity.set_time(time);
    }

    // We will use a Taylor-Hood function space of arbitrary order. This
    // parameter is used to initialize the FiniteElement space with the corret
    // FESystem object
    unsigned int velocity_degree = 2;

    // Instead of defining a time step increment, in this tutorial we prefer to
    // let the user choose a final simulation time, and the number of steps in
    // which we want to reach the final time
    unsigned int number_of_time_steps = 1;
    double       final_time           = 1.0;

    // Instead of producing an output at every time step, we allow the user to
    // select the frequency at which output is produced:
    unsigned int output_frequency = 1;

    // We allow every grid to be refined independently. In this tutorial, no
    // physics is resolved on the solid grid, and its velocity is given as a
    // datum. However it is relatively straightforward to incorporate some
    // elasticity model in this tutorial, and transform it into a fully fledged
    // FSI solver.
    unsigned int initial_fluid_refinement      = 3;
    unsigned int initial_solid_refinement      = 3;
    unsigned int particle_insertion_refinement = 1;

    // The only two parameters used in the equations are the viscosity of the
    // fluid, and the penalty term used in the Nitsche formulation:
    double viscosity    = 1.0;
    double penalty_term = 10;

    // By default, we create a hyper_cube without colorisation, and we use
    // homogenous Dirichlet boundary conditions. In this set we store the
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
    // maps, tuples) and basic dealii types (Point, Tensor, BoundingBox, etc.).
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
    // geometry of the domain. Every CAD file will be analysed and a Manifold of
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
    // compatible with the solid grid. Additionnaly, we perform local refinement
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
    int          max_level_refinement = 5;
    std::string  refinement_strategy  = "fixed_fraction";
    double       coarsening_fraction  = 0.3;
    double       refinement_fraction  = 0.3;
    unsigned int max_cells            = 1000;
    int          refinement_frequency = 5;

    // These two functions are used to control the source term of Stokes flow
    // and the angular velocity at which we move the solid body. In a more
    // realistic simulation, the solid velocity or its deformation would come
    // from the solution of an auxiliary problem on the solid domain. In this
    // example step we leave this part aside, and simply impose a fixed
    // rotational velocity field on the immersed solid, governed by function
    // that can be specified in the parameter file:
    mutable ParameterAcceptorProxy<Functions::ParsedFunction<spacedim>> rhs;
    mutable ParameterAcceptorProxy<Functions::ParsedFunction<spacedim>>
      angular_velocity;
  };


  // Once the angular velocity is provided as a Function object, we reconstruct
  // the pointwise solid velocity through the following class which derives
  // from the Function class. It derives the value of the velocity of
  // the solid body using the position of the points.
  template <int spacedim>
  class SolidVelocity : public Function<spacedim>
  {
  public:
    SolidVelocity(const Functions::ParsedFunction<spacedim> &angular_velocity)
      : angular_velocity(angular_velocity)
    {
      static_assert(spacedim > 1,
                    "Cannot instantiate SolidVelocity for spacedim == 1");
    }

    virtual double value(const Point<spacedim> &p,
                         unsigned int           component = 0) const override
    {
      Tensor<1, spacedim> velocity;
      if (spacedim == 3)
        {
          Tensor<1, spacedim> omega;
          for (unsigned int i = 0; i < spacedim; ++i)
            omega[i] = angular_velocity.value(p, i);

          velocity = cross_product_3d(p, omega);
        }
      else if (spacedim == 2)
        {
          const double omega = angular_velocity.value(p, 0);

          velocity[0] = -omega * p[1];
          velocity[1] = omega * p[0];
        }

      return velocity[component];
    }

  private:
    const Functions::ParsedFunction<spacedim> &angular_velocity;
  };

  // Similarly, we assume that the incremental solid displacement can be
  // computed simply by a one step time integration process (here using a
  // trivial forward Euler method), so that at each time step, the solid simply
  // displaces by `v*dt`.
  template <int spacedim>
  class SolidDisplacement : public Function<spacedim>
  {
  public:
    SolidDisplacement(
      const Functions::ParsedFunction<spacedim> &angular_velocity,
      const double                               time_step)
      : Function<spacedim>(spacedim)
      , angular_velocity(angular_velocity)
      , time_step(time_step)
    {
      static_assert(spacedim > 1,
                    "Cannot instantiate SolidDisplacement for spacedim == 1");
    }

    virtual double value(const Point<spacedim> &p,
                         unsigned int           component = 0) const override
    {
      Tensor<1, spacedim> displacement;

      double dtheta = angular_velocity.value(p, 0) * time_step;

      displacement[0] = std::cos(dtheta) * p[0] - std::sin(dtheta) * p[1];
      displacement[1] = std::sin(dtheta) * p[0] + std::cos(dtheta) * p[1];

      return displacement[component];
    }

    void set_time_step(const double new_time_step)
    {
      time_step = new_time_step;
    }

  private:
    const Functions::ParsedFunction<spacedim> &angular_velocity;
    double                                     time_step;
  };

  // We are now ready to introduce the main class of our tutorial program.
  template <int dim, int spacedim = dim>
  class StokesImmersedProblem
  {
  public:
    StokesImmersedProblem(
      const StokesImmersedProblemParameters<dim, spacedim> &par);

    // As usual, we leave a single public entry point to the user: the run
    // method. Everything else is left private, and accessed through the run
    // method itself.
    void run();

  private:
    // This method is similar to what is present in previous example. However
    // it not only takes care of generating the grid for the fluid, but also
    // the grid for the solid.
    void make_grid();

    // These two methods are new w.r.t. previous examples, and initialize the
    // Particles::ParticleHandler objects used in this class. We have two such
    // objects: one represents passive tracers, used to plot the trajectories
    // of fluid particles, while the the other represents material particles
    // of the solid, which are placed at quadrature points of the solid grid.
    void setup_tracer_particles();
    void setup_solid_particles();

    // The setup is split in two parts: create all objects that are needed once
    // per simulation,
    void initial_setup();
    // followed by all objects that need to be reinitialized at every refinement
    // step.
    void setup_dofs();

    // The assembly routine is very similar to other Stokes assembly routines,
    void assemble_stokes_system();
    // with the exception of the Nistche restriction part, which exploits one of
    // the particle handlers to integrate on a non-matching part of the fluid
    // domain, corresponding to the position of the solid.
    void assemble_nitsche_restriction();

    // Nothing new in the solve routine, which is almost identical to step-60
    void solve();

    // The refine_and_transfer() method is called only every
    // `refinement_frequency` steps, and makes sure that all the fields
    // that were computed on the time step before refinement are transfered
    // correctly to the new grid. This includes vector fields, as well as
    // particle information.
    void refine_and_transfer();

    // Similarly, we call the output_results() method only every
    // `output_frequency` steps. This method takes care of outputting both the
    // fields variables,
    void output_results(const unsigned int cycle, const double time) const;

    // and the tracers:
    void output_particles(const Particles::ParticleHandler<spacedim> &particles,
                          std::string                                 fprefix,
                          const unsigned int                          iter,
                          const double                                time) const;

    // As noted before, we make sure we cannot modify this object from within
    // this class, by making it a const reference.
    const StokesImmersedProblemParameters<dim, spacedim> &par;

    MPI_Comm mpi_communicator;

    // For the current implementation, only `fluid_fe` is really necessary.
    // For completeness, and to allow easy extensions, we keep also the
    // `solid_fe` around, which is however initialized to a FE_Nothing finite
    // element space, i.e., one that has no degrees of freedom.
    //
    // We declare both finite element spaces as unique pointers, to allow their
    // generation after StokesImmersedProblemParameters has been initialized. In
    // particular, they will be initialized in the initial_setup() method
    std::unique_ptr<FiniteElement<spacedim>>      fluid_fe;
    std::unique_ptr<FiniteElement<dim, spacedim>> solid_fe;

    // This is one of the main novelties w.r.t. the tutorial step-60. Here we
    // assume that both the solid and the fluid are fully distributed
    // triangulations. This allows the problem to scale to a very large number
    // of degrees of freedom, at the cost of communicating all the overlapping
    // regions between non matching triangulations. This is especially tricky,
    // since we make no assumptions on the relative position or distribution of
    // the various subdomains. In particular, we assume that every process owns
    // only a part of the solid_tria, and only a part of the fluid_tria, not
    // necessarily in the same physical region, and not necessarily overlapping.
    //
    // We could in principle try to create the initial subdivisions in such a
    // way that they overlap between the solid and the fluid regions. However,
    // this overlap would be destroyed during the simulation, and we would have
    // to redistribute the dofs again and again. The approach we follow in this
    // tutorial is more flexible, and not much more expensive. We make two
    // all-to-all communications at the beginning of the simulation to exchange
    // information about an (approximate) information of the geometrical
    // occupancy of each processor (done through a collection of bounding
    // boxes).
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

    AffineConstraints<double> constraints;

    LA::MPI::BlockSparseMatrix system_matrix;
    LA::MPI::BlockSparseMatrix coupling_matrix;

    LA::MPI::BlockSparseMatrix preconditioner_matrix;
    LA::MPI::BlockVector       solution;
    LA::MPI::BlockVector       locally_relevant_solution;
    LA::MPI::BlockVector       system_rhs;

    // For every tracer particle, we need to compute the velocity field in its
    // current position, and update its position using a discrete time stepping
    // scheme. We do this using distributed linear algebra objects, where the
    // owner of a particle is set to be equal to the process that generated that
    // particle at time t=0. This information is stored for every process in the
    // `owned_tracer_particles` IndexSet, that indicates which particles the
    // current process owns.
    //
    // Once the particles have been distributed around to match the process that
    // owns the region where the particle lives, we will need read access from
    // that process on the corresponding velocity field. We achieve this by
    // filling a read only velocity vector field, that contains the relevant
    // information in ghost entries. This is achieved using the
    // `relevant_tracer_particles` IndexSet, that keeps track of how things
    // change during the simulation, i.e., it keeps track of where particles
    // that I own have ended up being, and who owns the particles that ended up
    // in my subdomain.
    //
    // While this is not the most efficient strategy, we keep it this way to
    // illustrate how things would work in a real FSI problem. If a particle
    // is linked to a specific solid degree of freedom, we are not free to
    // choose who owns it, and we have to communicate this information around.
    // We illustrate this here, and show that the communication pattern is
    // point-to-point, and negligible in terms of total cost of the algorithm.
    IndexSet owned_tracer_particles;
    IndexSet relevant_tracer_particles;

    // These vectors are used to store the particles velocities (read-only, with
    // ghost entries) and their displacement (read/write, no ghost entries).
    LA::MPI::Vector tracer_particle_velocities;
    LA::MPI::Vector relevant_tracer_particle_displacements;

    // We fix once the quadrature formula that is used to integrate on the fluid
    // and on the solid domains.
    std::unique_ptr<Quadrature<spacedim>> fluid_quadrature_formula;
    std::unique_ptr<Quadrature<dim>>      solid_quadrature_formula;

    // Finally, these are the two Particles::ParticleHandler classes used to
    // couple the solid with the fluid, and to describe the passive tracers.
    Particles::ParticleHandler<spacedim> tracer_particle_handler;
    Particles::ParticleHandler<spacedim> solid_particle_handler;

    ConditionalOStream  pcout;
    mutable TimerOutput computing_timer;
  };


  // In the constructor, we create the mpi_communicator as well as
  // the triangulations and dof_handler for both the fluid and the solid.
  // Using the mpi_communicator, both the ConditionalOSStream and TimerOutput
  // are constructed.
  template <int dim, int spacedim>
  StokesImmersedProblem<dim, spacedim>::StokesImmersedProblem(
    const StokesImmersedProblemParameters<dim, spacedim> &par)
    : par(par)
    , mpi_communicator(MPI_COMM_WORLD)
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
    , pcout(std::cout,
            (Utilities::MPI::this_mpi_process(mpi_communicator) == 0))
    , computing_timer(mpi_communicator,
                      pcout,
                      TimerOutput::summary,
                      TimerOutput::wall_times)
  {}


  // In order to generate the grid, we first try to use the functions in the
  // deal.II GridGenerator namespace, by leveraging the
  // GridGenerator::generate_from_name_and_argument(), if this function fails,
  // then we use the following method, where the name is interpreted as a
  // filename, and the arguments are interpreted as a map from manifold ids to
  // CAD files, and are converted to Manifold descriptors using the OpenCASCADE
  // namespace facilities:
  template <int dim, int spacedim>
  void read_grid_and_cad_files(const std::string &grid_file_name,
                               const std::string &ids_and_cad_file_names,
                               Triangulation<dim, spacedim> &tria)
  {
    // Try to read the grid using GridIn facilities. Let the GridIn class
    // automatically detect the file format:
    GridIn<dim, spacedim> grid_in;
    grid_in.attach_triangulation(tria);
    grid_in.read(grid_file_name);

    // If we got to this point, then the Triangulation has been read, and we are
    // ready to attach to it the correct manifold descriptions. We perform the
    // next lines of codes only if deal.II has been built with OpenCASCADE
    // support. For each entry in the map, we try to open the corresponding CAD
    // file, we analyse it, and according to its content, opt for either a
    // ArchLengthProjectionLineManifold (if the CAD file contains a single
    // TopoDS_Edge or a single TopoDS_Wire) or a NURBSPatchManifold, if the file
    // contains a single face. Notice that if the CAD files do not contain
    // single wires, edges, or faces, an assertion will be throw in the
    // generation of the Manifold.
    //
    // We use the Patterns::Tools::Convert class to do the convertion from the
    // string to a map between manifold ids and file names for us:

#ifdef DEAL_II_WITH_OPENCASCADE
    using map_type  = std::map<types::manifold_id, std::string>;
    using Converter = Patterns::Tools::Convert<map_type>;

    for (const auto pair : Converter::to_value(ids_and_cad_file_names))
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

        // Now we check how many faces are contained in the Shape. OpenCASCADE
        // is intrinsically 3D, so if this number is zero, we interpret this as
        // a line manifold, otherwise as a NormalToMeshProjectionManifold in
        // spacedim = 3, or NURBSPatchManifold in spacedim = 2.
        const auto n_elements = OpenCASCADE::count_elements(shape);
        if ((std::get<0>(n_elements) == 0))
          tria.set_manifold(
            manifold_id,
            OpenCASCADE::ArclengthProjectionLineManifold<dim, spacedim>(shape));
        else if (spacedim == 3)
          {
            // We use this trick, because NormalToMeshProjectionManifold
            // is only implemented for spacedim = 3. The check above makes
            // sure that things actually work correctly.
            const auto t = reinterpret_cast<Triangulation<dim, 3> *>(&tria);
            t->set_manifold(manifold_id,
                            OpenCASCADE::NormalToMeshProjectionManifold<dim, 3>(
                              shape));
          }
        else
          // We also allow surface descriptions in two dimensional spaces based
          // on single NURBS patches. For this to work, the CAD file must
          // contain a single TopoDS_Face.
          tria.set_manifold(manifold_id,
                            OpenCASCADE::NURBSPatchManifold<dim, spacedim>(
                              TopoDS::Face(shape)));
      }
#else
    (void)ids_and_cad_file_names;
#endif
  }

  // Now let's put things together, and make all the necessary grids
  template <int dim, int spacedim>
  void StokesImmersedProblem<dim, spacedim>::make_grid()
  {
    try
      {
        // we first try to generate the grid internally
        GridGenerator::generate_from_name_and_arguments(
          fluid_tria, par.name_of_fluid_grid, par.arguments_for_fluid_grid);
      }
    catch (...)
      {
        // and if we fail, we proceed with the above function call
        read_grid_and_cad_files(par.name_of_fluid_grid,
                                par.arguments_for_fluid_grid,
                                fluid_tria);
      }
    fluid_tria.refine_global(par.initial_fluid_refinement);

    // The same is done for the solid grid:
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

  // Once the solid and fluid grids have been created, we start filling the
  // Particles::ParticleHandler objects. The first one we take care of is the
  // one we use to keep track of passive tracers in the fluid. No other role is
  // given to these particles, so their initialization is standard.
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
  // on an arbitrary grid (non-matching w.r.t. to the fluid grid) have no
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
  //
  // In the approach used here, ownership of the particles is assigned once at
  // the beginning, and one-to-one communication happen whenever the original
  // owner needs information from the process that owns the cell where the
  // particle lives. We make sure that we set ownership of the particles using
  // the initial particle distribution, and keep the same ownership throughout
  // the execution of the program.
  template <int dim, int spacedim>
  void StokesImmersedProblem<dim, spacedim>::setup_tracer_particles()
  {
    // Generate a triangulation that will be used to decide the position
    // of the particles to insert.
    parallel::distributed::Triangulation<spacedim> particle_insert_tria(
      mpi_communicator);
    GridGenerator::generate_from_name_and_arguments(
      particle_insert_tria,
      par.name_of_particle_grid,
      par.arguments_for_particle_grid);
    particle_insert_tria.refine_global(par.particle_insertion_refinement);

    // Generate the support point on the triangulation that will be used as
    // particle insertion point
    DoFHandler<spacedim> particles_dof_handler(particle_insert_tria);
    FE_Q<spacedim>       particles_fe(1);
    particles_dof_handler.distribute_dofs(particles_fe);

    // Create the particle handler associated with the fluid triangulation
    tracer_particle_handler.initialize(fluid_tria,
                                       StaticMappingQ1<spacedim>::mapping);

    // This is where things start to get complicated. In fully distributed
    // triangulations, the active process only knows about the locally owned
    // cells, and has no idea of how other processes have distributed their own
    // cells. On the other hand, by design we assign to the active process also
    // the particles that were generated as support points of the locally owned
    // parts of a non matching and arbitrary grid.
    //
    // The location of these particles is arbitrary, and may fall within a
    // region that we don't have access to (i.e., a region of the fluid domain
    // where cells are artificial). In order to understand who to send those
    // particles to, we need to have a (rough) idea of how the fluid grid is
    // distributed among processors.
    //
    // The following method computes a small collection of axis aligned
    // bounding boxes that represent well the locally owned part of the fluid
    // domain.
    auto bounding_boxes_of_locally_owned_cells =
      GridTools::compute_mesh_predicate_bounding_box(
        fluid_tria, IteratorFilters::LocallyOwnedCell());

    // These bounding boxes are then exchanged with all other processes, so
    // that now every active process has a rough idea (based on bounding boxes
    // that overlap with the locally owned part of the domain) of who they
    // should send particles to.
    auto global_bounding_boxes =
      Utilities::MPI::all_gather(MPI_COMM_WORLD,
                                 bounding_boxes_of_locally_owned_cells);


    // Finally generate the particles from the support points of the
    // tracer particles triangulation. This function call uses the
    // global_bounding_boxes object we just constructed. At the end of this
    // call, every particle will have been distributed to the correct process
    // (i.e., the process that owns the cell where the particle lives).
    Particles::Generators::dof_support_points(particles_dof_handler,
                                              global_bounding_boxes,
                                              tracer_particle_handler);

    // As soon as we have initialized the particles in each process, we set
    // their ownership to the current distribution. Since we will need this
    // information to initialize a vector containing velocity information for
    // each particle, we query the tracer_particle_handler for its locally
    // relevant ids, and construct the indices that would be needed to store in
    // a (parallel distributed) vector the position and velocity of all
    // particles.
    owned_tracer_particles =
      tracer_particle_handler.locally_relevant_ids().tensor_product(
        complete_index_set(spacedim));

    // At the beginning of the simulation, all particles are in their original
    // position. When particles move, they may jump to a part of the domain
    // which is owned by another process. If this happens, the current process
    // keeps formally "ownership" of the particles, but may need read access
    // from the process where the particle has landed. We keep this information
    // in another index set, which stores the indices of all particles that are
    // currently on my subdomain, independently if they have always been here or
    // not.
    //
    // Keeping this index set around allows us to leverage linear algebra
    // classes for all communications regarding positions and velocities of the
    // particles. This mimics what would happen in the case where another
    // problem was solved in the solid domain (as in fluid-structure
    // interaction. In this latter case, additional DOFs on the solid domain
    // would be coupled to what is occuring in the fluid domain.
    relevant_tracer_particles = owned_tracer_particles;

    // Now make sure that upon refinement, particles are correctly transferred.
    // When performing local refinement or coarsening, particles will land in
    // another cell. We could in principle redistribute all particles after
    // refining, however this would be overly expensive.
    //
    // The Particles::ParticleHandler class has a way to transfer information
    // from a cell to its children or to its parent upon refinement, without the
    // need to reconstruct the entire data structure. This is done by
    // "registering" two callback functions to the triangulation. These
    // functions will receive a signal when refinement is about to happen, and
    // when it has just happened, and will take care of transferring all
    // information to the newly refined grid with minimal computational cost.
    fluid_tria.signals.pre_distributed_refinement.connect(std::bind(
      &Particles::ParticleHandler<spacedim>::register_store_callback_function,
      &tracer_particle_handler));

    fluid_tria.signals.post_distributed_refinement.connect(std::bind(
      &Particles::ParticleHandler<spacedim>::register_load_callback_function,
      &tracer_particle_handler,
      false));

    // Finally, we display to the terminal the number of total tracer particles
    // that were generated
    pcout << "Tracer particles: "
          << tracer_particle_handler.n_global_particles() << std::endl;
  }


  // Similarly to what we have done for passive tracers, we now setup the solid
  // particles. The main difference here is that we also want to attach a weight
  // value to each of the quadrature points, so that we can compute integrals
  // even without direct access to the original solid grid.
  //
  // This is achieved by leveraging the "properties" concept of the
  // Particles::ParticleHandler class. It is possible to store (in a memory
  // efficient way) an arbitrary number of doubles alongside with a
  // Particles::Particle, inside a Particles::ParticleHandler object. We use
  // this possibility to store the JxW values of the quadrature points of the
  // solid grid.
  template <int dim, int spacedim>
  void StokesImmersedProblem<dim, spacedim>::setup_solid_particles()
  {
    QGauss<dim> quadrature(fluid_fe->degree + 1);

    // We only need to store one property per particle: the JxW value of the
    // integration on the solid grid. This is passed at construction time to the
    // solid_particle_handler object as the last argument
    const unsigned int n_properties = 1;
    solid_particle_handler.initialize(fluid_tria,
                                      StaticMappingQ1<spacedim>::mapping,
                                      n_properties);

    // The number of particles that we generate locally is equal to the total
    // number of cells times the number of quadrature points used in each cell.
    // We store all these points in a vector, and their corresponding properties
    // in a vector of vectors, which we fill with the actual position of the
    // quadrature points, and their JxW values (the only information that is
    // needed to perform integration on the solid cells, even with a
    // non-matching grid).
    std::vector<Point<spacedim>> quadrature_points_vec(
      quadrature.size() * solid_tria.n_locally_owned_active_cells());

    std::vector<std::vector<double>> properties(
      quadrature.size() * solid_tria.n_locally_owned_active_cells(),
      std::vector<double>(n_properties));

    FEValues<dim, spacedim> fe_v(*solid_fe,
                                 quadrature,
                                 update_JxW_values | update_quadrature_points);

    unsigned int point_index = 0;
    for (const auto &cell : solid_dh.active_cell_iterators())
      if (cell->is_locally_owned())
        {
          fe_v.reinit(cell);
          const auto &points = fe_v.get_quadrature_points();
          const auto &JxW    = fe_v.get_JxW_values();

          for (unsigned int q = 0; q < points.size(); ++q)
            {
              quadrature_points_vec[point_index] = points[q];
              properties[point_index][0]         = JxW[q];
              ++point_index;
            }
        }

    // We proceed in the same way we did with the tracer particles
    auto my_bounding_box = GridTools::compute_mesh_predicate_bounding_box(
      fluid_tria, IteratorFilters::LocallyOwnedCell());

    auto global_bounding_boxes =
      Utilities::MPI::all_gather(mpi_communicator, my_bounding_box);

    // Since we have already stored the position of the quadrature point,
    // we can use these positions to insert the particles directly using
    // the solid_particle_handler instead of having to go through a
    // Particles::Generators
    auto cpu_to_index =
      solid_particle_handler.insert_global_particles(quadrature_points_vec,
                                                     global_bounding_boxes,
                                                     properties);


    // Now make sure that upon refinement, particles are correctly transferred
    fluid_tria.signals.pre_distributed_refinement.connect(std::bind(
      &Particles::ParticleHandler<spacedim>::register_store_callback_function,
      &solid_particle_handler));

    fluid_tria.signals.post_distributed_refinement.connect(std::bind(
      &Particles::ParticleHandler<spacedim>::register_load_callback_function,
      &solid_particle_handler,
      false));

    pcout << "Solid particles: " << solid_particle_handler.n_global_particles()
          << std::endl;
  }

  // We set up the finite element space and the quadrature formula to be
  // used throughout the step. For the fluid, we use Taylor-Hood elements (e.g.
  // Q2-Q1). Since we do not solve any equation on the solid domain, an empty
  // finite element space is generated. A natural extension of this program
  // would be to solve a fluid structure interaction problem, which would
  // require that the solid_fe use a non-empty FiniteElement.
  template <int dim, int spacedim>
  void StokesImmersedProblem<dim, spacedim>::initial_setup()
  {
    // We store the time necessary to carry-out the initial_setup under the
    // label "initial setup" Numerous other calls to this timer are made in
    // various functions. They allow to monitor the absolute and relative load
    // of each individual function to identify the bottlenecks.
    TimerOutput::Scope t(computing_timer, "initial setup");

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


  // We construct the distributed block matrices and vectors which are used to
  // solve the linear equations that arise from the problem. This function is
  // adapted from step-55 and we refer to this step for a thorough explanation.
  template <int dim, int spacedim>
  void StokesImmersedProblem<dim, spacedim>::setup_dofs()
  {
    TimerOutput::Scope t(computing_timer, "setup dofs");

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
        ZeroFunction<spacedim>(spacedim + 1),
        constraints,
        fluid_fe->component_mask(velocities));
      constraints.close();
    }

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
        fluid_dh.compute_locally_owned_dofs_per_processor(),
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
        fluid_dh.compute_locally_owned_dofs_per_processor(),
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


  // We assemble the system matrix, the preconditioner matrix, and the right
  // hand side. The code is adapted from step-55 and is pretty standard.
  template <int dim, int spacedim>
  void StokesImmersedProblem<dim, spacedim>::assemble_stokes_system()
  {
    system_matrix         = 0;
    preconditioner_matrix = 0;
    system_rhs            = 0;

    TimerOutput::Scope t(computing_timer, "Stokes_assembly");


    FEValues<spacedim> fe_values(*fluid_fe,
                                 *fluid_quadrature_formula,
                                 update_values | update_gradients |
                                   update_quadrature_points |
                                   update_JxW_values);

    const unsigned int dofs_per_cell = fluid_fe->dofs_per_cell;
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


  // This method is the heart of the tutorial, but it is relatively
  // straightforward. Here we exploit the solid_particle_handler to compute the
  // Nitsche restriction or the penalization in the embedded domain.
  template <int dim, int spacedim>
  void StokesImmersedProblem<dim, spacedim>::assemble_nitsche_restriction()
  {
    TimerOutput::Scope t(computing_timer, "Nitsche_assembly");

    const FEValuesExtractors::Vector velocities(0);
    const FEValuesExtractors::Scalar pressure(spacedim);

    SolidVelocity<spacedim> solid_velocity(par.angular_velocity);

    std::vector<types::global_dof_index> fluid_dof_indices(
      fluid_fe->dofs_per_cell);

    FullMatrix<double>     local_matrix(fluid_fe->dofs_per_cell,
                                    fluid_fe->dofs_per_cell);
    dealii::Vector<double> local_rhs(fluid_fe->dofs_per_cell);

    const auto k = 1.0 / GridTools::minimal_cell_diameter(fluid_tria);

    // We loop over all the local particles. Although this could be achieved
    // directly by looping over all the cells, this would force us
    // to loop over numerous cells which do not contain particles.
    // Consequently, we loop over all the particles, but, we get the reference
    // of the cell in which the particle lies and then loop over all particles
    // within that cell. This enables us to skip the cells which do not contain
    // particles, yet to assemble the local matrix and rhs of each cell to apply
    // the Nitsche restriction
    auto particle = solid_particle_handler.begin();
    while (particle != solid_particle_handler.end())
      {
        local_matrix = 0;
        local_rhs    = 0;

        // We get the refence to the cell within which the particle lies from
        // the particle itself. Consequently, we can assemble the additional
        // terms in the system matrix the rhs as we would normally.
        const auto &cell = particle->get_surrounding_cell(fluid_tria);
        const auto &dh_cell =
          typename DoFHandler<spacedim>::cell_iterator(*cell, &fluid_dh);
        dh_cell->get_dof_indices(fluid_dof_indices);

        const auto pic = solid_particle_handler.particles_in_cell(cell);

        Assert(pic.begin() == particle, ExcInternalError());

        for (const auto &p : pic)
          {
            // From the particle we gather the location and the reference
            // location of the particle as well as the additional information
            // that is attached to the particle. In the present case, this
            // information is the JxW of the quadrature points which were used
            // to generate the particles.
            const auto &ref_q  = p.get_reference_location();
            const auto &real_q = p.get_location();
            const auto &JxW    = p.get_properties()[0];

            // We add the contribution of the quadrature point to the
            // local_matrix and local_rhs. We can evaluate the value of the
            // shape function at the position of each particle easily by using
            // its reference location.
            for (unsigned int i = 0; i < fluid_fe->dofs_per_cell; ++i)
              {
                const auto comp_i =
                  fluid_fe->system_to_component_index(i).first;
                if (comp_i < spacedim)
                  {
                    for (unsigned int j = 0; j < fluid_fe->dofs_per_cell; ++j)
                      {
                        const auto comp_j =
                          fluid_fe->system_to_component_index(j).first;
                        if (comp_i == comp_j)
                          local_matrix(i, j) +=
                            k * par.penalty_term *
                            fluid_fe->shape_value(i, ref_q) *
                            fluid_fe->shape_value(j, ref_q) * JxW;
                      }
                    local_rhs(i) += k * par.penalty_term *
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
  }


  // This function solves the linear system with MINRES with a block diagonal
  // preconditioner and AMG for the two diagonal blocks as used in step-55. The
  // preconditioner applies a v cycle to the 0,0 block and a CG with the mass
  // matrix for the 1,1 block (the Schur complement).
  template <int dim, int spacedim>
  void StokesImmersedProblem<dim, spacedim>::solve()
  {
    TimerOutput::Scope t(computing_timer, "solve");

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

    const auto P =
      block_diagonal_operator<2, LA::MPI::BlockVector>({amgA, amgS});

    SolverControl solver_control(system_matrix.m(),
                                 1e-10 * system_rhs.l2_norm());

    SolverMinRes<LA::MPI::BlockVector> solver(solver_control);

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


  // We deal with mesh refinement in a standard way.
  template <int dim, int spacedim>
  void StokesImmersedProblem<dim, spacedim>::refine_and_transfer()
  {
    TimerOutput::Scope               t(computing_timer, "refine");
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
      if (cell->refine_flag_set() && cell->level() == par.max_level_refinement)
        cell->clear_refine_flag();

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



  // We output the results (velocity and pressure) on the fluid domain
  // using the standard parallel capacities of deal.II. A single compressed vtu
  // file is written that agglomerates the information of all processors. An
  // additional .pvd record is written to associate the physical time to the vtu
  // files.
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

    LA::MPI::BlockVector interpolated;
    interpolated.reinit(fluid_owned_dofs, MPI_COMM_WORLD);
    VectorTools::interpolate(fluid_dh,
                             ConstantFunction<spacedim>(1.0, spacedim + 1),
                             interpolated);

    LA::MPI::BlockVector interpolated_relevant(fluid_owned_dofs,
                                               fluid_relevant_dofs,
                                               MPI_COMM_WORLD);
    interpolated_relevant = interpolated;
    {
      std::vector<std::string> solution_names(spacedim, "ref_u");
      solution_names.emplace_back("ref_p");
      data_out.add_data_vector(interpolated_relevant,
                               solution_names,
                               DataOut<spacedim>::type_dof_data,
                               data_component_interpretation);
    }


    Vector<float> subdomain(fluid_tria.n_active_cells());
    for (unsigned int i = 0; i < subdomain.size(); ++i)
      subdomain(i) = fluid_tria.locally_owned_subdomain();
    data_out.add_data_vector(subdomain, "subdomain");

    data_out.build_patches();

    const std::string filename =
      "solution-" + Utilities::int_to_string(cycle) + ".vtu";
    data_out.write_vtu_in_parallel(filename, mpi_communicator);

    static std::vector<std::pair<double, std::string>> times_and_names;
    times_and_names.push_back(std::make_pair(time, filename));
    std::ofstream ofile("solution.pvd");
    DataOutBase::write_pvd_record(ofile, times_and_names);
  }

  // We write the particles (either from the solid or the tracers)
  // as a single compressed vtu file through the Particles::DataOut object.
  // This simple object does not write the additional information
  // attached to the particles, but only writes their id.
  template <int dim, int spacedim>
  void StokesImmersedProblem<dim, spacedim>::output_particles(
    const Particles::ParticleHandler<spacedim> &particles,
    std::string                                 fprefix,
    const unsigned int                          iter,
    const double                                time) const
  {
    Particles::DataOut<spacedim, spacedim> particles_out;
    particles_out.build_patches(particles);
    const std::string filename =
      (fprefix + "-" + Utilities::int_to_string(iter) + ".vtu");
    particles_out.write_vtu_in_parallel(filename, mpi_communicator);


    static std::map<std::string, std::vector<std::pair<double, std::string>>>
      times_and_names;
    if (times_and_names.find(fprefix) != times_and_names.end())
      times_and_names[fprefix].push_back(std::make_pair(time, filename));
    else
      times_and_names[fprefix] = {std::make_pair(time, filename)};
    std::ofstream ofile(fprefix + ".pvd");
    DataOutBase::write_pvd_record(ofile, times_and_names[fprefix]);
  }


  // This function orchestrates the entire simulation. It is very similar
  // to the other transient steps.
  template <int dim, int spacedim>
  void StokesImmersedProblem<dim, spacedim>::run()
  {
#ifdef USE_PETSC_LA
    pcout << "Running using PETSc." << std::endl;
#else
    pcout << "Running using Trilinos." << std::endl;
#endif

    ComponentMask velocity_mask(spacedim + 1, true);
    velocity_mask.set(spacedim, false);

    const double time_step = par.final_time / (par.number_of_time_steps - 1);
    double       time      = 0;
    for (unsigned int cycle = 0; cycle < par.number_of_time_steps;
         ++cycle, time += time_step)
      {
        par.set_time(time);
        pcout << "Cycle " << cycle << ':' << std::endl
              << "Time : " << time << ", time step: " << time_step << std::endl;

        // We initialize all the elements of the simulation in the first cycle
        if (cycle == 0)
          {
            make_grid();
            initial_setup();
            setup_dofs();
            setup_tracer_particles();
            setup_solid_particles();
            tracer_particle_velocities.reinit(owned_tracer_particles,
                                              mpi_communicator);
            output_results(0, time);
            {
              TimerOutput::Scope t(computing_timer, "Output tracer particles");
              output_particles(tracer_particle_handler, "tracer", 0, time);
            }
            {
              TimerOutput::Scope t(computing_timer, "Output solid particles");
              output_particles(solid_particle_handler, "solid", 0, time);
            }
          }
        // On the other cycle, we displace the solid body to take into account
        // the fact that is has moved.
        else
          {
            TimerOutput::Scope t(computing_timer,
                                 "Set solid particle position");

            SolidDisplacement<spacedim> solid_displacement(par.angular_velocity,
                                                           time_step);
            solid_particle_handler.set_particle_positions(solid_displacement,
                                                          false);
          }
        {
          // We interpolate the fluid velocity at the position of the tracer
          // particles and, with a naive explicit Euler scheme, we advect the
          // massless tracer particles.
          TimerOutput::Scope t(computing_timer, "Set tracer particle motion");
          interpolate_field_on_particles(fluid_dh,
                                         tracer_particle_handler,
                                         locally_relevant_solution,
                                         tracer_particle_velocities,
                                         velocity_mask);

          tracer_particle_velocities *= time_step;

          relevant_tracer_particles =
            tracer_particle_handler.locally_relevant_ids().tensor_product(
              complete_index_set(spacedim));

          relevant_tracer_particle_displacements.reinit(
            owned_tracer_particles,
            relevant_tracer_particles,
            mpi_communicator);

          relevant_tracer_particle_displacements = tracer_particle_velocities;

          tracer_particle_handler.set_particle_positions(
            relevant_tracer_particle_displacements);
        }
        assemble_stokes_system();
        assemble_nitsche_restriction();
        solve();

        // At every output frequency, we write the information of the solid
        // particles, the tracer particles and the fluid domain.
        if (cycle % par.output_frequency == 0)
          {
            static unsigned int output_cycle = 0;
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


  template <int dim, int spacedim>
  StokesImmersedProblemParameters<dim,
                                  spacedim>::StokesImmersedProblemParameters()
    : ParameterAcceptor("Stokes Immersed Problem/")
    , rhs("Right hand side", spacedim + 1)
    , angular_velocity("Angular velocity", spacedim == 3 ? spacedim : 1)
  {
    // We split the parameters in various categories, by putting them in
    // different sections of the ParameterHandler class. We begin by
    // declaring all the global parameters used by StokesImmersedProblem
    // in the global scope:
    add_parameter(
      "Velocity degree", velocity_degree, "", this->prm, Patterns::Integer(1));

    add_parameter("Number of time steps", number_of_time_steps);
    add_parameter("Output frequency", output_frequency);

    add_parameter("Final time", final_time);

    add_parameter("Viscosity", viscosity);

    add_parameter("Nitsche penalty term", penalty_term);

    add_parameter("Initial fluid refinement",
                  initial_fluid_refinement,
                  "Initial mesh refinement used for the fluid domain Omega");

    add_parameter("Initial solid refinement",
                  initial_solid_refinement,
                  "Initial mesh refinement used for the solid domain Gamma");

    add_parameter(
      "Particle insertion refinement",
      particle_insertion_refinement,
      "Refinement of the volumetric mesh used to insert the particles");

    add_parameter(
      "Homogeneous Dirichlet boundary ids",
      homogeneous_dirichlet_ids,
      "Boundary Ids over which homogeneous Dirichlet boundary conditions are applied");

    // Next section is dedicated to the parameters used to create the
    // various grids. We will need three different triangulations: `Grid
    // one` is used to define the fluid domain, `Grid two` defines the
    // solid domain, and `Particle grid` is used to distribute some tracer
    // particles, that are advected with the velocity and only used as
    // passive tracers.
    enter_my_subsection(this->prm);
    this->prm.enter_subsection("Grid generation");
    this->prm.add_parameter("Grid one generator", name_of_fluid_grid);
    this->prm.add_parameter("Grid one generator arguments",
                            arguments_for_fluid_grid);

    this->prm.add_parameter("Grid two generator", name_of_solid_grid);
    this->prm.add_parameter("Grid two generator arguments",
                            arguments_for_solid_grid);

    this->prm.add_parameter("Particle grid generator", name_of_particle_grid);
    this->prm.add_parameter("Particle grid generator arguments",
                            arguments_for_particle_grid);
    this->prm.leave_subsection();

    leave_my_subsection(this->prm);



    enter_my_subsection(this->prm);
    this->prm.enter_subsection("Refinement and remeshing");
    this->prm.add_parameter("Refinement step frequency", refinement_frequency);
    this->prm.add_parameter("Refinement maximal level", max_level_refinement);
    this->prm.add_parameter("Refinement strategy",
                            refinement_strategy,
                            "",
                            Patterns::Selection("fixed_fraction|fixed_number"));
    this->prm.add_parameter("Refinement coarsening fraction",
                            coarsening_fraction);
    this->prm.add_parameter("Refinement fraction", refinement_fraction);
    this->prm.add_parameter("Maximum number of cells", max_cells);

    this->prm.leave_subsection();
    leave_my_subsection(this->prm);

    // correct the default dimension for the functions
    rhs.declare_parameters_call_back.connect([&]() {
      Functions::ParsedFunction<spacedim>::declare_parameters(this->prm,
                                                              spacedim + 1);
    });
    angular_velocity.declare_parameters_call_back.connect([&]() {
      Functions::ParsedFunction<spacedim>::declare_parameters(
        this->prm, spacedim == 3 ? spacedim : 1);
    });
  }

} // namespace Step70


// The remainder of the code, the main function, is standard.
int main(int argc, char *argv[])
{
  using namespace Step70;
  using namespace dealii;
  deallog.depth_console(1);
  try
    {
      Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

      StokesImmersedProblemParameters<2, 3> par;
      ParameterAcceptor::initialize("parameters_23.prm",
                                    "used_parameters_23.prm");

      StokesImmersedProblem<2, 3> problem(par);
      problem.run();
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
