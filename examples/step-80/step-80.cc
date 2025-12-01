/* ---------------------------------------------------------------------
 *
 * Copyright (C) 2025 by the deal.II authors
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
 * Authors: Luca Heltai, Bruno Blais, 2024
 */

// @sect3{Include files}
// We follow the step-70 documentation style: short section headers
// that guide readers through the main components. The headers below
// provide linear algebra back ends, particles, mappings, and CAD
// support needed for the immersed compressible solid/fluid coupling.
#include <deal.II/base/function.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/timer.h>

#include <deal.II/lac/generic_linear_algebra.h>
#include <deal.II/lac/linear_operator.h>
#include <deal.II/lac/linear_operator_tools.h>
#include <deal.II/lac/block_linear_operator.h>

#include <deal.II/lac/affine_constraints.templates.h>
#include <deal.II/fe/mapping_fe_field.h>

#include <boost/algorithm/string.hpp>
#include <deal.II/numerics/vector_tools_interpolate.h>

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
#include <deal.II/fe/fe_dgp.h>
#include <deal.II/fe/fe_dgp_nonparametric.h>
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
#include <deal.II/lac/block_linear_operator.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/error_estimator.h>
#include <deal.II/numerics/vector_tools.h>

#include <deal.II/particles/data_out.h>
#include <deal.II/particles/generators.h>
#include <deal.II/particles/particle_handler.h>
#include <deal.II/particles/utilities.h>

#include <deal.II/opencascade/manifold_lib.h>
#include <deal.II/opencascade/utilities.h>
#ifdef DEAL_II_WITH_OPENCASCADE
#  include <TopoDS.hxx>
#endif

#include <iostream>
#include <memory>

// @sect3{Run-time parameter handling}
// Helper functions and a ParameterAcceptor-derived structure collect all user
// options (dimensions, material properties, boundary data, mesh generators) and
// keep them synchronized with the parameter file, following the approach used
// in step-70.

namespace Step80
{
  using namespace dealii;

  template <typename VectorType, typename BlockVectorType>
  class StokesDLMALPreconditioner
  {
  public:
    using PayloadType = dealii::TrilinosWrappers::internal::
      LinearOperatorImplementation::TrilinosPayload;

    StokesDLMALPreconditioner(
      const LinearOperator<VectorType, VectorType> A11_inv_,
      const LinearOperator<VectorType, VectorType> A22_inv_,
      const LinearOperator<VectorType, VectorType> A12_,
      const LinearOperator<VectorType, VectorType> Bt_,
      const LinearOperator<VectorType, VectorType> Ct_,
      const LinearOperator<VectorType, VectorType> invW_,
      const LinearOperator<VectorType, VectorType> invMp_,
      const LinearOperator<VectorType, VectorType> M_,
      const double                                 gamma_fluid_,
      const double                                 gamma_solid_,
      const double                                 gamma_grad_div_)
    {
      // Aug_inv = Aug_inv_;
      A11_inv        = A11_inv_;
      A22_inv        = A22_inv_;
      A12            = A12_;
      Bt             = Bt_;
      Ct             = Ct_;
      invW           = invW_;
      invMp          = invMp_;
      M              = M_; // immersed
      gamma_fluid    = gamma_fluid_;
      gamma_solid    = gamma_solid_;
      gamma_grad_div = gamma_grad_div_;
    }

    void vmult(BlockVectorType &v, const BlockVectorType &u) const
    {
      v.block(0) = 0.;
      v.block(1) = 0.;
      v.block(2) = 0.;
      v.block(3) = 0.;

      v.block(3) = -gamma_grad_div * invMp * u.block(3);
      v.block(2) = -gamma_fluid * invW * u.block(2);
      v.block(1) = A22_inv * (u.block(1) - M * v.block(2));
      v.block(0) = A11_inv * (u.block(0) - A12 * v.block(1) - Ct * v.block(2) -
                              Bt * v.block(3));
    }

    // LinearOperator<Vector<double>> Aug_inv;
    LinearOperator<VectorType, VectorType> A11_inv;
    LinearOperator<VectorType, VectorType> A22_inv;
    LinearOperator<VectorType, VectorType> A12;
    LinearOperator<VectorType, VectorType> Ct;
    LinearOperator<VectorType, VectorType> Bt;
    LinearOperator<VectorType, VectorType> invW;
    LinearOperator<VectorType, VectorType> invMp;
    LinearOperator<VectorType, VectorType> M;
    double                                 gamma_fluid;
    double                                 gamma_solid;
    double                                 gamma_grad_div;
  };

  std::pair<unsigned int, unsigned int>
  get_dimension_and_spacedimension(const ParameterHandler &prm)
  {
    auto dim      = prm.get_integer("dimension");
    auto spacedim = prm.get_integer("space dimension");
    return {dim, spacedim};
  }

  // A free function to read the dimension and spacedimension from a parameter
  // file
  std::pair<unsigned int, unsigned int>
  get_dimension_and_spacedimension(const std::string &prm_file)
  {
    ParameterAcceptor::prm.declare_entry("dimension",
                                         "2",
                                         Patterns::Integer(2, 3));
    ParameterAcceptor::prm.declare_entry("space dimension",
                                         "2",
                                         Patterns::Integer(2, 3));
    // If reading of the input file fails, run by default in 2D-2D
    try
      {
        ParameterAcceptor::prm.parse_input(prm_file, "", true);
      }
    catch (std::exception &exc)
      {
        return {2, 2};
        throw;
      }
    return get_dimension_and_spacedimension(ParameterAcceptor::prm);
  }


  // @sect4{Parameter container}
  // NavierStokesImmersedProblemParameters bundles all inputs needed by the
  // coupled fluid/solid solver: polynomial degrees, material properties,
  // time-stepping controls, mesh generators, and right-hand sides. The class
  // mirrors the layout from step-70 but adds solid material constants and
  // multiple configuration fields for the immersed structure.

  template <int dim, int spacedim = dim>
  class NavierStokesImmersedProblemParameters : public ParameterAcceptor
  {
  public:
    NavierStokesImmersedProblemParameters();

    void set_time(const double &time) const;

    std::string output_directory = ".";

    unsigned int velocity_degree            = 2;
    unsigned int displacement_degree        = 1;
    unsigned int lagrange_multiplier_degree = 0;

    unsigned int number_of_time_steps = 501;
    double       final_time           = 1.0;

    unsigned int output_frequency = 1;

    unsigned int initial_fluid_refinement      = 5;
    unsigned int initial_solid_refinement      = 5;
    unsigned int particle_insertion_refinement = 3;

    unsigned int fluid_rtree_extraction_level = 1;

    double viscosity = 1.0;
    double density   = 1.0;

    double lame_mu     = 1.0;
    double lame_lambda = 1.0;

    std::list<types::boundary_id> dirichlet_ids{0};

    std::string name_of_fluid_grid           = "hyper_cube";
    std::string arguments_for_fluid_grid     = "-1: 1: false";
    std::string name_of_solid_grid           = "hyper_rectangle";
    std::string arguments_for_solid_grid     = spacedim == 2 ?
                                                 "-.5, -.1: .5, .1: false" :
                                                 "-.5, -.1, -.1: .5, .1, .1: false";
    std::string name_of_tracer_particle_grid = "hyper_ball";
    std::string arguments_for_tracer_particle_grid =
      spacedim == 2 ? "0.3, 0.3: 0.1: false" : "0.3, 0.3, 0.3 : 0.1: false";

    int          max_level_refinement = 8;
    int          min_level_refinement = 5;
    std::string  refinement_strategy  = "fixed_fraction";
    double       coarsening_fraction  = 0.3;
    double       refinement_fraction  = 0.3;
    unsigned int max_cells            = 20000;
    int          refinement_frequency = 5;

    double gamma_AL_background = 100;
    double gamma_AL_immersed   = 1e-2;

    unsigned int inner_max_iterations            = 100;
    double       inner_tolerance                 = 1e-14;
    unsigned int outer_max_iterations            = 1000;
    double       outer_tolerance                 = 1e-5;
    unsigned int inner_lagrangian_max_iterations = 1000;
    double       inner_lagrangian_tolerance      = 1e-2;

    using PrmFunction =
      ParameterAcceptorProxy<Functions::ParsedFunction<spacedim>>;

    // These functions depend on time, so we need to be able to set their
    // internal time. Make them mutable for that reason.
    mutable PrmFunction navier_stokes_rhs;
    mutable PrmFunction solid_rhs;
    mutable PrmFunction navier_stokes_bc;

    // These functions do not depend on time, so they don't need to be mutable
    PrmFunction navier_stokes_initial_conditions;
    PrmFunction solid_initial_displacement;
    PrmFunction solid_reference_configuration;
  };



  template <int dim, int spacedim>
  NavierStokesImmersedProblemParameters<dim, spacedim>::
    NavierStokesImmersedProblemParameters()
    : ParameterAcceptor("Navier-Stokes Immersed Problem/")
    , navier_stokes_rhs("Navier-Stokes right hand side", spacedim + 1)
    , solid_rhs("Solid right hand side", 2 * spacedim)
    , navier_stokes_bc("Navier-Stokes boundary conditions", spacedim + 1)
    , navier_stokes_initial_conditions("Navier-Stokes initial conditions",
                                       spacedim + 1)
    , solid_initial_displacement("Initial displacement and multiplier",
                                 2 * spacedim)
    , solid_reference_configuration("Reference configuration and multiplier",
                                    2 * spacedim)
  {
    add_parameter(
      "Velocity degree", velocity_degree, "", this->prm, Patterns::Integer(1));
    add_parameter("Displacement degree",
                  displacement_degree,
                  "",
                  this->prm,
                  Patterns::Integer(1));
    add_parameter("Lagrange multiplier degree",
                  lagrange_multiplier_degree,
                  "",
                  this->prm,
                  Patterns::Integer(0));

    add_parameter("Number of time steps", number_of_time_steps);
    add_parameter("Output frequency", output_frequency);

    add_parameter("Output directory", output_directory);

    add_parameter("Final time", final_time);

    add_parameter("Viscosity", viscosity);

    add_parameter("Density", density);

    add_parameter("Lame mu", lame_mu);

    add_parameter("Lame lambda", lame_lambda);

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
      "Tracer particle insertion refinement",
      particle_insertion_refinement,
      "Refinement of the volumetric mesh used to insert the particles");

    add_parameter(
      "Dirichlet boundary ids",
      dirichlet_ids,
      "Boundary Ids over which Dirichlet boundary conditions are applied");

    enter_subsection("Grid generation");
    {
      add_parameter("Fluid grid generator", name_of_fluid_grid);
      add_parameter("Fluid grid generator arguments", arguments_for_fluid_grid);

      add_parameter("Solid grid generator", name_of_solid_grid);
      add_parameter("Solid grid generator arguments", arguments_for_solid_grid);

      add_parameter("Tracer particle grid generator",
                    name_of_tracer_particle_grid);
      add_parameter("Tracer particle grid generator arguments",
                    arguments_for_tracer_particle_grid);
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

    enter_subsection("Solver parameters");
    {
      add_parameter("Gamma AL background", gamma_AL_background);
      add_parameter("Gamma AL immersed", gamma_AL_immersed);
      add_parameter("Inner solver maximum iterations", inner_max_iterations);
      add_parameter("Inner solver tolerance", inner_tolerance);
      add_parameter("Outer solver maximum iterations", outer_max_iterations);
      add_parameter("Outer solver tolerance", outer_tolerance);
      add_parameter("Inner Lagrangian solver maximum iterations",
                    inner_lagrangian_max_iterations);
      add_parameter("Inner Lagrangian solver tolerance",
                    inner_lagrangian_tolerance);
    }
    leave_subsection();

    // Make sure that dim and spacedim are actually declared also in the
    // application
    declare_parameters_call_back.connect([&]() {
      this->leave_my_subsection(this->prm);
      this->prm.declare_entry("dimension", "2", Patterns::Integer(2, 3));
      this->prm.declare_entry("space dimension", "2", Patterns::Integer(2, 3));
      this->enter_my_subsection(this->prm);
    });

    // And make sure we check that the dimension and space dimension in the file
    // match the ones of the application
    parse_parameters_call_back.connect([&]() {
      this->leave_my_subsection(this->prm);
      const auto [dim_, spacedim_] =
        get_dimension_and_spacedimension(this->prm);
      AssertThrow(dim_ == dim && spacedim_ == spacedim,
                  ExcMessage(
                    "The dimension and space dimension in the parameter "
                    "file do not match the ones of the running application."
                    " This should not happen: aborting."));
      this->enter_my_subsection(this->prm);
    });
  }



  template <int dim, int spacedim>
  inline void NavierStokesImmersedProblemParameters<dim, spacedim>::set_time(
    const double &time) const
  {
    navier_stokes_rhs.set_time(time);
    solid_rhs.set_time(time);
    navier_stokes_bc.set_time(time);
  }



  // @sect3{The NavierStokesImmersedProblem class declaration}
  // The driver holds distributed triangulations for fluid and solid, mapping
  // objects to deform the solid mesh, particle handlers for coupling quadrature
  // transfer, and block matrices/vectors for the coupled Navier-Stokes /
  // compressible elasticity system with a distributed Lagrange multiplier.
  template <int dim, int spacedim = dim>
  class NavierStokesImmersedProblem
  {
  public:
    NavierStokesImmersedProblem(
      const NavierStokesImmersedProblemParameters<dim, spacedim> &par);

    void run();

  private:
    void make_grid();

    double compute_time_step() const;

    void setup_solid_particles();

    void initial_setup();
    void setup_dofs();
    void interpolate_initial_conditions();
    void setup_coupling();

    void assemble_navier_stokes_system(const double &time_step);
    void assemble_navier_stokes_rhs(const double &time_step);

    void assemble_elasticity_system(const double &time_step);
    void assemble_elasticity_rhs(const double &time_step);

    void assemble_coupling_sparsity(BlockDynamicSparsityPattern &dsp);
    void assemble_coupling();

    void solve(const double time_step);

    void refine_and_transfer();

    void output_results(const unsigned int cycle, const double time) const;
    void output_particles(const Particles::ParticleHandler<spacedim> &particles,
                          const std::string                          &fprefix,
                          const unsigned int                          iter,
                          const double time) const;

    void update_particle_positions();

    const NavierStokesImmersedProblemParameters<dim, spacedim> &par;

    MPI_Comm mpi_communicator;

    ConditionalOStream pcout;

    mutable TimerOutput computing_timer;

    parallel::distributed::Triangulation<spacedim>      fluid_tria;
    parallel::distributed::Triangulation<dim, spacedim> solid_tria;

    std::unique_ptr<FiniteElement<spacedim>>      fluid_fe;
    std::unique_ptr<FiniteElement<dim, spacedim>> solid_fe;

    DoFHandler<spacedim>      fluid_dh;
    DoFHandler<dim, spacedim> solid_dh;

    std::unique_ptr<MappingFEField<dim, spacedim, LA::MPI::BlockVector>>
      solid_mapping;

    std::vector<types::global_dof_index> fluid_dofs_per_block;
    std::vector<types::global_dof_index> solid_dofs_per_block;

    std::vector<IndexSet> fluid_owned_dofs;
    std::vector<IndexSet> solid_owned_dofs;

    std::vector<IndexSet> fluid_relevant_dofs;
    std::vector<IndexSet> solid_relevant_dofs;

    AffineConstraints<double> fluid_constraints;
    AffineConstraints<double> solid_constraints;

    LA::MPI::BlockSparseMatrix fluid_matrix; // velocity and pressure
    LA::MPI::BlockSparseMatrix fluid_mass_matrix;
    LA::MPI::BlockSparseMatrix fluid_preconditioner;

    LA::MPI::BlockSparseMatrix
      solid_matrix; // displacement and Lagrange multiplier
    LA::MPI::BlockSparseMatrix
      coupling_interpolation_matrix; // between displacement and velocity

    LA::MPI::BlockVector fluid_solution;
    LA::MPI::BlockVector fluid_locally_relevant_solution;
    LA::MPI::BlockVector fluid_locally_relevant_solution_old;
    LA::MPI::BlockVector fluid_system_rhs;
    LA::MPI::BlockVector fluid_dual_of_constant_pressure;
    LA::MPI::BlockVector fluid_constant_pressure;

    LA::MPI::BlockVector solid_solution;
    LA::MPI::BlockVector solid_locally_relevant_solution;
    LA::MPI::BlockVector solid_locally_relevant_solution_old;
    LA::MPI::BlockVector solid_system_rhs;

    LA::MPI::BlockVector solid_reference_configuration;
    LA::MPI::BlockVector solid_current_position;

    Particles::ParticleHandler<spacedim> tracer_particle_handler;
    Particles::ParticleHandler<spacedim> solid_particle_handler;

    IndexSet locally_owned_tracer_particle_coordinates;
    IndexSet locally_relevant_tracer_particle_coordinates;

    LA::MPI::Vector tracer_particle_velocity;
    LA::MPI::Vector relevant_tracer_particle_displacements;

    std::vector<std::vector<BoundingBox<spacedim>>> global_fluid_bounding_boxes;

    FEValuesExtractors::Vector velocity;
    FEValuesExtractors::Scalar pressure;

    FEValuesExtractors::Vector displacement;
    FEValuesExtractors::Vector lagrange_multiplier;
  };



  // @sect3{The NavierStokesImmersedProblem class implementation}
  // The following functions are organized in the same spirit as step-70:
  // construction and mesh initialization, particle setup, DoF initialization,
  // assembly of fluid/solid/coupling blocks, solver, refinement, output, and
  // the main time loop.


  // @sect4{Construction and mesh initialization}
  template <int dim, int spacedim>
  NavierStokesImmersedProblem<dim, spacedim>::NavierStokesImmersedProblem(
    const NavierStokesImmersedProblemParameters<dim, spacedim> &par)
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
    , velocity(0)
    , pressure(spacedim)
    , displacement(0)
    , lagrange_multiplier(spacedim)
  {}


  template <int dim, int spacedim>
  void read_grid_and_cad_files(const std::string &grid_file_name,
                               const std::string &ids_and_cad_file_names,
                               Triangulation<dim, spacedim> &tria)
  {
    GridIn<dim, spacedim> grid_in;
    grid_in.attach_triangulation(tria);
    grid_in.read(grid_file_name);

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

        const auto n_elements = OpenCASCADE::count_elements(shape);
        if ((std::get<0>(n_elements) == 0))
          tria.set_manifold(
            manifold_id,
            OpenCASCADE::ArclengthProjectionLineManifold<dim, spacedim>(shape));
        else if constexpr (spacedim == 3)
          {
            tria.set_manifold(
              manifold_id,
              OpenCASCADE::NormalToMeshProjectionManifold<dim, 3>(shape));
          }
        else
          tria.set_manifold(manifold_id,
                            OpenCASCADE::NURBSPatchManifold<dim, spacedim>(
                              TopoDS::Face(shape)));
      }
#else
    (void)ids_and_cad_file_names;
    AssertThrow(false, ExcNotImplemented("Generation of the grid failed."));
#endif
  }



  // @sect4{Mesh generation}
  template <int dim, int spacedim>
  void NavierStokesImmersedProblem<dim, spacedim>::make_grid()
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


  // @sect4{Particle initialization}
  template <int dim, int spacedim>
  void NavierStokesImmersedProblem<dim, spacedim>::setup_solid_particles()
  {
    const unsigned int n_properties = 0;
    solid_particle_handler.initialize(fluid_tria,
                                      StaticMappingQ1<spacedim>::mapping,
                                      n_properties);

    std::vector<BoundingBox<spacedim>> all_boxes;
    all_boxes.reserve(fluid_tria.n_locally_owned_active_cells());
    for (const auto &cell : fluid_tria.active_cell_iterators())
      if (cell->is_locally_owned())
        all_boxes.emplace_back(cell->bounding_box());

    const auto tree = pack_rtree(all_boxes);
    const auto local_boxes =
      extract_rtree_level(tree, par.fluid_rtree_extraction_level);

    global_fluid_bounding_boxes =
      Utilities::MPI::all_gather(mpi_communicator, local_boxes);

    std::vector<bool> components(2 * spacedim, false);
    components[0] = true;

    Particles::Generators::dof_support_points(solid_dh,
                                              global_fluid_bounding_boxes,
                                              solid_particle_handler,
                                              *solid_mapping,
                                              ComponentMask(components));

    pcout << "   Number of particles (solid support points): "
          << solid_particle_handler.n_global_particles() << " ("
          << solid_particle_handler.n_global_particles() * spacedim
          << " displacement dofs)" << std::endl;

    AssertDimension(solid_particle_handler.n_global_particles() * spacedim,
                    solid_solution.block(0).size());
  }



  // @sect4{Finite element and DoF initialization}
  template <int dim, int spacedim>
  void NavierStokesImmersedProblem<dim, spacedim>::initial_setup()
  {
    TimerOutput::Scope t(computing_timer, "Initial setup");

    fluid_fe =
      std::make_unique<FESystem<spacedim>>(FE_Q<spacedim>(par.velocity_degree),
                                           spacedim,
                                           FE_Q<spacedim>(
                                             //  FE_DGPNonparametric<spacedim>(
                                             par.velocity_degree - 1),
                                           1);

    // Solid displacement and Lagrange multiplier (same degree for both)
    solid_fe = std::make_unique<FESystem<spacedim>>(
      FE_Q<dim, spacedim>(par.displacement_degree),
      spacedim,
      par.lagrange_multiplier_degree == 0 ?
        *FE_DGQ<dim, spacedim>(0).clone() :
        *FE_Q<dim, spacedim>(par.lagrange_multiplier_degree).clone(),
      spacedim);
  }


  template <int dim, int spacedim>
  void NavierStokesImmersedProblem<dim, spacedim>::setup_dofs()
  {
    TimerOutput::Scope t(computing_timer, "Setup dofs");

    fluid_dh.distribute_dofs(*fluid_fe);

    std::vector<unsigned int> navier_stokes_sub_blocks(spacedim + 1, 0);
    navier_stokes_sub_blocks[spacedim] = 1;
    DoFRenumbering::component_wise(fluid_dh, navier_stokes_sub_blocks);

    fluid_dofs_per_block =
      DoFTools::count_dofs_per_fe_block(fluid_dh, navier_stokes_sub_blocks);

    const unsigned int n_u = fluid_dofs_per_block[0],
                       n_p = fluid_dofs_per_block[1];

    pcout << "   Number of degrees of freedom for Navie-Stokes equation: "
          << fluid_dh.n_dofs() << " (" << n_u << '+' << n_p << ')' << std::endl;

    fluid_owned_dofs.resize(2);
    fluid_owned_dofs[0] = fluid_dh.locally_owned_dofs().get_view(0, n_u);
    fluid_owned_dofs[1] =
      fluid_dh.locally_owned_dofs().get_view(n_u, n_u + n_p);

    const IndexSet locally_relevant_fluid_dofs =
      DoFTools::extract_locally_relevant_dofs(fluid_dh);
    fluid_relevant_dofs.resize(2);
    fluid_relevant_dofs[0] = locally_relevant_fluid_dofs.get_view(0, n_u);
    fluid_relevant_dofs[1] =
      locally_relevant_fluid_dofs.get_view(n_u, n_u + n_p);

    {
      fluid_constraints.reinit(locally_relevant_fluid_dofs,
                               locally_relevant_fluid_dofs);

      DoFTools::make_hanging_node_constraints(fluid_dh, fluid_constraints);
      for (const auto id : par.dirichlet_ids)
        VectorTools::interpolate_boundary_values(fluid_dh,
                                                 id,
                                                 par.navier_stokes_bc,
                                                 fluid_constraints,
                                                 fluid_fe->component_mask(
                                                   velocity));

      fluid_constraints.close();
    }

    auto locally_owned_fluid_dofs_per_processor =
      Utilities::MPI::all_gather(mpi_communicator,
                                 fluid_dh.locally_owned_dofs());
    {
      fluid_matrix.clear();

      Table<2, DoFTools::Coupling> coupling(spacedim + 1, spacedim + 1);
      for (unsigned int c = 0; c < spacedim + 1; ++c)
        for (unsigned int d = 0; d < spacedim + 1; ++d)
          if (c == spacedim && d == spacedim)
            coupling[c][d] = DoFTools::none;
          else if (c == spacedim || d == spacedim || c == d)
            coupling[c][d] = DoFTools::always;
          else
            coupling[c][d] = DoFTools::none;

      BlockDynamicSparsityPattern dsp(fluid_dofs_per_block,
                                      fluid_dofs_per_block);

      DoFTools::make_sparsity_pattern(
        fluid_dh, coupling, dsp, fluid_constraints, false);

      SparsityTools::distribute_sparsity_pattern(
        dsp,
        locally_owned_fluid_dofs_per_processor,
        mpi_communicator,
        locally_relevant_fluid_dofs);

      fluid_matrix.reinit(fluid_owned_dofs, dsp, mpi_communicator);
    }

    {
      fluid_preconditioner.clear();

      Table<2, DoFTools::Coupling> coupling(spacedim + 1, spacedim + 1);
      for (unsigned int c = 0; c < spacedim + 1; ++c)
        for (unsigned int d = 0; d < spacedim + 1; ++d)
          if (c == spacedim && d == spacedim)
            coupling[c][d] = DoFTools::always;
          else
            coupling[c][d] = DoFTools::none;

      BlockDynamicSparsityPattern dsp(fluid_dofs_per_block,
                                      fluid_dofs_per_block);

      DoFTools::make_sparsity_pattern(
        fluid_dh, coupling, dsp, fluid_constraints, false);
      SparsityTools::distribute_sparsity_pattern(
        dsp,
        locally_owned_fluid_dofs_per_processor,
        mpi_communicator,
        locally_relevant_fluid_dofs);
      fluid_preconditioner.reinit(fluid_owned_dofs, dsp, mpi_communicator);
    }

    {
      fluid_mass_matrix.clear();

      Table<2, DoFTools::Coupling> coupling(spacedim + 1, spacedim + 1);
      for (unsigned int c = 0; c < spacedim + 1; ++c)
        for (unsigned int d = 0; d < spacedim + 1; ++d)
          if (c == spacedim && d == spacedim)
            coupling[c][d] = DoFTools::none;
          else if (c == spacedim || d == spacedim || c == d)
            coupling[c][d] = DoFTools::always;
          else
            coupling[c][d] = DoFTools::none;

      BlockDynamicSparsityPattern dsp(fluid_dofs_per_block,
                                      fluid_dofs_per_block);

      DoFTools::make_sparsity_pattern(
        fluid_dh, coupling, dsp, fluid_constraints, false);

      SparsityTools::distribute_sparsity_pattern(
        dsp,
        locally_owned_fluid_dofs_per_processor,
        mpi_communicator,
        locally_relevant_fluid_dofs);

      fluid_mass_matrix.reinit(fluid_owned_dofs, dsp, mpi_communicator);
    }



    fluid_locally_relevant_solution.reinit(fluid_owned_dofs,
                                           fluid_relevant_dofs,
                                           mpi_communicator);
    fluid_locally_relevant_solution_old.reinit(fluid_owned_dofs,
                                               fluid_relevant_dofs,
                                               mpi_communicator);
    fluid_system_rhs.reinit(fluid_owned_dofs, mpi_communicator);
    fluid_solution.reinit(fluid_owned_dofs, mpi_communicator);
    fluid_dual_of_constant_pressure.reinit(fluid_owned_dofs, mpi_communicator);

    VectorTools::create_right_hand_side(
      fluid_dh,
      QGauss<spacedim>(fluid_fe->degree + 1),
      ComponentSelectFunction<spacedim>(spacedim, 1.0, spacedim + 1),
      fluid_dual_of_constant_pressure,
      fluid_constraints);


    // Setup system for the solid mechanics component
    {
      solid_dh.distribute_dofs(*solid_fe);

      std::vector<unsigned int> solid_sub_blocks(2 * spacedim, 0);
      for (unsigned int d = spacedim; d < 2 * spacedim; ++d)
        solid_sub_blocks[d] = 1;
      DoFRenumbering::component_wise(solid_dh, solid_sub_blocks);
      solid_dofs_per_block =
        DoFTools::count_dofs_per_fe_block(solid_dh, solid_sub_blocks);

      const unsigned int n_disp = solid_dofs_per_block[0],
                         n_lag  = solid_dofs_per_block[1];

      pcout << "   Number of degrees of freedom for solid mechanics: "
            << solid_dh.n_dofs() << " (" << n_disp << '+' << n_lag << ')'
            << std::endl;

      solid_owned_dofs.resize(2);
      solid_owned_dofs[0] = solid_dh.locally_owned_dofs().get_view(0, n_disp);
      solid_owned_dofs[1] =
        solid_dh.locally_owned_dofs().get_view(n_disp, n_disp + n_lag);

      const IndexSet locally_relevant_solid_dofs =
        DoFTools::extract_locally_relevant_dofs(solid_dh);
      solid_relevant_dofs.resize(2);
      solid_relevant_dofs[0] = locally_relevant_solid_dofs.get_view(0, n_disp);
      solid_relevant_dofs[1] =
        locally_relevant_solid_dofs.get_view(n_disp, n_disp + n_lag);

      solid_constraints.reinit(locally_relevant_solid_dofs,
                               locally_relevant_solid_dofs);


      solid_locally_relevant_solution.reinit(solid_owned_dofs,
                                             solid_relevant_dofs,
                                             mpi_communicator);

      solid_locally_relevant_solution_old.reinit(solid_owned_dofs,
                                                 solid_relevant_dofs,
                                                 mpi_communicator);

      solid_system_rhs.reinit(solid_owned_dofs, mpi_communicator);
      solid_solution.reinit(solid_owned_dofs, mpi_communicator);

      DoFTools::make_hanging_node_constraints(solid_dh, solid_constraints);
      solid_constraints.close();

      BlockDynamicSparsityPattern dsp(solid_dofs_per_block,
                                      solid_dofs_per_block);

      DoFTools::make_sparsity_pattern(solid_dh, dsp, solid_constraints, false);

      auto locally_owned_solid_dofs_per_processor =
        Utilities::MPI::all_gather(mpi_communicator,
                                   solid_dh.locally_owned_dofs());

      SparsityTools::distribute_sparsity_pattern(
        dsp,
        locally_owned_solid_dofs_per_processor,
        mpi_communicator,
        locally_relevant_solid_dofs);

      solid_matrix.reinit(solid_owned_dofs, dsp, mpi_communicator);

      solid_reference_configuration.reinit(solid_locally_relevant_solution);
      solid_current_position.reinit(solid_locally_relevant_solution);
    }
  }



  template <int dim, int spacedim>
  void
  NavierStokesImmersedProblem<dim, spacedim>::interpolate_initial_conditions()
  {
    VectorTools::interpolate(fluid_dh,
                             par.navier_stokes_initial_conditions,
                             fluid_solution,
                             ComponentMask(fluid_fe->component_mask(velocity)));
    fluid_locally_relevant_solution_old = fluid_solution;

    VectorTools::interpolate(solid_dh,
                             par.solid_reference_configuration,
                             solid_solution,
                             ComponentMask(
                               solid_fe->component_mask(displacement)));
    solid_reference_configuration = solid_solution;

    VectorTools::interpolate(solid_dh,
                             par.solid_initial_displacement,
                             solid_solution,
                             ComponentMask(
                               solid_fe->component_mask(displacement)));
    solid_locally_relevant_solution_old = solid_solution;

    solid_current_position = solid_reference_configuration;
    solid_current_position += solid_locally_relevant_solution_old;
    solid_mapping =
      std::make_unique<MappingFEField<dim, spacedim, LA::MPI::BlockVector>>(
        solid_dh,
        solid_current_position,
        ComponentMask(solid_fe->component_mask(displacement)));
  }



  template <int dim, int spacedim>
  void NavierStokesImmersedProblem<dim, spacedim>::setup_coupling()
  {
    BlockDynamicSparsityPattern dsp(fluid_dofs_per_block, solid_dofs_per_block);
    assemble_coupling_sparsity(dsp);

    coupling_interpolation_matrix.clear();
    coupling_interpolation_matrix.reinit(fluid_dofs_per_block.size(),
                                         solid_dofs_per_block.size());
    for (unsigned int i = 0; i < fluid_dofs_per_block.size(); ++i)
      for (unsigned int j = 0; j < solid_dofs_per_block.size(); ++j)
        coupling_interpolation_matrix.block(i, j).reinit(fluid_owned_dofs[i],
                                                         solid_owned_dofs[j],
                                                         dsp.block(i, j),
                                                         mpi_communicator,
                                                         true);
    coupling_interpolation_matrix.collect_sizes();
  }



  // @sect4{Assembly of fluid, solid, and coupling blocks}
  template <int dim, int spacedim>
  void
  NavierStokesImmersedProblem<dim, spacedim>::assemble_navier_stokes_system(
    const double &time_step)
  {
    fluid_matrix = 0;

    TimerOutput::Scope t(computing_timer, "Assemble Navier-Stokes terms");

    QGauss<spacedim>   quadrature_formula(fluid_fe->degree + 1);
    FEValues<spacedim> fe_values(*fluid_fe,
                                 quadrature_formula,
                                 update_values | update_gradients |
                                   update_quadrature_points |
                                   update_JxW_values);

    const unsigned int dofs_per_cell = fluid_fe->n_dofs_per_cell();
    const unsigned int n_q_points    = quadrature_formula.size();

    FullMatrix<double> cell_fluid_matrix(dofs_per_cell, dofs_per_cell);
    FullMatrix<double> cell_fluid_preconditioner(dofs_per_cell, dofs_per_cell);
    FullMatrix<double> cell_fluid_mass_matrix(dofs_per_cell, dofs_per_cell);

    std::vector<Tensor<1, spacedim>> phi_u(dofs_per_cell);
    std::vector<Tensor<2, spacedim>> grad_phi_u(dofs_per_cell);
    std::vector<double>              div_phi_u(dofs_per_cell);
    std::vector<double>              phi_p(dofs_per_cell);
    std::vector<Tensor<1, spacedim>> u(n_q_points);
    std::vector<Tensor<2, spacedim>> grad_u(n_q_points);

    std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

    for (const auto &cell : fluid_dh.active_cell_iterators())
      if (cell->is_locally_owned())
        {
          cell_fluid_matrix         = 0;
          cell_fluid_preconditioner = 0;
          cell_fluid_mass_matrix    = 0;

          fe_values.reinit(cell);
          for (unsigned int q = 0; q < n_q_points; ++q)
            {
              for (unsigned int k = 0; k < dofs_per_cell; ++k)
                {
                  phi_u[k]      = fe_values[velocity].value(k, q);
                  grad_phi_u[k] = fe_values[velocity].gradient(k, q);
                  div_phi_u[k]  = fe_values[velocity].divergence(k, q);

                  phi_p[k] = fe_values[pressure].value(k, q);
                }

              for (unsigned int i = 0; i < dofs_per_cell; ++i)
                {
                  for (unsigned int j = 0; j < dofs_per_cell; ++j)
                    {
                      cell_fluid_matrix(i, j) +=
                        ((par.density / time_step) * phi_u[i] * phi_u[j] +
                         par.viscosity *
                           scalar_product(grad_phi_u[i], grad_phi_u[j]) -
                         div_phi_u[i] * phi_p[j] - phi_p[i] * div_phi_u[j]) *
                        fe_values.JxW(q);

                      cell_fluid_preconditioner(i, j) +=
                        phi_p[i] * phi_p[j] * fe_values.JxW(q);

                      cell_fluid_mass_matrix(i, j) +=
                        (par.density / time_step) * phi_u[i] * phi_u[j] *
                        fe_values.JxW(q);
                    }
                }
            }


          cell->get_dof_indices(local_dof_indices);
          fluid_constraints.distribute_local_to_global(cell_fluid_matrix,
                                                       local_dof_indices,
                                                       fluid_matrix);

          fluid_constraints.distribute_local_to_global(
            cell_fluid_preconditioner, local_dof_indices, fluid_preconditioner);

          fluid_constraints.distribute_local_to_global(cell_fluid_mass_matrix,
                                                       local_dof_indices,
                                                       fluid_mass_matrix);
        }

    fluid_matrix.compress(VectorOperation::add);
    fluid_preconditioner.compress(VectorOperation::add);
    fluid_mass_matrix.compress(VectorOperation::add);
  }


  template <int dim, int spacedim>
  void NavierStokesImmersedProblem<dim, spacedim>::assemble_navier_stokes_rhs(
    const double &time_step)
  {
    fluid_system_rhs = 0;

    TimerOutput::Scope t(computing_timer, "Assemble Navier-Stokes rhs terms");

    QGauss<spacedim>   quadrature_formula(fluid_fe->degree + 1);
    FEValues<spacedim> fe_values(*fluid_fe,
                                 quadrature_formula,
                                 update_values | update_gradients |
                                   update_quadrature_points |
                                   update_JxW_values);

    const unsigned int dofs_per_cell = fluid_fe->n_dofs_per_cell();
    const unsigned int n_q_points    = quadrature_formula.size();

    FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);
    Vector<double>     cell_rhs(dofs_per_cell);

    std::vector<Vector<double>> rhs_values(n_q_points,
                                           Vector<double>(spacedim + 1));

    std::vector<Tensor<1, spacedim>> phi_u(dofs_per_cell);
    std::vector<Tensor<2, spacedim>> grad_phi_u(dofs_per_cell);
    std::vector<double>              div_phi_u(dofs_per_cell);
    std::vector<double>              phi_p(dofs_per_cell);
    std::vector<Tensor<1, spacedim>> u(n_q_points);
    std::vector<Tensor<2, spacedim>> grad_u(n_q_points);

    std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

    for (const auto &cell : fluid_dh.active_cell_iterators())
      if (cell->is_locally_owned())
        {
          cell_matrix = 0;
          cell_rhs    = 0;

          fe_values.reinit(cell);
          par.navier_stokes_rhs.vector_value_list(
            fe_values.get_quadrature_points(), rhs_values);
          fe_values[velocity].get_function_values(
            fluid_locally_relevant_solution_old, u);
          fe_values[velocity].get_function_gradients(
            fluid_locally_relevant_solution_old, grad_u);
          for (unsigned int q = 0; q < n_q_points; ++q)
            {
              for (unsigned int k = 0; k < dofs_per_cell; ++k)
                {
                  phi_u[k]      = fe_values[velocity].value(k, q);
                  grad_phi_u[k] = fe_values[velocity].gradient(k, q);
                  div_phi_u[k]  = fe_values[velocity].divergence(k, q);

                  phi_p[k] = fe_values[pressure].value(k, q);
                }

              for (unsigned int i = 0; i < dofs_per_cell; ++i)
                {
                  if (cell->at_boundary())
                    {
                      for (unsigned int j = 0; j < dofs_per_cell; ++j)
                        {
                          cell_matrix(i, j) +=
                            (par.density / time_step * phi_u[i] * phi_u[j] +
                             par.viscosity *
                               scalar_product(grad_phi_u[i], grad_phi_u[j]) -
                             div_phi_u[i] * phi_p[j] -
                             phi_p[i] * div_phi_u[j]) *
                            fe_values.JxW(q);
                        }
                    }

                  const unsigned int component_i =
                    fluid_fe->system_to_component_index(i).first;
                  cell_rhs(i) +=
                    (par.density *
                     (fe_values.shape_value(i, q) * rhs_values[q](component_i) +
                      (u[q] / time_step - grad_u[q] * u[q]) * phi_u[i])) *
                    fe_values.JxW(q);
                }
            }


          cell->get_dof_indices(local_dof_indices);
          fluid_constraints.distribute_local_to_global(cell_rhs,
                                                       local_dof_indices,
                                                       fluid_system_rhs,
                                                       cell_matrix);
        }

    fluid_system_rhs.compress(VectorOperation::add);
  }


  template <int dim, int spacedim>
  void NavierStokesImmersedProblem<dim, spacedim>::assemble_elasticity_system(
    const double &time_step)
  {
    solid_matrix = 0;

    TimerOutput::Scope t(computing_timer, "Assemble Elasticity terms");

    QGauss<spacedim>   quadrature_formula(solid_fe->degree + 1);
    FEValues<spacedim> fe_values(*solid_fe,
                                 quadrature_formula,
                                 update_values | update_gradients |
                                   update_quadrature_points |
                                   update_JxW_values);

    const unsigned int dofs_per_cell = solid_fe->n_dofs_per_cell();
    const unsigned int n_q_points    = quadrature_formula.size();

    FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);

    std::vector<SymmetricTensor<2, spacedim>> grad_eps_phi_w(dofs_per_cell);
    std::vector<Tensor<1, spacedim>>          phi_w(dofs_per_cell);
    std::vector<double>                       div_phi_w(dofs_per_cell);
    std::vector<Tensor<1, spacedim>>          phi_lagrange(dofs_per_cell);

    std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

    for (const auto &cell : solid_dh.active_cell_iterators())
      if (cell->is_locally_owned())
        {
          cell_matrix = 0;

          fe_values.reinit(cell);

          for (unsigned int q = 0; q < n_q_points; ++q)
            {
              for (unsigned int k = 0; k < dofs_per_cell; ++k)
                {
                  grad_eps_phi_w[k] =
                    fe_values[displacement].symmetric_gradient(k, q);
                  div_phi_w[k]    = fe_values[displacement].divergence(k, q);
                  phi_w[k]        = fe_values[displacement].value(k, q);
                  phi_lagrange[k] = fe_values[lagrange_multiplier].value(k, q);
                }

              for (unsigned int i = 0; i < dofs_per_cell; ++i)
                {
                  for (unsigned int j = 0; j < dofs_per_cell; ++j)
                    {
                      cell_matrix(i, j) +=
                        (2 * par.lame_mu *
                           scalar_product(grad_eps_phi_w[i],
                                          grad_eps_phi_w[j]) +
                         par.lame_lambda * div_phi_w[i] * div_phi_w[j] -
                         // lagrange * disp
                         (phi_lagrange[i] * phi_w[j] / time_step) -
                         // disp * lagrange
                         (phi_w[i] * phi_lagrange[j] / time_step) +
                         // lagr * lagr
                         phi_lagrange[i] * phi_lagrange[j]) *
                        // JxW
                        fe_values.JxW(q);
                    }
                }
            }

          cell->get_dof_indices(local_dof_indices);
          solid_constraints.distribute_local_to_global(cell_matrix,
                                                       local_dof_indices,
                                                       solid_matrix);
        }
    solid_matrix.compress(VectorOperation::add);
  }



  template <int dim, int spacedim>
  void NavierStokesImmersedProblem<dim, spacedim>::assemble_elasticity_rhs(
    const double &time_step)
  {
    solid_system_rhs = 0;

    TimerOutput::Scope t(computing_timer, "Assemble Elastic rhs terms");

    QGauss<spacedim>   quadrature_formula(solid_fe->degree + 1);
    FEValues<spacedim> fe_values_rhs(*solid_fe,
                                     quadrature_formula,
                                     update_values | update_quadrature_points |
                                       update_JxW_values);

    FEValues<spacedim> fe_values_matrix(*solid_fe,
                                        quadrature_formula,
                                        update_values | update_gradients |
                                          update_JxW_values);

    const unsigned int dofs_per_cell = solid_fe->n_dofs_per_cell();
    const unsigned int n_q_points    = quadrature_formula.size();

    FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);
    Vector<double>     cell_rhs(dofs_per_cell);

    std::vector<SymmetricTensor<2, spacedim>> grad_eps_phi_w(dofs_per_cell);
    std::vector<Tensor<1, spacedim>>          phi_w(dofs_per_cell);
    std::vector<double>                       div_phi_w(dofs_per_cell);
    std::vector<Tensor<1, spacedim>>          phi_lagrange(dofs_per_cell);

    std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

    std::vector<Tensor<1, spacedim>> w_old(n_q_points);
    std::vector<Vector<double>>      solid_rhs_values(n_q_points,
                                                 Vector<double>(2 * spacedim));

    for (const auto &cell : solid_dh.active_cell_iterators())
      if (cell->is_locally_owned())
        {
          fe_values_rhs.reinit(cell);

          fe_values_rhs[displacement].get_function_values(
            solid_locally_relevant_solution_old, w_old);

          par.solid_rhs.vector_value_list(fe_values_rhs.get_quadrature_points(),
                                          solid_rhs_values);

          cell_rhs = 0;
          for (unsigned int q = 0; q < n_q_points; ++q)
            {
              for (unsigned int i = 0; i < dofs_per_cell; ++i)
                {
                  const auto &phi_w =
                    fe_values_rhs[lagrange_multiplier].value(i, q);

                  const auto comp_i =
                    solid_fe->system_to_component_index(i).first;

                  cell_rhs(i) += (-w_old[q] / time_step * phi_w +
                                  solid_rhs_values[q](comp_i) *
                                    fe_values_rhs.shape_value(i, q)) *
                                 fe_values_rhs.JxW(q);
                }
            }
          cell->get_dof_indices(local_dof_indices);

          if (cell->at_boundary())
            {
              fe_values_matrix.reinit(cell);
              cell_matrix = 0;
              for (unsigned int q = 0; q < n_q_points; ++q)
                {
                  for (unsigned int k = 0; k < dofs_per_cell; ++k)
                    {
                      grad_eps_phi_w[k] =
                        fe_values_matrix[displacement].symmetric_gradient(k, q);
                      div_phi_w[k] =
                        fe_values_matrix[displacement].divergence(k, q);
                      phi_w[k] = fe_values_matrix[displacement].value(k, q);
                      phi_lagrange[k] =
                        fe_values_matrix[lagrange_multiplier].value(k, q);
                    }

                  for (unsigned int i = 0; i < dofs_per_cell; ++i)
                    {
                      for (unsigned int j = 0; j < dofs_per_cell; ++j)
                        {
                          cell_matrix(i, j) +=
                            (2 * par.lame_mu *
                               scalar_product(grad_eps_phi_w[i],
                                              grad_eps_phi_w[j]) +
                             par.lame_lambda * div_phi_w[i] * div_phi_w[j] -
                             // lagrange * disp
                             (phi_lagrange[i] * phi_w[j] / time_step) -
                             // disp * lagrange
                             (phi_w[i] * phi_lagrange[j] / time_step) +
                             // lagr * lagr
                             phi_lagrange[i] * phi_lagrange[j]) *
                            // JxW
                            fe_values_matrix.JxW(q);
                        }
                    }
                }
              solid_constraints.distribute_local_to_global(cell_rhs,
                                                           local_dof_indices,
                                                           solid_system_rhs,
                                                           cell_matrix);
            }
          else
            {
              solid_constraints.distribute_local_to_global(cell_rhs,
                                                           local_dof_indices,
                                                           solid_system_rhs);
            }
        }
    solid_system_rhs.compress(VectorOperation::add);
  }


  template <int dim, int spacedim>
  void NavierStokesImmersedProblem<dim, spacedim>::assemble_coupling_sparsity(
    BlockDynamicSparsityPattern &dsp)
  {
    TimerOutput::Scope t(computing_timer, "Assemble Coupling sparsity");

    std::vector<types::global_dof_index> fluid_dof_indices(
      fluid_fe->n_dofs_per_cell());

    FullMatrix<double> local_matrix;

    auto particle = solid_particle_handler.begin();
    while (particle != solid_particle_handler.end())
      {
        const auto &cell = particle->get_surrounding_cell();
        const auto &dh_cell =
          typename DoFHandler<spacedim>::cell_iterator(*cell, &fluid_dh);
        dh_cell->get_dof_indices(fluid_dof_indices);


        const auto pic = solid_particle_handler.particles_in_cell(cell);

        Assert(pic.begin() == particle, ExcInternalError());
        for (const auto &p : pic)
          {
            const auto global_particle_id = p.get_id();

            for (unsigned int i = 0; i < fluid_fe->n_dofs_per_cell(); ++i)
              {
                const auto comp_i =
                  fluid_fe->system_to_component_index(i).first;

                if (comp_i < spacedim)
                  dsp.add(fluid_dof_indices[i],
                          global_particle_id * spacedim + comp_i);
              }
          }
        particle = pic.end();
      }
  }



  template <int dim, int spacedim>
  void NavierStokesImmersedProblem<dim, spacedim>::assemble_coupling()
  {
    TimerOutput::Scope t(computing_timer, "Assemble Coupling terms");

    std::vector<types::global_dof_index> fluid_dof_indices(
      fluid_fe->n_dofs_per_cell());
    std::vector<types::global_dof_index> solid_dof_indices;

    FullMatrix<double> local_matrix;

    auto particle = solid_particle_handler.begin();
    while (particle != solid_particle_handler.end())
      {
        const auto &cell = particle->get_surrounding_cell();
        const auto &dh_cell =
          typename DoFHandler<spacedim>::cell_iterator(*cell, &fluid_dh);
        dh_cell->get_dof_indices(fluid_dof_indices);


        const auto pic = solid_particle_handler.particles_in_cell(cell);

        const auto n_particles_in_cell =
          solid_particle_handler.n_particles_in_cell(cell);

        local_matrix.reinit(fluid_fe->n_dofs_per_cell(),
                            n_particles_in_cell * spacedim);

        solid_dof_indices.resize(n_particles_in_cell * spacedim);

        unsigned int local_solid_dof_index = 0;
        unsigned int local_pic_index       = 0;
        Assert(pic.begin() == particle, ExcInternalError());
        for (const auto &p : pic)
          {
            const auto &ref_q              = p.get_reference_location();
            const auto  global_particle_id = p.get_id();

            for (unsigned int d = 0; d < spacedim; ++d)
              solid_dof_indices[local_solid_dof_index++] =
                global_particle_id * spacedim + d;

            for (unsigned int i = 0; i < fluid_fe->n_dofs_per_cell(); ++i)
              {
                const auto comp_i =
                  fluid_fe->system_to_component_index(i).first;

                if (comp_i < spacedim)
                  local_matrix(i, local_pic_index * spacedim + comp_i) =
                    fluid_fe->shape_value(i, ref_q);
              }
            local_pic_index++;
          }

        fluid_constraints.distribute_local_to_global(
          local_matrix,
          fluid_dof_indices,
          solid_constraints,
          solid_dof_indices,
          coupling_interpolation_matrix);
        particle = pic.end();
      }

    coupling_interpolation_matrix.compress(VectorOperation::add);
  }



  // @sect4{Solving the coupled system}
  template <int dim, int spacedim>
  void NavierStokesImmersedProblem<dim, spacedim>::solve(const double time_step)
  {
    TimerOutput::Scope t(computing_timer, "Solve");

    using Vec   = LA::MPI::Vector;
    using LinOp = LinearOperator<Vec>;

    const auto A  = LinOp(fluid_matrix.block(0, 0));
    const auto Bt = LinOp(fluid_matrix.block(0, 1));
    const auto B  = LinOp(fluid_matrix.block(1, 0));
    const auto Z6 = 0.0 * LinOp(fluid_matrix.block(1, 1));
    const auto Mp = LinOp(fluid_preconditioner.block(1, 1));

    const auto K  = LinOp(solid_matrix.block(0, 0));
    const auto Dt = LinOp(solid_matrix.block(0, 1));
    const auto D  = LinOp(solid_matrix.block(1, 0));
    const auto M  = LinOp(solid_matrix.block(1, 1));

    const auto Pt  = LinOp(coupling_interpolation_matrix.block(0, 0));
    const auto Z1t = 0.0 * LinOp(coupling_interpolation_matrix.block(0, 1));
    const auto Z2t = 0.0 * LinOp(coupling_interpolation_matrix.block(1, 0));
    const auto Z3t = 0.0 * LinOp(coupling_interpolation_matrix.block(1, 1));

    const auto P  = transpose_operator(Pt);
    const auto Z1 = transpose_operator(Z1t);
    const auto Z2 = transpose_operator(Z2t);
    const auto Z3 = transpose_operator(Z3t);

    const auto Z4 = 0.0 * M;
    using BVec    = typename LA::MPI::BlockVector;

    auto MC  = -time_step * D * P;
    auto CtM = -time_step * Pt * Dt;

    // Inversion of the mass matrices
    SolverControl inner_solver_control(par.inner_max_iterations,
                                       par.inner_tolerance,
                                       false,
                                       false);
    SolverCG<TrilinosWrappers::MPI::Vector> cg_solver(inner_solver_control);

    TrilinosWrappers::PreconditionILU Mp_inv_ilu;
    Mp_inv_ilu.initialize(fluid_preconditioner.block(1, 1));
    auto invMp = inverse_operator(Mp, cg_solver, Mp_inv_ilu);

    TrilinosWrappers::PreconditionILU M_inv_ilu;
    M_inv_ilu.initialize(solid_matrix.block(1, 1));
    auto         invM = inverse_operator(M, cg_solver, M_inv_ilu);
    const double h    = GridTools::maximal_cell_diameter(solid_tria);
    const auto   invW = invM;

    const auto gamma1 = par.gamma_AL_background;
    const auto gamma2 = par.gamma_AL_immersed * time_step;

    auto A11_aug = A + gamma1 * Pt * M * P + gamma1 * Bt * invMp * B;
    auto A22_aug = (1 / time_step) * K + gamma2 * Dt * invW * D;
    auto A12_aug = gamma1 * CtM * invW * D;
    auto A21_aug = gamma2 * Dt * invW * MC;

    SolverControl inner_solver_control_lagrangian(
      par.inner_lagrangian_max_iterations,
      par.inner_lagrangian_tolerance,
      false,
      false);
    SolverCG<Vec> cg_solver_lagrangian(inner_solver_control_lagrangian);

    TrilinosWrappers::PreconditionAMG amg_A;
    amg_A.initialize(fluid_matrix.block(0, 0));
    auto A11_aug_inv =
      inverse_operator(A11_aug,
                       cg_solver_lagrangian,
                       amg_A); // inverse of fluid velocity block

    TrilinosWrappers::PreconditionAMG amg_K;
    amg_K.initialize(solid_matrix.block(0, 0));
    auto A22_aug_inv =
      inverse_operator(A22_aug,
                       cg_solver_lagrangian,
                       amg_K); // inverse of solid displacement block

    StokesDLMALPreconditioner<Vec, BVec> prec_AL(A11_aug_inv,
                                                 A22_aug_inv,
                                                 A12_aug,
                                                 Bt,
                                                 CtM,
                                                 invW,
                                                 invMp,
                                                 Dt,
                                                 gamma1,
                                                 gamma2,
                                                 gamma1);


    // std::array<std::array<LinOp, 4>, 4> system_array = {
    //   {{{A, Z1t, CtM, Bt}},  // vel
    //    {{Z1, K, Dt, Z2t}},   // disp
    //    {{MC, D, Z4, Z3}},    // lagr
    //    {{B, Z2, Z3t, Z6}}}}; // pres

    std::array<std::array<LinOp, 4>, 4> system_array = {
      {{{A11_aug, A12_aug, CtM, Bt}}, // vel
       {{A21_aug, A22_aug, Dt, Z2t}}, // disp
       {{MC, D, Z4, Z3}},             // lagr
       {{B, Z2, Z3t, Z6}}}};          // pres

    const auto system = block_operator<4, 4, BVec>(system_array);
    BVec       block_system_rhs(4), block_system_solution(4);
    block_system_rhs.block(0) = fluid_system_rhs.block(0);
    block_system_rhs.block(1) = solid_system_rhs.block(0);
    block_system_rhs.block(2) = solid_system_rhs.block(1);
    block_system_rhs.block(3) = fluid_system_rhs.block(1);

    block_system_rhs.block(0) +=
      gamma1 * CtM * invW * block_system_rhs.block(2);
    block_system_rhs.block(1) += gamma2 * Dt * invW * block_system_rhs.block(2);


    block_system_solution.reinit(block_system_rhs);
    block_system_solution.block(0) =
      fluid_locally_relevant_solution.block(0); // u
    block_system_solution.block(1) =
      solid_locally_relevant_solution.block(0); // w
    block_system_solution.block(2) =
      solid_locally_relevant_solution.block(1); // lambda
    block_system_solution.block(3) =
      fluid_locally_relevant_solution.block(1); // p

    SolverControl                      solver_control(par.outer_max_iterations,
                                 par.outer_tolerance,
                                 true,
                                 false);
    SolverFGMRES<LA::MPI::BlockVector> solver(solver_control);

    fluid_constraints.set_zero(fluid_solution);

    solver.solve(system, block_system_solution, block_system_rhs, prec_AL);
    fluid_solution.block(0) = block_system_solution.block(0);
    solid_solution.block(0) = block_system_solution.block(1);
    solid_solution.block(1) = block_system_solution.block(2);
    fluid_solution.block(1) = block_system_solution.block(3);

    pcout << "   Solved in " << solver_control.last_step() << " iterations."
          << std::endl;

    fluid_constraints.distribute(fluid_solution);

    if (fluid_constant_pressure.size() == 0)
      {
        fluid_constant_pressure.reinit(fluid_solution);
        fluid_constant_pressure.block(1) =
          invMp * fluid_dual_of_constant_pressure.block(1);
      }

    const auto avg_pressure = fluid_dual_of_constant_pressure * fluid_solution;

    fluid_solution.block(1).add(-avg_pressure,
                                fluid_constant_pressure.block(1));

    fluid_locally_relevant_solution = fluid_solution;
    solid_locally_relevant_solution = solid_solution;
  }



  // @sect4{Mesh refinement}
  template <int dim, int spacedim>
  void NavierStokesImmersedProblem<dim, spacedim>::refine_and_transfer()
  {
    TimerOutput::Scope t(computing_timer, "Refine");

    Vector<float> error_per_cell(fluid_tria.n_active_cells());
    KellyErrorEstimator<spacedim>::estimate(fluid_dh,
                                            QGauss<spacedim - 1>(
                                              par.velocity_degree + 1),
                                            {},
                                            fluid_locally_relevant_solution,
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

    SolutionTransfer<spacedim, LA::MPI::BlockVector> transfer(fluid_dh);

    fluid_tria.prepare_coarsening_and_refinement();
    transfer.prepare_for_coarsening_and_refinement(
      fluid_locally_relevant_solution);
    tracer_particle_handler.prepare_for_coarsening_and_refinement();
    solid_particle_handler.prepare_for_coarsening_and_refinement();

    fluid_tria.execute_coarsening_and_refinement();

    setup_dofs();

    transfer.interpolate(fluid_solution);
    tracer_particle_handler.unpack_after_coarsening_and_refinement();
    solid_particle_handler.unpack_after_coarsening_and_refinement();

    fluid_constraints.distribute(fluid_solution);
    fluid_locally_relevant_solution = fluid_solution;
  }



  // @sect4{Output and visualization}
  template <int dim, int spacedim>
  void NavierStokesImmersedProblem<dim, spacedim>::output_results(
    const unsigned int cycle,
    double             time) const
  {
    // Fluid solution
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
      data_out.add_data_vector(fluid_locally_relevant_solution,
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
      times_and_names.emplace_back(std::make_pair(time, filename));
      std::ofstream ofile(par.output_directory + "/" + "solution.pvd");
      DataOutBase::write_pvd_record(ofile, times_and_names);
    }

    // Solid solution
    {
      TimerOutput::Scope t(computing_timer, "Output solid");

      std::vector<std::string> solution_names(2 * spacedim);
      for (unsigned int i = 0; i < spacedim; ++i)
        solution_names[i] = "displacement";
      for (unsigned int i = spacedim; i < 2 * spacedim; ++i)
        solution_names[i] = "lagrange_multiplier";

      std::vector<DataComponentInterpretation::DataComponentInterpretation>
        data_component_interpretation(
          2 * spacedim,
          DataComponentInterpretation::component_is_part_of_vector);

      DataOut<spacedim> data_out;
      data_out.attach_dof_handler(solid_dh);
      data_out.add_data_vector(solid_locally_relevant_solution,
                               solution_names,
                               DataOut<spacedim>::type_dof_data,
                               data_component_interpretation);


      Vector<float> subdomain(solid_tria.n_active_cells());
      for (unsigned int i = 0; i < subdomain.size(); ++i)
        subdomain(i) = solid_tria.locally_owned_subdomain();
      data_out.add_data_vector(subdomain, "subdomain");

      data_out.build_patches();

      const std::string filename =
        "solution-solid-" + Utilities::int_to_string(cycle) + ".vtu";
      data_out.write_vtu_in_parallel(par.output_directory + "/" + filename,
                                     mpi_communicator);

      static std::vector<std::pair<double, std::string>> times_and_names;
      times_and_names.emplace_back(std::make_pair(time, filename));
      std::ofstream ofile(par.output_directory + "/" + "solution-solid.pvd");
      DataOutBase::write_pvd_record(ofile, times_and_names);
    }
    output_particles(solid_particle_handler,
                     "solution-solid-particles",
                     cycle,
                     time);
  }


  template <int dim, int spacedim>
  void NavierStokesImmersedProblem<dim, spacedim>::output_particles(
    const Particles::ParticleHandler<spacedim> &particles,
    const std::string                          &fprefix,
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
      times_and_names[fprefix].emplace_back(std::make_pair(time, filename));
    else
      times_and_names[fprefix] = {std::make_pair(time, filename)};
    std::ofstream ofile(par.output_directory + "/" + fprefix + ".pvd");
    DataOutBase::write_pvd_record(ofile, times_and_names[fprefix]);
  }

  template <int dim, int spacedim>
  inline void
  NavierStokesImmersedProblem<dim, spacedim>::update_particle_positions()
  {
    // After the first time step, we displace the solid body to take
    // into account the fact it has moved
    pcout << "Number of particles "
          << solid_particle_handler.n_global_particles() << std::endl;

    solid_current_position = solid_reference_configuration;
    solid_current_position += solid_locally_relevant_solution;

    // Create an intermediate vector to store particle positions
    // The issue is that solid_current_position.block(0) is distributed
    // according to the solid mesh, but particles are distributed according
    // to the fluid mesh. We need an intermediate vector with locally_relevant
    // index sets that match the particle distribution.
    const IndexSet locally_owned_particle_coords =
      solid_particle_handler.locally_owned_particle_ids().tensor_product(
        complete_index_set(spacedim));

    const IndexSet locally_relevant_particle_coords =
      locally_owned_particle_coords;

    LA::MPI::Vector particle_positions;
    particle_positions.reinit(locally_owned_particle_coords,
                              locally_relevant_particle_coords,
                              mpi_communicator);

    // Copy the current positions into the particle positions vector
    particle_positions = solid_current_position.block(0);

    solid_particle_handler.set_particle_positions(particle_positions, false);
  };



  // @sect4{Time loop}
  template <int dim, int spacedim>
  void NavierStokesImmersedProblem<dim, spacedim>::run()
  {
#ifdef USE_PETSC_LA
    pcout << "Running NavierStokesImmersedProblem<"
          << Utilities::dim_string(dim, spacedim) << "> using PETSc."
          << std::endl;
#else
    pcout << "Running NavierStokesImmersedProblem<"
          << Utilities::dim_string(dim, spacedim) << "> using Trilinos."
          << std::endl;
#endif
    par.prm.print_parameters(par.output_directory + "/" + "used_parameters_" +
                               std::to_string(dim) + std::to_string(spacedim) +
                               ".prm",
                             ParameterHandler::Short);

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
            interpolate_initial_conditions();
            setup_solid_particles();
            output_results(output_cycle, time);

            assemble_navier_stokes_system(time_step);
            assemble_elasticity_system(time_step);

            assemble_navier_stokes_rhs(time_step);
            assemble_elasticity_rhs(time_step);
          }
        else
          {
            assemble_navier_stokes_rhs(time_step);
            assemble_elasticity_rhs(time_step);
          }

        setup_coupling();
        assemble_coupling();

        solve(time_step);

        update_particle_positions();

        if (cycle % par.output_frequency == 0)
          {
            output_results(output_cycle + 1, time + time_step);
            ++output_cycle;
          }
        // if (cycle % par.refinement_frequency == 0 &&
        //     cycle != par.number_of_time_steps - 1)
        //   refine_and_transfer();

        fluid_locally_relevant_solution_old = fluid_locally_relevant_solution;
        solid_locally_relevant_solution_old = solid_locally_relevant_solution;
      }
  }
} // namespace Step80



// @sect3{The main() function}
int main(int argc, char *argv[])
{
  using namespace Step80;
  using namespace dealii;
  deallog.depth_console(10);
  try
    {
      Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

      std::string prm_file;
      if (argc > 1)
        prm_file = argv[1];
      else
        prm_file = "parameters.prm";

      // Extract the dimension from the parameter file.
      auto [dim, spacedim] = get_dimension_and_spacedimension(prm_file);

      if (dim == 2 && spacedim == 2)
        {
          NavierStokesImmersedProblemParameters<2> par;
          ParameterAcceptor::initialize(prm_file);

          NavierStokesImmersedProblem<2> problem(par);
          problem.run();
        }
      else if (dim == 3 && spacedim == 3)
        {
          NavierStokesImmersedProblemParameters<3> par;
          ParameterAcceptor::initialize(prm_file);

          NavierStokesImmersedProblem<3> problem(par);
          problem.run();
        }
      else
        {
          AssertThrow(false,
                      ExcNotImplemented(
                        "The combination of dimension " + std::to_string(dim) +
                        " and spacedimension " + std::to_string(spacedim) +
                        " is not implemented."));
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
