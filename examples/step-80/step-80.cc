/* ---------------------------------------------------------------------
 *
 * Copyright (C) 2024 by the deal.II authors
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


#include <deal.II/base/function.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/timer.h>

#include <deal.II/lac/block_linear_operator.h>
#include <deal.II/lac/generic_linear_algebra.h>
#include <deal.II/lac/linear_operator.h>
#include <deal.II/lac/linear_operator_tools.h>

#include <boost/algorithm/string.hpp>

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

#include <deal.II/particles/data_out.h>
#include <deal.II/particles/generators.h>
#include <deal.II/particles/particle_handler.h>
#include <deal.II/particles/utilities.h>

#include <deal.II/opencascade/manifold_lib.h>
#include <deal.II/opencascade/utilities.h>
#ifdef DEAL_II_WITH_OPENCASCADE
#  include <TopoDS.hxx>
#endif

#include <cmath>
#include <fstream>
#include <iostream>
#include <memory>

namespace Step80
{
  using namespace dealii;


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
    // If reading of the input file fails, run by default in 1D-2D
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


  template <int dim, int spacedim = dim>
  class StokesImmersedProblemParameters : public ParameterAcceptor
  {
  public:
    StokesImmersedProblemParameters();

    void set_time(const double &time) const
    {
      rhs.set_time(time);
      angular_velocity.set_time(time);
    }

    std::string output_directory = ".";

    unsigned int velocity_degree = 2;

    unsigned int number_of_time_steps = 501;
    double       final_time           = 1.0;

    unsigned int output_frequency = 1;

    unsigned int initial_fluid_refinement      = 5;
    unsigned int initial_solid_refinement      = 5;
    unsigned int particle_insertion_refinement = 3;

    unsigned int fluid_rtree_extraction_level = 1;

    double viscosity    = 1.0;
    double penalty_term = 100;

    std::list<types::boundary_id> homogeneous_dirichlet_ids{0};

    std::string name_of_fluid_grid       = "hyper_cube";
    std::string arguments_for_fluid_grid = "-1: 1: false";
    std::string name_of_solid_grid       = "hyper_rectangle";
    std::string arguments_for_solid_grid = spacedim == 2 ?
                                             "-.5, -.1: .5, .1: false" :
                                             "-.5, -.1, -.1: .5, .1, .1: false";
    std::string name_of_particle_grid    = "hyper_ball";
    std::string arguments_for_particle_grid =
      spacedim == 2 ? "0.3, 0.3: 0.1: false" : "0.3, 0.3, 0.3 : 0.1: false";

    int          max_level_refinement = 8;
    int          min_level_refinement = 5;
    std::string  refinement_strategy  = "fixed_fraction";
    double       coarsening_fraction  = 0.3;
    double       refinement_fraction  = 0.3;
    unsigned int max_cells            = 20000;
    int          refinement_frequency = 5;

    mutable ParameterAcceptorProxy<Functions::ParsedFunction<spacedim>> rhs;
    mutable ParameterAcceptorProxy<Functions::ParsedFunction<spacedim>>
      angular_velocity;
  };



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

    rhs.declare_parameters_call_back.connect([&]() {
      Functions::ParsedFunction<spacedim>::declare_parameters(this->prm,
                                                              spacedim + 1);
    });
    angular_velocity.declare_parameters_call_back.connect([&]() {
      this->prm.set("Function expression",
                    "t < .500001 ? 6.283185 : -6.283185");
    });

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



  template <int dim, int spacedim = dim>
  class StokesImmersedProblem
  {
  public:
    StokesImmersedProblem(
      const StokesImmersedProblemParameters<dim, spacedim> &par);

    void run();

  private:
    void make_grid();

    double compute_time_step() const;

    void setup_tracer_particles();
    void setup_solid_particles();

    void initial_setup();
    void setup_dofs();

    void assemble_stokes_system();
    void assemble_nitsche_restriction();

    void solve();

    void refine_and_transfer();

    void output_results(const unsigned int cycle, const double time) const;
    void output_particles(const Particles::ParticleHandler<spacedim> &particles,
                          std::string                                 fprefix,
                          const unsigned int                          iter,
                          const double time) const;

    const StokesImmersedProblemParameters<dim, spacedim> &par;

    MPI_Comm mpi_communicator;

    ConditionalOStream pcout;

    mutable TimerOutput computing_timer;

    parallel::distributed::Triangulation<spacedim>      fluid_tria;
    parallel::distributed::Triangulation<dim, spacedim> solid_tria;

    std::unique_ptr<FiniteElement<spacedim>>      fluid_fe;
    std::unique_ptr<FiniteElement<dim, spacedim>> solid_fe;

    DoFHandler<spacedim>      fluid_dh;
    DoFHandler<dim, spacedim> solid_dh;

    std::unique_ptr<MappingFEField<dim, spacedim>> solid_mapping;

    std::vector<IndexSet> fluid_owned_dofs;
    std::vector<IndexSet> solid_owned_dofs;

    std::vector<IndexSet> fluid_relevant_dofs;
    std::vector<IndexSet> solid_relevant_dofs;

    AffineConstraints<double> constraints;

    LA::MPI::BlockSparseMatrix fluid_matrix;
    LA::MPI::BlockSparseMatrix fluid_preconditioner;

    LA::MPI::BlockSparseMatrix solid_matrix;
    LA::MPI::BlockSparseMatrix coupling_matrix;

    LA::MPI::BlockVector solution;
    LA::MPI::BlockVector locally_relevant_solution;
    LA::MPI::BlockVector system_rhs;

    Particles::ParticleHandler<spacedim> tracer_particle_handler;
    Particles::ParticleHandler<spacedim> solid_particle_handler;

    IndexSet locally_owned_tracer_particle_coordinates;
    IndexSet locally_relevant_tracer_particle_coordinates;

    LA::MPI::Vector tracer_particle_velocities;
    LA::MPI::Vector relevant_tracer_particle_displacements;

    std::vector<std::vector<BoundingBox<spacedim>>> global_fluid_bounding_boxes;
  };



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
        else if (spacedim == 3)
          {
            const auto t = reinterpret_cast<Triangulation<dim, 3> *>(&tria);
            t->set_manifold(manifold_id,
                            OpenCASCADE::NormalToMeshProjectionManifold<dim, 3>(
                              shape));
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

    tracer_particle_handler.initialize(fluid_tria,
                                       StaticMappingQ1<spacedim>::mapping);

    Particles::Generators::dof_support_points(particles_dof_handler,
                                              global_fluid_bounding_boxes,
                                              tracer_particle_handler);

    pcout << "Tracer particles: "
          << tracer_particle_handler.n_global_particles() << std::endl;

    locally_owned_tracer_particle_coordinates =
      tracer_particle_handler.locally_owned_particle_ids().tensor_product(
        complete_index_set(spacedim));

    locally_relevant_tracer_particle_coordinates =
      locally_owned_tracer_particle_coordinates;
  }


  template <int dim, int spacedim>
  void StokesImmersedProblem<dim, spacedim>::setup_solid_particles()
  {
    QGauss<dim> quadrature(fluid_fe->degree + 1);

    const unsigned int n_properties = 1;
    solid_particle_handler.initialize(fluid_tria,
                                      StaticMappingQ1<spacedim>::mapping,
                                      n_properties);

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

    Assert(!global_fluid_bounding_boxes.empty(),
           ExcInternalError(
             "I was expecting the "
             "global_fluid_bounding_boxes to be filled at this stage. "
             "Make sure you fill this vector before trying to use it "
             "here. Bailing out."));

    solid_particle_handler.insert_global_particles(quadrature_points_vec,
                                                   global_fluid_bounding_boxes,
                                                   properties);

    pcout << "Solid particles: " << solid_particle_handler.n_global_particles()
          << std::endl;
  }



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
  }


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

    const IndexSet locally_relevant_dofs =
      DoFTools::extract_locally_relevant_dofs(fluid_dh);
    fluid_relevant_dofs.resize(2);
    fluid_relevant_dofs[0] = locally_relevant_dofs.get_view(0, n_u);
    fluid_relevant_dofs[1] = locally_relevant_dofs.get_view(n_u, n_u + n_p);

    {
      constraints.reinit(locally_relevant_dofs);

      const FEValuesExtractors::Vector velocities(0);
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

      BlockDynamicSparsityPattern dsp(dofs_per_block, dofs_per_block);

      DoFTools::make_sparsity_pattern(
        fluid_dh, coupling, dsp, constraints, false);

      SparsityTools::distribute_sparsity_pattern(
        dsp,
        locally_owned_dofs_per_processor,
        mpi_communicator,
        locally_relevant_dofs);

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

      BlockDynamicSparsityPattern dsp(dofs_per_block, dofs_per_block);

      DoFTools::make_sparsity_pattern(
        fluid_dh, coupling, dsp, constraints, false);
      SparsityTools::distribute_sparsity_pattern(
        dsp,
        locally_owned_dofs_per_processor,
        mpi_communicator,
        locally_relevant_dofs);
      fluid_preconditioner.reinit(fluid_owned_dofs, dsp, mpi_communicator);
    }

    locally_relevant_solution.reinit(fluid_owned_dofs,
                                     fluid_relevant_dofs,
                                     mpi_communicator);
    system_rhs.reinit(fluid_owned_dofs, mpi_communicator);
    solution.reinit(fluid_owned_dofs, mpi_communicator);
  }



  template <int dim, int spacedim>
  void StokesImmersedProblem<dim, spacedim>::assemble_stokes_system()
  {
    fluid_matrix = 0;
    system_rhs   = 0;

    TimerOutput::Scope t(computing_timer, "Assemble Stokes terms");

    QGauss<spacedim>   quadrature_formula(fluid_fe->degree + 1);
    FEValues<spacedim> fe_values(*fluid_fe,
                                 quadrature_formula,
                                 update_values | update_gradients |
                                   update_quadrature_points |
                                   update_JxW_values);

    const unsigned int dofs_per_cell = fluid_fe->n_dofs_per_cell();
    const unsigned int n_q_points    = quadrature_formula.size();

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
          constraints.distribute_local_to_global(
            cell_matrix, cell_rhs, local_dof_indices, fluid_matrix, system_rhs);

          constraints.distribute_local_to_global(cell_matrix2,
                                                 local_dof_indices,
                                                 fluid_preconditioner);
        }

    fluid_matrix.compress(VectorOperation::add);
    fluid_preconditioner.compress(VectorOperation::add);
    system_rhs.compress(VectorOperation::add);
  }


  template <int dim, int spacedim>
  void StokesImmersedProblem<dim, spacedim>::assemble_nitsche_restriction()
  {
    TimerOutput::Scope t(computing_timer, "Assemble Nitsche terms");

    const FEValuesExtractors::Vector velocities(0);
    const FEValuesExtractors::Scalar pressure(spacedim);

    // SolidVelocity<spacedim> solid_velocity(par.angular_velocity);

    std::vector<types::global_dof_index> fluid_dof_indices(
      fluid_fe->n_dofs_per_cell());

    FullMatrix<double>     local_matrix(fluid_fe->n_dofs_per_cell(),
                                    fluid_fe->n_dofs_per_cell());
    dealii::Vector<double> local_rhs(fluid_fe->n_dofs_per_cell());

    const auto penalty_parameter =
      1.0 / GridTools::minimal_cell_diameter(fluid_tria);

    auto particle = solid_particle_handler.begin();
    while (particle != solid_particle_handler.end())
      {
        local_matrix = 0;
        local_rhs    = 0;

        const auto &cell = particle->get_surrounding_cell();
        const auto &dh_cell =
          typename DoFHandler<spacedim>::cell_iterator(*cell, &fluid_dh);
        dh_cell->get_dof_indices(fluid_dof_indices);

        const auto pic = solid_particle_handler.particles_in_cell(cell);
        Assert(pic.begin() == particle, ExcInternalError());
        for (const auto &p : pic)
          {
            const auto &ref_q = p.get_reference_location();
            // const auto &real_q = p.get_location();
            const auto &JxW = p.get_properties()[0];

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
                                    1 * // [TODO]
                                    fluid_fe->shape_value(i, ref_q) * JxW;
                  }
              }
          }

        constraints.distribute_local_to_global(
          local_matrix, local_rhs, fluid_dof_indices, fluid_matrix, system_rhs);
        particle = pic.end();
      }

    fluid_matrix.compress(VectorOperation::add);
    system_rhs.compress(VectorOperation::add);
  }



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
      prec_A.initialize(fluid_matrix.block(0, 0), data);
    }

    LA::MPI::PreconditionAMG prec_S;
    {
      LA::MPI::PreconditionAMG::AdditionalData data;

#ifdef USE_PETSC_LA
      data.symmetric_operator = true;
#endif
      prec_S.initialize(fluid_preconditioner.block(1, 1), data);
    }

    const auto A = linear_operator<LA::MPI::Vector>(fluid_matrix.block(0, 0));
    const auto amgA = linear_operator(A, prec_A);

    const auto S =
      linear_operator<LA::MPI::Vector>(fluid_preconditioner.block(1, 1));
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

    SolverControl solver_control(fluid_matrix.m(),
                                 1e-10 * system_rhs.l2_norm());

    SolverFGMRES<LA::MPI::BlockVector> solver(solver_control);

    constraints.set_zero(solution);

    solver.solve(fluid_matrix, solution, system_rhs, P);


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
    tracer_particle_handler.prepare_for_coarsening_and_refinement();
    solid_particle_handler.prepare_for_coarsening_and_refinement();

    fluid_tria.execute_coarsening_and_refinement();

    setup_dofs();

    transfer.interpolate(solution);
    tracer_particle_handler.unpack_after_coarsening_and_refinement();
    solid_particle_handler.unpack_after_coarsening_and_refinement();

    constraints.distribute(solution);
    locally_relevant_solution = solution;
  }



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
        else
          {
            TimerOutput::Scope t(computing_timer,
                                 "Set solid particle position");

            // SolidPosition<spacedim> solid_position(par.angular_velocity,
            //                                        time_step);
            // solid_particle_handler.set_particle_positions(solid_position,
            //                                               false);
            // [TODO]
          }

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

        assemble_stokes_system();
        assemble_nitsche_restriction();
        solve();

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
} // namespace Step80



int main(int argc, char *argv[])
{
  using namespace Step80;
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

      // Extract the dimension from the parameter file.
      auto [dim, spacedim] = get_dimension_and_spacedimension(prm_file);

      if (dim == 2 && spacedim == 2)
        {
          StokesImmersedProblemParameters<2> par;
          ParameterAcceptor::initialize(prm_file);

          StokesImmersedProblem<2> problem(par);
          problem.run();
        }
      else if (dim == 2 && spacedim == 3)
        {
          StokesImmersedProblemParameters<2, 3> par;
          ParameterAcceptor::initialize(prm_file);

          StokesImmersedProblem<2, 3> problem(par);
          problem.run();
        }
      else if (dim == 3 && spacedim == 3)
        {
          StokesImmersedProblemParameters<3> par;
          ParameterAcceptor::initialize(prm_file);

          StokesImmersedProblem<3> problem(par);
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
