/* ------------------------------------------------------------------------
 *
 * SPDX-License-Identifier: LGPL-2.1-or-later
 * Copyright (C) 2016 - 2024 by the deal.II authors
 *
 * This file is part of the deal.II library.
 *
 * Part of the source code is dual licensed under Apache-2.0 WITH
 * LLVM-exception OR LGPL-2.1-or-later. Detailed license information
 * governing the source code and code contributions can be found in
 * LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
 *
 * ------------------------------------------------------------------------
 *
 * Author: Wolfgang Bangerth, Colorado State University, 2024
 *         Jean-Paul Pelteret, 2024.
 */

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function.h>
#include <deal.II/base/timer.h>
#include <deal.II/base/function_lib.h>

#include <deal.II/differentiation/ad.h>

#include <deal.II/lac/block_linear_operator.h>
#include <deal.II/lac/generic_linear_algebra.h>
#include <deal.II/lac/linear_operator.h>

/* #define FORCE_USE_OF_TRILINOS */

#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/solver_gmres.h>
#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>

#include <deal.II/lac/petsc_sparse_matrix.h>
#include <deal.II/lac/petsc_vector.h>
#include <deal.II/lac/petsc_solver.h>
#include <deal.II/lac/petsc_precondition.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_renumbering.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/error_estimator.h>

#include <deal.II/base/utilities.h>
#include <deal.II/base/conditional_ostream.h>
#include <deal.II/base/index_set.h>
#include <deal.II/lac/sparsity_tools.h>
#include <deal.II/distributed/tria.h>
#include <deal.II/distributed/grid_refinement.h>

#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/device/file.hpp>

#include <deal.II/sundials/ida.h>

#include <cmath>
#include <fstream>
#include <iostream>

namespace Step55
{
  using namespace dealii;

  namespace LA
  {
#if defined(DEAL_II_WITH_PETSC) && !defined(DEAL_II_PETSC_WITH_COMPLEX) && \
  !(defined(DEAL_II_WITH_TRILINOS) && defined(FORCE_USE_OF_TRILINOS))
    using namespace LinearAlgebraPETSc;
#  define USE_PETSC_LA
#elif defined(DEAL_II_WITH_TRILINOS)
    using namespace LinearAlgebraTrilinos;
#else
#  error DEAL_II_WITH_PETSC or DEAL_II_WITH_TRILINOS required
#endif
  } // namespace LA

  using VectorType = LA::MPI::BlockVector;
  using MatrixType = LA::MPI::BlockSparseMatrix;


  namespace ModelParameters
  {
    constexpr double stream_power_exponent_m    = 0.4;
    constexpr double stream_power_exponent_n    = 1;
    constexpr double stream_power_coefficient_k = 2e-5;
    constexpr double diffusion_coefficient_Kd   = 0.01;
    constexpr double rainfall_rate_p            = 0.6;
    constexpr double regularization_epsilon     = 0.0001;
    constexpr double stabilization_c            = 0.1;
  } // namespace ModelParameters


  class ColoradoTopography : public Function<2>
  {
  public:
    ColoradoTopography(const MPI_Comm mpi_communicator);

    virtual double value(const Point<2> &p,
                         const unsigned int /*component*/ = 0) const override
    {
      return data->value(p);
    }

  private:
    std::unique_ptr<Functions::InterpolatedUniformGridData<2>> data;
  };


  // TODO: share data among processors
  // Exact: -109 to -102, 37 to 41
  ColoradoTopography::ColoradoTopography(const MPI_Comm mpi_communicator)
  {
    unsigned int n_rows;
    unsigned int n_columns;
    Point<2>     lower_left_corner;
    double       pixel_size;

    Table<2, double> elevation_data;

    if (Utilities::MPI::this_mpi_process(mpi_communicator) == 0)
      {
        const std::string filename = "colorado-topography-1800m.txt.gz";

        boost::iostreams::filtering_istream in;
        in.push(boost::iostreams::basic_gzip_decompressor<>());
        in.push(boost::iostreams::file_source(filename));

        std::string word;

        in >> word;
        AssertThrow(word == "ncols",
                    ExcMessage(
                      "The first line of the input file needs to start "
                      "with the word 'ncols', but starts with '" +
                      word + "'."));
        in >> n_columns;

        in >> word;
        AssertThrow(word == "nrows",
                    ExcMessage(
                      "The second line of the input file needs to start "
                      "with the word 'nrows', but starts with '" +
                      word + "'."));
        in >> n_rows;

        in >> word;
        AssertThrow(word == "xllcorner",
                    ExcMessage(
                      "The third line of the input file needs to start "
                      "with the word 'xllcorner', but starts with '" +
                      word + "'."));
        in >> lower_left_corner[0];

        in >> word;
        AssertThrow(word == "yllcorner",
                    ExcMessage(
                      "The fourth line of the input file needs to start "
                      "with the word 'yllcorner', but starts with '" +
                      word + "'."));
        in >> lower_left_corner[1];


        in >> word;
        AssertThrow(word == "cellsize",
                    ExcMessage(
                      "The fourth line of the input file needs to start "
                      "with the word 'cellsize', but starts with '" +
                      word + "'."));
        in >> pixel_size;

        elevation_data.reinit(n_columns, n_rows);
        for (unsigned int row = 0; row < n_rows; ++row)
          for (unsigned int column = 0; column < n_columns; ++column)
            {
              try
                {
                  double elevation;
                  in >> elevation;

                  elevation_data(column, n_rows - row - 1) = elevation;
                }
              catch (...)
                {
                  AssertThrow(false,
                              ExcMessage(
                                "Could not read all expected data points "
                                "from the file <" +
                                filename + ">!"));
                }
            }
      }

    n_rows    = Utilities::MPI::broadcast(mpi_communicator, n_rows, 0);
    n_columns = Utilities::MPI::broadcast(mpi_communicator, n_columns, 0);
    lower_left_corner =
      Utilities::MPI::broadcast(mpi_communicator, lower_left_corner, 0);
    pixel_size = Utilities::MPI::broadcast(mpi_communicator, pixel_size, 0);

    elevation_data.replicate_across_communicator(mpi_communicator, 0);

    data = std::make_unique<Functions::InterpolatedUniformGridData<2>>(
      std::array<std::pair<double, double>, 2>{
        {std::make_pair(lower_left_corner[0],
                        lower_left_corner[0] + n_columns * pixel_size),
         std::make_pair(lower_left_corner[1],
                        lower_left_corner[1] + n_rows * pixel_size)}},
      std::array<unsigned int, 2>{{n_columns - 1, n_rows - 1}},
      std::move(elevation_data));
  }



  template <int dim>
  class RainFallRate : public Function<dim>
  {
  public:
    virtual double value(const Point<dim>  &p,
                         const unsigned int component) const override;
  };


  template <int dim>
  double RainFallRate<dim>::value(const Point<dim> & /*p*/,
                                  const unsigned int /*component*/) const
  {
    return ModelParameters::rainfall_rate_p;
  }



  template <int dim>
  class StreamPowerErosionProblem
  {
  public:
    StreamPowerErosionProblem();

    void run();

  private:
    void make_grid();
    void setup_system();

    void interpolate_initial_elevation();
    void compute_initial_constraints();
    void assemble_initial_waterflow_system();
    void solve_initial_waterflow_system();

    template <typename NumberType>
    void compute_local_residual(
      const FEValues<dim>           &fe_values,
      const std::vector<NumberType> &local_solution_elevation_at_q_points,
      const std::vector<NumberType> &local_solution_water_at_q_points,
      const std::vector<Tensor<1, dim, NumberType>>
                                    &local_gradient_elevation_at_q_points,
      const std::vector<NumberType> &div_Ih_d_wh_at_q_points,
      const std::vector<double>     &local_solution_dot_elevation_at_q_points,
      const std::vector<double>     &rain_fall_rate_rhs_values,
      const double                   cell_diameter,
      const double                   time,
      std::vector<NumberType>       &cell_residual) const;

    void assemble_residual(const double      time,
                           const VectorType &locally_relevant_solution,
                           const VectorType &locally_relevant_solution_dot,
                           VectorType       &residual);
    void assemble_jacobian(const double      time,
                           const VectorType &locally_relevant_solution,
                           const VectorType &locally_relevant_solution_dot,
                           const double      alpha,
                           MatrixType       &jacobian);
    void solve_with_jacobian(const VectorType &rhs,
                             VectorType       &dst,
                             const double      tolerance);

    void refine_grid();
    void output_results(const double       time,
                        const VectorType  &locally_relevant_solution,
                        const VectorType  &locally_relevant_solution_dot,
                        const unsigned int step_number);

    MPI_Comm mpi_communicator;

    const FESystem<dim>                       fe;
    parallel::distributed::Triangulation<dim> triangulation;
    DoFHandler<dim>                           dof_handler;

    std::vector<IndexSet> owned_partitioning;
    std::vector<IndexSet> relevant_partitioning;

    AffineConstraints<double> constraints;

    MatrixType system_matrix;
    VectorType locally_relevant_solution;
    VectorType locally_relevant_solution_dot;
    VectorType system_rhs;

    ConditionalOStream pcout;
    TimerOutput        computing_timer;
  };



  template <int dim>
  StreamPowerErosionProblem<dim>::StreamPowerErosionProblem()
    : mpi_communicator(MPI_COMM_WORLD)
    , fe(FE_Q<dim>(1), FE_Q<dim>(1))
    , triangulation(mpi_communicator,
                    typename Triangulation<dim>::MeshSmoothing(
                      Triangulation<dim>::smoothing_on_refinement |
                      Triangulation<dim>::smoothing_on_coarsening))
    , dof_handler(triangulation)
    , pcout(std::cout,
            (Utilities::MPI::this_mpi_process(mpi_communicator) == 0))
    , computing_timer(mpi_communicator,
                      pcout,
                      TimerOutput::never,
                      TimerOutput::wall_times)
  {}


  template <int dim>
  void StreamPowerErosionProblem<dim>::make_grid()
  {
    pcout << "Make grid... " << std::flush;

    GridGenerator::subdivided_hyper_rectangle(triangulation,
                                              {7, 4},
                                              Point<2>(-109., 37.),
                                              Point<2>(-102., 41.));
    triangulation.refine_global(4);

    pcout << "done. " << std::endl;
  }


  template <int dim>
  void StreamPowerErosionProblem<dim>::setup_system()
  {
    TimerOutput::Scope t(computing_timer, "setup");
    pcout << "Setup system... " << std::flush;

    dof_handler.distribute_dofs(fe);
    DoFRenumbering::component_wise(dof_handler);

    const std::vector<types::global_dof_index> dofs_per_block =
      DoFTools::count_dofs_per_fe_block(dof_handler);

    const unsigned int n_elevation      = dofs_per_block[0];
    const unsigned int n_waterflow_rate = dofs_per_block[1];

    pcout << "   Number of degrees of freedom: " << dof_handler.n_dofs() << " ("
          << n_elevation << '+' << n_waterflow_rate << ')' << std::endl;

    const IndexSet &locally_owned_dofs = dof_handler.locally_owned_dofs();
    owned_partitioning                 = {
      locally_owned_dofs.get_view(0, n_elevation),
      locally_owned_dofs.get_view(n_elevation, n_elevation + n_waterflow_rate)};

    const IndexSet locally_relevant_dofs =
      DoFTools::extract_locally_relevant_dofs(dof_handler);
    relevant_partitioning = {locally_relevant_dofs.get_view(0, n_elevation),
                             locally_relevant_dofs.get_view(
                               n_elevation, n_elevation + n_waterflow_rate)};

    {
      system_matrix.clear();

      BlockDynamicSparsityPattern dsp(relevant_partitioning);
      DoFTools::make_sparsity_pattern(dof_handler, dsp, constraints);

      SparsityTools::distribute_sparsity_pattern(
        dsp,
        dof_handler.locally_owned_dofs(),
        mpi_communicator,
        locally_relevant_dofs);

      system_matrix.reinit(owned_partitioning, dsp, mpi_communicator);
    }

    locally_relevant_solution.reinit(owned_partitioning,
                                     relevant_partitioning,
                                     mpi_communicator);
    locally_relevant_solution_dot.reinit(owned_partitioning,
                                         relevant_partitioning,
                                         mpi_communicator);
    system_rhs.reinit(owned_partitioning, mpi_communicator);

    pcout << "done. " << std::endl;
  }


  template <int dim>
  void StreamPowerErosionProblem<dim>::interpolate_initial_elevation()
  {
    TimerOutput::Scope t(computing_timer,
                         "initial conditions: interpolate elevation");
    pcout << "Interpolate elevation... " << std::flush;

    const ColoradoTopography colorado_topography(mpi_communicator);
    const VectorFunctionFromScalarFunctionObject<dim> initial_values(
      [&](const Point<dim> &p) { return colorado_topography.value(p); }, 0, 2);

    VectorType interpolated;
    interpolated.reinit(owned_partitioning, MPI_COMM_WORLD);
    VectorTools::interpolate(dof_handler, initial_values, interpolated);

    locally_relevant_solution = interpolated;

    pcout << "done. " << std::endl;
  }



  template <int dim>
  void StreamPowerErosionProblem<dim>::compute_initial_constraints()
  {
    pcout << "Compute initial constraints... " << std::flush;

    const IndexSet &locally_owned_dofs = dof_handler.locally_owned_dofs();
    const IndexSet  locally_relevant_dofs =
      DoFTools::extract_locally_relevant_dofs(dof_handler);
    constraints.reinit(locally_owned_dofs, locally_relevant_dofs);

    Assert(Utilities::MPI::n_mpi_processes(mpi_communicator) == 1,
           ExcNotImplemented());

    std::vector<bool> locally_relevant_elevation_dofs_that_are_maxima(
      locally_relevant_dofs.n_elements(), true);
    std::vector<types::global_dof_index> local_dof_indices(fe.dofs_per_cell);
    for (const auto &cell : dof_handler.active_cell_iterators())
      if (cell->is_locally_owned() || cell->is_ghost())
        {
          cell->get_dof_indices(local_dof_indices);

          for (unsigned int i = 0; i < fe.dofs_per_cell; ++i)
            if (fe.system_to_component_index(i).first == 0) // is elevation
              {
                const types::global_dof_index elevation_dof_index =
                  local_dof_indices[i];

                // Test this DoF if we haven't already determined that it is
                // not a high point
                if (locally_relevant_elevation_dofs_that_are_maxima
                      [locally_relevant_dofs.nth_index_in_set(
                        elevation_dof_index)] == true)
                  {
                    const double elevation =
                      locally_relevant_solution(elevation_dof_index);

                    for (unsigned int j = 0; j < fe.dofs_per_cell; ++j)
                      if (j != i)
                        if (fe.system_to_component_index(j).first ==
                            0) // is also elevation
                          {
                            const types::global_dof_index
                              other_elevation_dof_index = local_dof_indices[j];
                            const double other_elevation =
                              locally_relevant_solution(
                                other_elevation_dof_index);
                            if (other_elevation > elevation)
                              {
                                locally_relevant_elevation_dofs_that_are_maxima
                                  [locally_relevant_dofs.nth_index_in_set(
                                    elevation_dof_index)] = false;
                                break;
                              }
                          }
                  }
              }
        }

    IndexSet locally_relevant_water_dofs_at_inflow_boundaries(
      dof_handler.n_dofs());
    {
      const Quadrature<dim - 1> face_node_points(
        fe.get_unit_face_support_points());
      FEFaceValues<dim> fe_face_values(
        fe, face_node_points, update_gradients | update_normal_vectors);

      std::vector<Tensor<1, dim>> elevation_gradients_at_face_nodes(
        face_node_points.size());
      const FEValuesExtractors::Scalar elevation(0);

      for (const auto &cell : dof_handler.active_cell_iterators())
        if (cell->is_locally_owned() || cell->is_ghost())
          for (const auto &face : cell->face_iterators())
            if (face->at_boundary())
              {
                fe_face_values.reinit(cell, face);
                fe_face_values[elevation].get_function_gradients(
                  locally_relevant_solution, elevation_gradients_at_face_nodes);
                face->get_dof_indices(local_dof_indices);

                for (unsigned int i = 0; i < fe.dofs_per_face; ++i)
                  if (fe.face_system_to_component_index(i).first ==
                      1) // is water
                    {
                      const types::global_dof_index water_dof_index =
                        local_dof_indices[i];

                      // If we haven't already determined that this water dof is
                      // at an inflow boundary, then check so:
                      if (locally_relevant_water_dofs_at_inflow_boundaries
                            .is_element(water_dof_index) == false)
                        {
                          const bool is_inflow =
                            (elevation_gradients_at_face_nodes[i] *
                               fe_face_values.normal_vector(i) >
                             0);
                          if (is_inflow)
                            locally_relevant_water_dofs_at_inflow_boundaries
                              .add_index(water_dof_index);
                        }
                    }
              }
    }

    // Convert the maps above into an AffineConstraints object.
    VectorType distributed_solution(owned_partitioning, mpi_communicator);
    for (unsigned int i = 0; i < dof_handler.n_dofs() / 2;
         ++i) // could be done better
      if (locally_relevant_elevation_dofs_that_are_maxima
            [locally_owned_dofs.nth_index_in_set(i)])
        {
          constraints.add_constraint(system_rhs.block(0).size() + i,
                                     {},
                                     0.); // could probably do better
          // For debugging: Set the value of the solution vector for the second
          // component (water) for constrained DoFs to a value that can be
          // visualized if we don't overwrite this vector during solving linear
          // systems.
          distributed_solution(i + system_rhs.block(0).size()) = 1;
        }
    for (const auto index : locally_relevant_water_dofs_at_inflow_boundaries)
      {
        if (constraints.is_constrained(index) == false)
          {
            constraints.add_constraint(index, {}, 0.);
            // For debugging: Set the value of the solution vector
            // to a value that can be visualized if
            // we don't overwrite this vector during solving linear systems.
            distributed_solution(index) = 2;
          }
        else
          distributed_solution(index) = 3;
      }
    constraints.close();
    constraints.make_consistent_in_parallel(
      locally_owned_dofs,
      DoFTools::extract_locally_relevant_dofs(dof_handler),
      mpi_communicator);

    distributed_solution.compress(VectorOperation::insert);

    // For debugging
    locally_relevant_solution.block(1) = distributed_solution.block(1);

    pcout << "done. " << std::endl;
  }



  template <int dim>
  void StreamPowerErosionProblem<dim>::assemble_initial_waterflow_system()
  {
    TimerOutput::Scope t(computing_timer,
                         "initial conditions: assemble waterflow system");
    pcout << "Assemble initial waterflow system... " << std::flush;

    system_matrix = 0;
    system_rhs    = 0;

    const QGauss<dim> quadrature_formula(fe.degree + 1);

    FEValues<dim> fe_values(fe,
                            quadrature_formula,
                            update_values | update_gradients |
                              update_quadrature_points | update_JxW_values);
    FEValues<dim> fe_values_at_node_points(
      fe,
      Quadrature<dim>(
        fe.get_unit_support_points()), // could be made more efficient
      update_gradients);

    const unsigned int dofs_per_cell = fe.n_dofs_per_cell();
    const unsigned int n_q_points    = quadrature_formula.size();

    FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);
    Vector<double>     cell_rhs(dofs_per_cell);

    const RainFallRate<dim> rain_fall_rate_rhs;
    std::vector<double>     rain_fall_rate_rhs_values(n_q_points);

    std::vector<Tensor<1, dim>> elevation_grad_at_q_points(
      fe_values.n_quadrature_points);
    std::vector<Tensor<1, dim>> d_at_q_points(fe_values.n_quadrature_points);

    std::vector<Tensor<1, dim>> elevation_grad_at_node_points(dofs_per_cell);
    std::vector<Tensor<1, dim>> d_at_node_points(dofs_per_cell);

    std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);
    const FEValuesExtractors::Scalar     elevation(0);
    const FEValuesExtractors::Scalar     water_flow_rate(1);

    for (const auto &cell : dof_handler.active_cell_iterators())
      if (cell->is_locally_owned())
        {
          cell_matrix = 0;
          cell_rhs    = 0;

          fe_values.reinit(cell);
          fe_values_at_node_points.reinit(cell);

          rain_fall_rate_rhs.value_list(fe_values.get_quadrature_points(),
                                        rain_fall_rate_rhs_values);

          fe_values[elevation].get_function_gradients(
            locally_relevant_solution, elevation_grad_at_q_points);
          for (unsigned int q = 0; q < fe_values.n_quadrature_points; ++q)
            d_at_q_points[q] =
              -elevation_grad_at_q_points[q] /
              std::sqrt(elevation_grad_at_q_points[q] *
                          elevation_grad_at_q_points[q] +
                        ModelParameters::regularization_epsilon *
                          ModelParameters::regularization_epsilon);

          fe_values_at_node_points[elevation].get_function_gradients(
            locally_relevant_solution, elevation_grad_at_node_points);
          Assert(fe_values_at_node_points.n_quadrature_points ==
                   elevation_grad_at_node_points.size(),
                 ExcInternalError()); // remove after the PR for this is merged
          for (unsigned int j = 0; j < dofs_per_cell; ++j)
            {
              Assert(j < d_at_node_points.size(), ExcInternalError());
              d_at_node_points[j] = // TODO: get j for just the w component
                -elevation_grad_at_node_points[j] /
                std::sqrt(elevation_grad_at_node_points[j] *
                            elevation_grad_at_node_points[j] +
                          ModelParameters::regularization_epsilon *
                            ModelParameters::regularization_epsilon);
            }

          for (unsigned int q = 0; q < n_q_points; ++q)
            {
              for (unsigned int i = 0; i < dofs_per_cell; ++i)
                {
                  const double phi_i_w_at_q =
                    fe_values[water_flow_rate].value(i, q);
                  const Tensor<1, dim> grad_phi_i_w_at_q =
                    fe_values[water_flow_rate].gradient(i, q);

                  for (unsigned int j = 0; j < dofs_per_cell; ++j)
                    {
                      const Tensor<1, dim> grad_phi_j_w_at_q =
                        fe_values[water_flow_rate].gradient(j, q);

                      cell_matrix(i, j) +=
                        ((phi_i_w_at_q +
                          ModelParameters::stabilization_c * cell->diameter() *
                            (d_at_q_points[q] * grad_phi_i_w_at_q)) *
                         (d_at_node_points[j] *
                          grad_phi_j_w_at_q) * // TODO: get j for just the w
                                               // component
                         fe_values.JxW(q));
                    }

                  cell_rhs(i) +=
                    ((phi_i_w_at_q + ModelParameters::stabilization_c *
                                       cell->diameter() *
                                       (d_at_q_points[q] * grad_phi_i_w_at_q)) *
                     rain_fall_rate_rhs_values[q] * fe_values.JxW(q));
                }
            }

          cell->get_dof_indices(local_dof_indices);
          constraints.distribute_local_to_global(cell_matrix,
                                                 cell_rhs,
                                                 local_dof_indices,
                                                 system_matrix,
                                                 system_rhs);
        }

    system_matrix.compress(VectorOperation::add);
    system_rhs.compress(VectorOperation::add);

    pcout << "done. " << std::endl;
  }


  template <int dim>
  void StreamPowerErosionProblem<dim>::solve_initial_waterflow_system()
  {
    TimerOutput::Scope t(computing_timer,
                         "initial conditions: solve waterflow system");
    pcout << "Solve initial waterflow system... " << std::endl;

    pcout << "    Initial residual=" << system_rhs.block(1).l2_norm()
          << std::endl;
    SolverControl solver_control(system_matrix.block(1, 1).m(),
                                 1e-6 * system_rhs.block(1).l2_norm());

    SolverGMRES<LA::MPI::Vector> solver(solver_control);
    LA::MPI::PreconditionILU     preconditioner;
    preconditioner.initialize(system_matrix.block(1, 1));

    VectorType distributed_solution(owned_partitioning, mpi_communicator);

    constraints.set_zero(distributed_solution);

    solver.solve(system_matrix.block(1, 1),
                 distributed_solution.block(1),
                 system_rhs.block(1),
                 preconditioner);

    pcout << "    Solved in " << solver_control.last_step() << " iterations."
          << std::endl;

    constraints.distribute(distributed_solution);

    locally_relevant_solution.block(1) = distributed_solution.block(1);

    pcout << "    Solve done." << std::endl;
  }



  template <int dim>
  template <typename NumberType>
  void StreamPowerErosionProblem<dim>::compute_local_residual(
    const FEValues<dim>           &fe_values,
    const std::vector<NumberType> &local_solution_elevation_at_q_points,
    const std::vector<NumberType> &local_solution_water_at_q_points,
    const std::vector<Tensor<1, dim, NumberType>>
                                  &local_gradient_elevation_at_q_points,
    const std::vector<NumberType> &div_Ih_d_wh_at_q_points,
    const std::vector<double>     &local_solution_dot_elevation_at_q_points,
    const std::vector<double>     &rain_fall_rate_rhs_values,
    const double                   cell_diameter,
    const double                   time,
    std::vector<NumberType>       &cell_residual) const
  {
    (void)time; // TODO

    const FEValuesExtractors::Scalar elevation(0);
    const FEValuesExtractors::Scalar water_flow_rate(1);

    const unsigned int H_dof = 0;
    const unsigned int w_dof = 1;

    for (const unsigned int q : fe_values.quadrature_point_indices())
      {
        const NumberType &H     = local_solution_elevation_at_q_points[q];
        const NumberType &H_dot = local_solution_dot_elevation_at_q_points[q];
        const NumberType &w     = local_solution_water_at_q_points[q];
        const Tensor<1, dim, NumberType> &grad_H =
          local_gradient_elevation_at_q_points[q];
        const NumberType &div_Ih_d_wh = div_Ih_d_wh_at_q_points[q];
        const double      p           = rain_fall_rate_rhs_values[q];
        const double     &JxW         = fe_values.JxW(q);

        (void)H;

        const NumberType S =
          std::sqrt(local_gradient_elevation_at_q_points[q] *
                      local_gradient_elevation_at_q_points[q] +
                    ModelParameters::regularization_epsilon *
                      ModelParameters::regularization_epsilon);
        const Tensor<1, dim, NumberType> d =
          -local_gradient_elevation_at_q_points[q] / S;

        constexpr double m  = ModelParameters::stream_power_exponent_m;
        constexpr double n  = ModelParameters::stream_power_exponent_n;
        constexpr double k  = ModelParameters::stream_power_coefficient_k;
        constexpr double Kd = ModelParameters::diffusion_coefficient_Kd;
        constexpr double c  = ModelParameters::stabilization_c;

        for (const unsigned int i : fe_values.dof_indices())
          {
            const unsigned int i_group = fe.system_to_base_index(i).first.first;

            if (i_group == H_dof)
              {
                const double         Nx_i = fe_values[elevation].value(i, q);
                const Tensor<1, dim> grad_Nx_i =
                  fe_values[elevation].gradient(i, q);

                cell_residual[i] -=
                  (Nx_i * (H_dot + k * std::pow(w, m) * std::pow(S, n)) +
                   grad_Nx_i * (Kd * grad_H)) *
                  JxW;
              }
            else if (i_group == w_dof)
              {
                const double Nx_i = fe_values[water_flow_rate].value(i, q);
                const Tensor<1, dim> grad_Nx_i =
                  fe_values[water_flow_rate].gradient(i, q);
                const auto stabNx_i =
                  (Nx_i + c * cell_diameter * d * grad_Nx_i);

                // TODO: Get this right with the I_h part
                cell_residual[i] -= (stabNx_i * (p - div_Ih_d_wh)) * JxW;
              }
            else
              {
                AssertThrow(i_group <= w_dof, ExcMessage("Unknown DoF group"));
              }
          }
      }
  }



  template <int dim>
  void StreamPowerErosionProblem<dim>::assemble_residual(
    const double      time,
    const VectorType &locally_relevant_solution,
    const VectorType &locally_relevant_solution_dot,
    VectorType       &residual)
  {
    TimerOutput::Scope t(computing_timer, "linear system: assemble residual");

    residual = 0;

    const QGauss<dim> quadrature_formula(fe.degree + 1);

    FEValues<dim> fe_values(fe,
                            quadrature_formula,
                            update_values | update_gradients |
                              update_quadrature_points | update_JxW_values);
    FEValues<dim> fe_values_at_node_points(
      fe,
      Quadrature<dim>(
        fe.get_unit_support_points()), // could be made more efficient
      update_gradients);

    const unsigned int dofs_per_cell = fe.n_dofs_per_cell();
    const unsigned int n_q_points    = quadrature_formula.size();

    const RainFallRate<dim> rain_fall_rate_rhs;
    std::vector<double>     rain_fall_rate_rhs_values(n_q_points);

    std::vector<double> elevation_at_q_points(fe_values.n_quadrature_points);
    std::vector<double> elevation_dot_at_q_points(
      fe_values.n_quadrature_points);
    std::vector<double> water_at_q_points(fe_values.n_quadrature_points);
    std::vector<Tensor<1, dim>> elevation_grad_at_q_points(
      fe_values.n_quadrature_points);

    std::vector<Tensor<1, dim>> elevation_grad_at_node_points(dofs_per_cell);
    std::vector<double> div_Ih_d_wh_at_q_points(fe_values.n_quadrature_points);

    std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);
    const FEValuesExtractors::Scalar     elevation(0);
    const FEValuesExtractors::Scalar     water_flow_rate(1);

    std::vector<double> local_dof_values(dofs_per_cell);
    std::vector<double> cell_residual(dofs_per_cell);

    for (const auto &cell : dof_handler.active_cell_iterators())
      if (cell->is_locally_owned())
        {
          for (const unsigned int i : fe_values.dof_indices())
            cell_residual[i] = 0.0;

          fe_values.reinit(cell);
          fe_values_at_node_points.reinit(cell);

          cell->get_dof_indices(local_dof_indices);
          locally_relevant_solution.extract_subvector_to(local_dof_indices,
                                                         local_dof_values);

          rain_fall_rate_rhs.value_list(fe_values.get_quadrature_points(),
                                        rain_fall_rate_rhs_values);

          fe_values[elevation].get_function_values(locally_relevant_solution,
                                                   elevation_at_q_points);
          fe_values[elevation].get_function_values(
            locally_relevant_solution_dot, elevation_dot_at_q_points);
          fe_values[elevation].get_function_gradients(
            locally_relevant_solution, elevation_grad_at_q_points);
          fe_values[water_flow_rate].get_function_values(
            locally_relevant_solution, water_at_q_points);

          fe_values_at_node_points[elevation].get_function_gradients(
            locally_relevant_solution, elevation_grad_at_node_points);
          Assert(fe_values_at_node_points.n_quadrature_points ==
                   elevation_grad_at_node_points.size(),
                 ExcInternalError());
          for (const unsigned int j : fe_values.dof_indices())
            {
              // TODO: We need some sort of offsets here, because the
              // support points of W can coincide with support points
              // for H, but ultimately the solution coefficients are
              // never both non-zero for the same index j.
              // As an intermediate step, we've hack something up to
              // work around this, under the assumption that all H dofs
              // are enumerated before all w dofs, and that they're
              // enumerated in a way such that traversal of these two
              // subblocks in lock-step implies visiting the same support
              // point.

              if (j % 2 == 1) // if j is a water DoF    if
                              // (fe.shape_function_belongs_to(water_flow_rate))
                {
                  const unsigned int jj = j / 2; // TODO: do this better

                  const double         Wj = local_dof_values[j];
                  const Tensor<1, dim> d_j =
                    -elevation_grad_at_node_points[jj] /
                    std::sqrt(elevation_grad_at_node_points[jj] *
                                elevation_grad_at_node_points[jj] +
                              ModelParameters::regularization_epsilon *
                                ModelParameters::regularization_epsilon);

                  for (const unsigned int q :
                       fe_values.quadrature_point_indices())
                    {
                      const Tensor<1, dim> grad_Nx_j =
                        fe_values[water_flow_rate].gradient(jj, q);

                      div_Ih_d_wh_at_q_points[q] += Wj * d_j * grad_Nx_j;
                    }
                }
            }

          compute_local_residual(fe_values,
                                 elevation_at_q_points,
                                 water_at_q_points,
                                 elevation_grad_at_q_points,
                                 div_Ih_d_wh_at_q_points,
                                 elevation_dot_at_q_points,
                                 rain_fall_rate_rhs_values,
                                 cell->diameter(),
                                 time,
                                 cell_residual);

          constraints.distribute_local_to_global(cell_residual,
                                                 local_dof_indices,
                                                 residual);
        }

    residual.compress(VectorOperation::add);
  }



  template <int dim>
  void StreamPowerErosionProblem<dim>::assemble_jacobian(
    const double      time,
    const VectorType &locally_relevant_solution,
    const VectorType &locally_relevant_solution_dot,
    const double      alpha,
    MatrixType       &jacobian)
  {
    using ADHelper = Differentiation::AD::ResidualLinearization<
      Differentiation::AD::NumberTypes::sacado_dfad,
      double>;
    using ADNumberType = typename ADHelper::ad_type;

    TimerOutput::Scope t(computing_timer, "linear system: assemble jacobian");

    jacobian = 0;

    const QGauss<dim> quadrature_formula(fe.degree + 1);

    FEValues<dim> fe_values(fe,
                            quadrature_formula,
                            update_values | update_gradients |
                              update_quadrature_points | update_JxW_values);
    FEValues<dim> fe_values_at_node_points(
      fe,
      Quadrature<dim>(
        fe.get_unit_support_points()), // could be made more efficient
      update_gradients);

    const unsigned int dofs_per_cell = fe.n_dofs_per_cell();
    const unsigned int n_q_points    = quadrature_formula.size();

    const unsigned int n_independent_variables = dofs_per_cell;
    const unsigned int n_dependent_variables   = dofs_per_cell;
    ADHelper ad_helper(n_independent_variables, n_dependent_variables);

    FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);

    const RainFallRate<dim> rain_fall_rate_rhs;
    std::vector<double>     rain_fall_rate_rhs_values(n_q_points);

    std::vector<ADNumberType> elevation_at_q_points(
      fe_values.n_quadrature_points);
    std::vector<double> elevation_dot_at_q_points(
      fe_values.n_quadrature_points);
    std::vector<ADNumberType> water_at_q_points(fe_values.n_quadrature_points);
    std::vector<Tensor<1, dim, ADNumberType>> elevation_grad_at_q_points(
      fe_values.n_quadrature_points);

    std::vector<Tensor<1, dim, ADNumberType>> elevation_grad_at_node_points(
      dofs_per_cell);
    std::vector<ADNumberType> div_Ih_d_wh_at_q_points(
      fe_values.n_quadrature_points);

    std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);
    const FEValuesExtractors::Scalar     elevation(0);
    const FEValuesExtractors::Scalar     water_flow_rate(1);

    for (const auto &cell : dof_handler.active_cell_iterators())
      if (cell->is_locally_owned())
        {
          cell_matrix = 0;
          ad_helper.reset();

          fe_values.reinit(cell);
          fe_values_at_node_points.reinit(cell);

          cell->get_dof_indices(local_dof_indices);
          ad_helper.register_dof_values(locally_relevant_solution,
                                        local_dof_indices);
          const std::vector<ADNumberType> &dof_values_ad =
            ad_helper.get_sensitive_dof_values();

          rain_fall_rate_rhs.value_list(fe_values.get_quadrature_points(),
                                        rain_fall_rate_rhs_values);

          fe_values[elevation].get_function_values_from_local_dof_values(
            dof_values_ad, elevation_at_q_points);
          fe_values[elevation].get_function_values(
            locally_relevant_solution_dot, elevation_dot_at_q_points);
          fe_values[elevation].get_function_gradients_from_local_dof_values(
            dof_values_ad, elevation_grad_at_q_points);
          fe_values[water_flow_rate].get_function_values_from_local_dof_values(
            dof_values_ad, water_at_q_points);

          fe_values_at_node_points[elevation]
            .get_function_gradients_from_local_dof_values(
              dof_values_ad, elevation_grad_at_node_points);
          Assert(fe_values_at_node_points.n_quadrature_points ==
                   elevation_grad_at_node_points.size(),
                 ExcInternalError());
          for (const unsigned int j : fe_values.dof_indices())
            {
              // TODO: We need some sort of offsets here, because the
              // support points of W can coincide with support points
              // for H, but ultimately the solution coefficients are
              // never both non-zero for the same index j.
              // As an intermediate step, we've hack something up to
              // work around this, under the assumption that all H dofs
              // are enumerated before all w dofs, and that they're
              // enumerated in a way such that traveral of these two
              // subblocks in lock-step implies visiting the same support
              // point.
              const unsigned int jj = j % dofs_per_cell / 2;

              const ADNumberType                 Wj = dof_values_ad[jj];
              const Tensor<1, dim, ADNumberType> d_j =
                -elevation_grad_at_node_points[jj] /
                std::sqrt(elevation_grad_at_node_points[jj] *
                            elevation_grad_at_node_points[jj] +
                          ModelParameters::regularization_epsilon *
                            ModelParameters::regularization_epsilon);

              for (const unsigned int q : fe_values.quadrature_point_indices())
                {
                  const Tensor<1, dim> grad_Nx_j =
                    fe_values[water_flow_rate].gradient(jj, q);

                  div_Ih_d_wh_at_q_points[q] += Wj * d_j * grad_Nx_j;
                }
            }

          std::vector<ADNumberType> residual_ad(n_dependent_variables,
                                                ADNumberType(0.0));
          compute_local_residual(fe_values,
                                 elevation_at_q_points,
                                 water_at_q_points,
                                 elevation_grad_at_q_points,
                                 div_Ih_d_wh_at_q_points,
                                 elevation_dot_at_q_points,
                                 rain_fall_rate_rhs_values,
                                 cell->diameter(),
                                 time,
                                 residual_ad);

          ad_helper.register_residual_vector(residual_ad);
          ad_helper.compute_linearization(cell_matrix);

          // Assemble the local contribution to the Jacobian that accounts
          // for the time integration scheme adopted by SUNDIALS:
          // J = K + alpha K'
          for (const unsigned int q : fe_values.quadrature_point_indices())
            for (const unsigned int i : fe_values.dof_indices())
              for (const unsigned int j : fe_values.dof_indices())
                cell_matrix(i, j) += alpha * fe_values[elevation].value(i, q) *
                                     fe_values[elevation].value(j, q) *
                                     fe_values.JxW(q);

          constraints.distribute_local_to_global(cell_matrix,
                                                 local_dof_indices,
                                                 jacobian);
        }

    jacobian.compress(VectorOperation::add);
  }



  template <int dim>
  void
  StreamPowerErosionProblem<dim>::solve_with_jacobian(const VectorType &rhs,
                                                      VectorType       &dst,
                                                      const double tolerance)
  {
    TimerOutput::Scope t(computing_timer, "solve");

    // SolverControl solver_control(system_matrix.m(),
    //                              1e-6 * system_rhs.l2_norm());

    SolverControl solver_control(system_matrix.m(), tolerance * rhs.l2_norm());

    SolverGMRES<VectorType>  solver(solver_control);
    LA::MPI::PreconditionILU preconditioner_H;
    LA::MPI::PreconditionILU preconditioner_w;
    preconditioner_H.initialize(system_matrix.block(0, 0));
    preconditioner_w.initialize(system_matrix.block(1, 1));

    const auto lo_prec_H =
      linear_operator<typename LA::MPI::BlockVector::BlockType>(
        system_matrix.block(0, 0), preconditioner_H);
    const auto lo_prec_w =
      linear_operator<typename LA::MPI::BlockVector::BlockType>(
        system_matrix.block(1, 1), preconditioner_w);

    const auto preconditioner =
      block_diagonal_operator<2, LA::MPI::BlockVector>(
        std::array<LinearOperator<typename LA::MPI::BlockVector::BlockType>, 2>{
          {lo_prec_H, lo_prec_w}});

    VectorType distributed_solution(owned_partitioning, mpi_communicator);

    constraints.set_zero(distributed_solution);

    solver.solve(system_matrix,
                 distributed_solution,
                 system_rhs,
                 preconditioner);

    pcout << "   Solved in " << solver_control.last_step() << " iterations."
          << std::endl;

    constraints.distribute(distributed_solution);

    dst = distributed_solution;
  }



  template <int dim>
  void StreamPowerErosionProblem<dim>::refine_grid()
  {
    TimerOutput::Scope t(computing_timer, "refine");

    triangulation.refine_global();
  }


  template <int dim>
  class DownhillFlowPostprocessor : public DataPostprocessorVector<dim>
  {
  public:
    DownhillFlowPostprocessor()
      : DataPostprocessorVector<dim>("downhill_direction",
                                     update_values | update_gradients)
    {}

    virtual void evaluate_vector_field(
      const DataPostprocessorInputs::Vector<dim> &input_data,
      std::vector<Vector<double>> &computed_quantities) const override
    {
      AssertDimension(input_data.solution_gradients.size(),
                      computed_quantities.size());

      for (unsigned int p = 0; p < input_data.solution_gradients.size(); ++p)
        {
          AssertDimension(computed_quantities[p].size(), dim);

          const double catchment_area = 1.; // input_data.solution_values[p][1];
          for (unsigned int d = 0; d < dim; ++d)
            computed_quantities[p][d] =
              -catchment_area *
              input_data.solution_gradients[p][/* vector component= */ 0][d] /
              input_data.solution_gradients[p][0].norm();
        }
    }
  };


  template <int dim>
  void StreamPowerErosionProblem<dim>::output_results(
    const double       time,
    const VectorType  &locally_relevant_solution,
    const VectorType  &locally_relevant_solution_dot,
    const unsigned int step_number)
  {
    (void)time; // TODO: Incorporate this into the PVTU data?
    TimerOutput::Scope t(computing_timer, "output");

    const std::vector<std::string> solution_names     = {"elevation",
                                                         "water_flow_rate"};
    const std::vector<std::string> solution_dot_names = {"elevation_rate",
                                                         "water_flow_rate_dot"};
    const std::vector<DataComponentInterpretation::DataComponentInterpretation>
      data_component_interpretation = {
        DataComponentInterpretation::component_is_scalar,
        DataComponentInterpretation::component_is_scalar};

    DownhillFlowPostprocessor<dim> downhill_flow_postprocessor;
    DataOut<dim>                   data_out;
    data_out.attach_dof_handler(dof_handler);
    data_out.add_data_vector(locally_relevant_solution,
                             solution_names,
                             DataOut<dim>::type_dof_data,
                             data_component_interpretation);
    data_out.add_data_vector(
      locally_relevant_solution_dot,
      solution_dot_names,
      DataOut<dim>::type_dof_data,
      data_component_interpretation); // TODO: Limit this to the elevation
                                      // component only
    data_out.add_data_vector(locally_relevant_solution,
                             downhill_flow_postprocessor);

    Vector<float> subdomain(triangulation.n_active_cells());
    for (unsigned int i = 0; i < subdomain.size(); ++i)
      subdomain(i) = triangulation.locally_owned_subdomain();
    data_out.add_data_vector(subdomain, "subdomain");

    data_out.build_patches();

    data_out.write_vtu_with_pvtu_record(
      "./", "solution", step_number, mpi_communicator, 2);
  }



  template <int dim>
  void StreamPowerErosionProblem<dim>::run()
  {
    make_grid();
    setup_system();
    interpolate_initial_elevation();
    compute_initial_constraints();

    // Get rid of these three lines eventually (?):
    assemble_initial_waterflow_system();
    output_results(/* time= */ 0,
                   locally_relevant_solution,
                   locally_relevant_solution_dot,
                   /*cycle*/ 0);
    solve_initial_waterflow_system();
    output_results(/* time= */ 0,
                   locally_relevant_solution,
                   locally_relevant_solution_dot,
                   /*cycle*/ 1);

    const SUNDIALS::IDA<VectorType>::AdditionalData time_integrator_parameters(
      /* initial_time     */ 0.0, // Initial time is today
      /* final_time       */ 1e6, // End time is one million years from now
      /* initial_step_size*/ 10,  // Start with a time step of ten years
      /* output_period    */ 100, // Produce output every 100 years
                                  // Running parameters
      /* minimum_step_size            */ 1,
      /* maximum_order                */ 5, // Is this the time stepping order?
      /* maximum_non_linear_iterations*/ 10,
      /* ls_norm_factor               */ 0,       // What is this??
                                                  // Error parameters
      /* absolute_tolerance               */ 0.1, // 0.1 meters is a good
                                                  // absolute tolerance
      /* relative_tolerance               */ 1e-5,
      /* ignore_algebraic_terms_for_errors*/ true,
      // Initial conditions parameters
      /* const InitialConditionCorrection */
      SUNDIALS::IDA<VectorType>::AdditionalData::use_y_diff, // ???
      /* const InitialConditionCorrection */
      SUNDIALS::IDA<VectorType>::AdditionalData::use_y_diff, // ???
      /* maximum_non_linear_iterations_ic*/ 5);
    SUNDIALS::IDA<VectorType> time_integrator(time_integrator_parameters);


    time_integrator.reinit_vector = [this](VectorType &vector) {
      vector.reinit(owned_partitioning, mpi_communicator);
    };

    time_integrator.residual =
      [this](const double      time,
             const VectorType &locally_relevant_solution,
             const VectorType &locally_relevant_solution_dot,
             VectorType       &residual) {
        this->assemble_residual(time,
                                locally_relevant_solution,
                                locally_relevant_solution_dot,
                                residual);
      };

    time_integrator.setup_jacobian =
      [this](const double      time,
             const VectorType &locally_relevant_solution,
             const VectorType &locally_relevant_solution_dot,
             const double      alpha) {
        this->assemble_jacobian(time,
                                locally_relevant_solution,
                                locally_relevant_solution_dot,
                                alpha,
                                system_matrix);
      };

    time_integrator.solve_with_jacobian =
      [this](const VectorType &rhs, VectorType &dst, const double tolerance) {
        this->solve_with_jacobian(rhs, dst, tolerance);
      };

    time_integrator.output_step =
      [this](const double       time,
             const VectorType  &locally_relevant_solution,
             const VectorType  &locally_relevant_solution_dot,
             const unsigned int step_number) {
        (void)time;
        (void)locally_relevant_solution;
        this->output_results(time,
                             locally_relevant_solution,
                             locally_relevant_solution_dot,
                             step_number);
      };

    //    time_integrator.differential_components;

    // TODO: Figure out whether the following vectors need to be fully
    // distributed, ghosted, or otherwise. Set to the correct vector as computed
    // above.
    pcout << "Solve time-dependent system... " << std::flush;
    VectorType solution;
    VectorType solution_dot;
    solution.reinit(owned_partitioning, mpi_communicator);
    solution_dot.reinit(owned_partitioning, mpi_communicator);

    solution = locally_relevant_solution;

    time_integrator.solve_dae(solution, solution_dot);
    pcout << "done." << std::endl;

    computing_timer.print_summary();
  }
} // namespace Step55



int main(int argc, char *argv[])
{
  try
    {
      using namespace dealii;
      using namespace Step55;

      Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

      StreamPowerErosionProblem<2> problem;
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
