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

// #define DO_FULL_3D

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
    constexpr double stabilization_c            = 1; // TODO: This should be 0.1
  } // namespace ModelParameters


  template <int spacedim, typename NumberType>
  NumberType slope_from_elevation_gradient(
    const Tensor<1, spacedim, NumberType> &elevation_gradient)
  {
    return std::sqrt(elevation_gradient * elevation_gradient +
                     ModelParameters::regularization_epsilon *
                       ModelParameters::regularization_epsilon);
  }


  template <int spacedim, typename NumberType>
  Tensor<1, spacedim, NumberType> downhill_direction_from_elevation_gradient(
    const Tensor<1, spacedim, NumberType> &elevation_gradient)
  {
    return -elevation_gradient /
           slope_from_elevation_gradient(elevation_gradient);
  }


  template <int dim, int spacedim, typename NumberType>
  void compute_div_Ih_d_wh_at_q_points_from_nodal_points(
    const FEValues<dim, spacedim> &fe_values,
    const std::vector<NumberType> &local_dof_values,
    const std::vector<Tensor<1, spacedim, NumberType>>
                            &elevation_grad_at_node_points,
    std::vector<NumberType> &div_Ih_d_wh_at_q_points)
  {
    const unsigned int               w_dof = 1;
    const FEValuesExtractors::Scalar water_flow_rate(1);

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

        const unsigned int j_group =
          fe_values.get_fe().system_to_component_index(j).first;

        // TODO: if (fe.shape_function_belongs_to(water_flow_rate))
        if (j_group == w_dof)
          {
            // TODO: do this better and check for correctness.
            const unsigned int jj = (j % 2 == 0 ? j : j - 1);
            // We could write this:
            //   Assert(j % 2 == 1, ExcMessage("Wrong component"));
            //   const unsigned int jj = j - 1;
            // assuming that the H, w components are interleaved.

            const NumberType                      Wj = local_dof_values[j];
            const Tensor<1, spacedim, NumberType> d_j =
              downhill_direction_from_elevation_gradient(
                elevation_grad_at_node_points[jj]);

            for (const unsigned int q : fe_values.quadrature_point_indices())
              {
                const Tensor<1, spacedim> grad_Nx_j =
                  fe_values[water_flow_rate].gradient(j, q);

                div_Ih_d_wh_at_q_points[q] += Wj * d_j * grad_Nx_j;
              }
          }
      }
  }


  class ColoradoTopography : public Function<3>
  {
  public:
    ColoradoTopography(const MPI_Comm mpi_communicator);

    virtual double value(const Point<3> &p,
                         const unsigned int /*component*/ = 0) const override
    {
#ifdef DO_FULL_3D
      // First pull back p to longitude/latitude, expressed in degrees
      const Point<2> p_long_lat(
        std::atan2(p[1], p[0]) * 360 / (2 * numbers::PI),

        std::atan2(p[2], std::sqrt(p[0] * p[0] + p[1] * p[1])) * 360 /
          (2 * numbers::PI));

      // TODO: This is of course just a dummy elevation:
      return 4000 *
             (1 -
              std::pow(Point<2>(-105.5, 39).distance(p_long_lat), 2) /
                std::pow(Point<2>(-105.5, 39).distance(Point<2>(-109, 37)), 2));
#else
      // Convert p to long-lat by scaling to what a degree corresponds
      // to in meters
      const Point<2> p_long_lat(p[0] / 111000,
                                p[1] / 111000); // about 111km per arc degree
      return 4000 * (1 - (p_long_lat[0] + 109) / 7);
#endif

      //      return data->value(p_long_lat);
    }

  private:
    std::unique_ptr<Functions::InterpolatedUniformGridData<2>> data;
  };


  // Exact values for Colorado: -109 to -102, 37 to 41
  // But the data set is slightly larger
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



  template <int spacedim>
  class RainFallRate : public Function<spacedim>
  {
  public:
    virtual double value(const Point<spacedim> &p,
                         const unsigned int     component) const override;
  };


  template <int spacedim>
  double RainFallRate<spacedim>::value(const Point<spacedim> & /*p*/,
                                       const unsigned int /*component*/) const
  {
    return ModelParameters::rainfall_rate_p;
  }



  class StreamPowerErosionProblem
  {
  public:
    static constexpr int dim      = 2;
    static constexpr int spacedim = 3;

    StreamPowerErosionProblem();

    void run();

  private:
    void make_grid();
    void setup_system();

    void interpolate_initial_elevation();
    void compute_initial_constraints();
    void assemble_initial_waterflow_system();
    void solve_initial_waterflow_system();
    void check_conservation_for_waterflow_system(const VectorType &solution);

    template <typename NumberType>
    void compute_local_residual(
      const FEValues<dim, spacedim> &fe_values,
      const std::vector<NumberType> &local_solution_water_at_q_points,
      const std::vector<Tensor<1, spacedim, NumberType>>
                                    &local_gradient_elevation_at_q_points,
      const std::vector<NumberType> &div_Ih_d_wh_at_q_points,
      const std::vector<double>     &local_solution_dot_elevation_at_q_points,
      const std::vector<double>     &rain_fall_rate_rhs_values,
      const double                   cell_diameter,
      std::vector<NumberType>       &cell_residual,
      const bool                     debug = false) const;

    void assemble_residual(const VectorType &locally_relevant_solution,
                           const VectorType &locally_relevant_solution_dot,
                           VectorType       &residual);
    void assemble_jacobian(const VectorType &locally_relevant_solution,
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

    const FESystem<dim, spacedim>                       fe;
    parallel::distributed::Triangulation<dim, spacedim> triangulation;
    DoFHandler<dim, spacedim>                           dof_handler;

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



  StreamPowerErosionProblem::StreamPowerErosionProblem()
    : mpi_communicator(MPI_COMM_WORLD)
    // TODO
    // cf. https://dealii.org/developer/doxygen/deal.II/classFiniteElement.html
    // Caution! How one sets up the FESystem has a big
    // effect on how
    // fe.system_to_component_index(i)
    // and
    // fe.system_to_base_index(i) work.
    //
    , fe(FE_Q<dim, spacedim>(1) ^ 2) // All base elements are the same
    // , fe(FE_Q<dim, spacedim>(1),
    //      1,
    //      FE_Q<dim, spacedim>(1),
    //      1) // Two different base elements for two different fields
    , triangulation(mpi_communicator,
                    typename Triangulation<dim, spacedim>::MeshSmoothing(
                      Triangulation<dim, spacedim>::smoothing_on_refinement |
                      Triangulation<dim, spacedim>::smoothing_on_coarsening))
    , dof_handler(triangulation)
    , pcout(std::cout,
            (Utilities::MPI::this_mpi_process(mpi_communicator) == 0))
    , computing_timer(mpi_communicator,
                      pcout,
                      TimerOutput::never,
                      TimerOutput::wall_times)
  {}


  void StreamPowerErosionProblem::make_grid()
  {
    pcout << "Make grid... " << std::flush;

    GridGenerator::subdivided_hyper_rectangle(
      triangulation,
      {7, 4}, // number of subdivisions to make cells square
      Point<2>(-109., 37.),
      Point<2>(-102., 41.));
    triangulation.refine_global(4);

#ifdef DO_FULL_3D
    GridTools::transform(
      [](const Point<spacedim> &p_long_lat_degrees) {
        const Point<2> p_long_lat(p_long_lat_degrees[0] / 360 *
                                    (2 * numbers::PI),
                                  p_long_lat_degrees[1] / 360 *
                                    (2 * numbers::PI));
        const double   R = 6371000;
        return Point<spacedim>(R * std::cos(p_long_lat[1]) *
                                 std::cos(p_long_lat[0]), // X
                               R * std::cos(p_long_lat[1]) *
                                 std::sin(p_long_lat[0]),    // Y
                               R * std::sin(p_long_lat[1])); // Z
      },
      triangulation);
#else
    // 111km per degree on the earth surface
    GridTools::scale(111000., triangulation);
#endif

    pcout << "done. " << std::endl;
  }


  void StreamPowerErosionProblem::setup_system()
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


  void StreamPowerErosionProblem::interpolate_initial_elevation()
  {
    TimerOutput::Scope t(computing_timer,
                         "initial conditions: interpolate elevation");
    pcout << "Interpolate elevation... " << std::flush;

    const ColoradoTopography colorado_topography(mpi_communicator);
    const VectorFunctionFromScalarFunctionObject<spacedim> initial_values(
      [&](const Point<spacedim> &p) { return colorado_topography.value(p); },
      0,
      2);

    VectorType interpolated;
    interpolated.reinit(owned_partitioning, MPI_COMM_WORLD);
    VectorTools::interpolate(dof_handler, initial_values, interpolated);

    locally_relevant_solution = interpolated;

    pcout << "done. " << std::endl;
  }



  void StreamPowerErosionProblem::compute_initial_constraints()
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
      FEFaceValues<dim, spacedim> fe_face_values(
        fe, face_node_points, update_gradients | update_normal_vectors);

      std::vector<Tensor<1, spacedim>> elevation_gradients_at_face_nodes(
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



  void StreamPowerErosionProblem::assemble_initial_waterflow_system()
  {
    TimerOutput::Scope t(computing_timer,
                         "initial conditions: assemble waterflow system");
    pcout << "Assemble initial waterflow system... " << std::flush;

    system_matrix = 0;
    system_rhs    = 0;

    const QGauss<dim> quadrature_formula(fe.degree + 1);

    FEValues<dim, spacedim> fe_values(fe,
                                      quadrature_formula,
                                      update_values | update_gradients |
                                        update_quadrature_points |
                                        update_JxW_values);
    FEValues<dim, spacedim> fe_values_at_node_points(
      fe,
      Quadrature<dim>(
        fe.get_unit_support_points()), // could be made more efficient
      update_gradients);

    const unsigned int dofs_per_cell = fe.n_dofs_per_cell();
    const unsigned int n_q_points    = quadrature_formula.size();

    FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);
    Vector<double>     cell_rhs(dofs_per_cell);

    const RainFallRate<spacedim> rain_fall_rate_rhs;
    std::vector<double>          rain_fall_rate_rhs_values(n_q_points);

    std::vector<Tensor<1, spacedim>> elevation_grad_at_q_points(
      fe_values.n_quadrature_points);
    std::vector<Tensor<1, spacedim>> d_at_q_points(
      fe_values.n_quadrature_points);

    std::vector<Tensor<1, spacedim>> elevation_grad_at_node_points(
      dofs_per_cell);
    std::vector<Tensor<1, spacedim>> d_at_node_points(dofs_per_cell);

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
            d_at_q_points[q] = downhill_direction_from_elevation_gradient(
              elevation_grad_at_q_points[q]);

          fe_values_at_node_points[elevation].get_function_gradients(
            locally_relevant_solution, elevation_grad_at_node_points);
          for (unsigned int j = 0; j < dofs_per_cell; ++j)
            if (fe.system_to_component_index(j).first == 1) // j is water DoF
              {
                d_at_node_points[j] =
                  downhill_direction_from_elevation_gradient(
                    elevation_grad_at_node_points[j]);
              }
            else
              d_at_node_points[j] =
                numbers::signaling_nan<Tensor<1, spacedim>>();

          for (unsigned int q = 0; q < n_q_points; ++q)
            {
              for (unsigned int i = 0; i < dofs_per_cell; ++i)
                if (fe.system_to_component_index(i).first ==
                    1) // i is water DoF
                  {
                    const double phi_i_w_at_q =
                      fe_values[water_flow_rate].value(i, q);
                    const Tensor<1, spacedim> grad_phi_i_w_at_q =
                      fe_values[water_flow_rate].gradient(i, q);

                    for (unsigned int j = 0; j < dofs_per_cell; ++j)
                      if (fe.system_to_component_index(j).first ==
                          1) // j is water DoF
                        {
                          const Tensor<1, spacedim> grad_phi_j_w_at_q =
                            fe_values[water_flow_rate].gradient(j, q);

                          cell_matrix(i, j) +=
                            ((phi_i_w_at_q +
                              ModelParameters::stabilization_c *
                                cell->diameter() *
                                (d_at_q_points[q] * grad_phi_i_w_at_q)) *
                             (d_at_node_points[j] * grad_phi_j_w_at_q) *
                             fe_values.JxW(q));

                          Assert(numbers::is_finite(cell_matrix(i, j)),
                                 ExcInternalError());
                        }

                    cell_rhs(i) +=
                      ((phi_i_w_at_q +
                        ModelParameters::stabilization_c * cell->diameter() *
                          (d_at_q_points[q] * grad_phi_i_w_at_q)) *
                       rain_fall_rate_rhs_values[q] * fe_values.JxW(q));
                    Assert(numbers::is_finite(cell_rhs(i)), ExcInternalError());
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


  void StreamPowerErosionProblem::solve_initial_waterflow_system()
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


  void StreamPowerErosionProblem::check_conservation_for_waterflow_system(
    const VectorType &solution)
  {
    TimerOutput::Scope t(computing_timer,
                         "conservation check: waterflow system");

    const QGauss<dim>     quadrature_formula_cell(fe.degree + 1);
    const QGauss<dim - 1> quadrature_formula_face(fe.degree + 1);

    FEValues<dim, spacedim>     fe_values(fe,
                                      quadrature_formula_cell,
                                      update_values | update_quadrature_points |
                                        update_JxW_values);
    FEFaceValues<dim, spacedim> fe_face_values(fe,
                                               quadrature_formula_face,
                                               update_values |
                                                 update_gradients |
                                                 update_quadrature_points |
                                                 update_normal_vectors |
                                                 update_JxW_values);

    const unsigned int n_q_points_cell = quadrature_formula_cell.size();

    const RainFallRate<spacedim> rain_fall_rate_rhs;
    std::vector<double>          rain_fall_rate_rhs_values(n_q_points_cell);

    std::vector<Tensor<1, spacedim>> elevation_grad_at_fq_points(
      fe_face_values.n_quadrature_points);
    std::vector<double> water_at_fq_points(fe_face_values.n_quadrature_points);

    const FEValuesExtractors::Scalar elevation(0);
    const FEValuesExtractors::Scalar water_flow_rate(1);

    double input_from_rain_rate   = 0.0;
    double output_from_water_rate = 0.0;

    for (const auto &cell : dof_handler.active_cell_iterators())
      if (cell->is_locally_owned())
        {
          fe_values.reinit(cell);

          rain_fall_rate_rhs.value_list(fe_values.get_quadrature_points(),
                                        rain_fall_rate_rhs_values);

          for (const unsigned int q : fe_values.quadrature_point_indices())
            {
              const double  p   = rain_fall_rate_rhs_values[q];
              const double &JxW = fe_values.JxW(q);

              input_from_rain_rate += p * JxW;
            }

          for (const auto &face : cell->face_iterators())
            if (face->at_boundary())
              {
                fe_face_values.reinit(cell, face);

                fe_face_values[elevation].get_function_gradients(
                  solution, elevation_grad_at_fq_points);
                fe_face_values[water_flow_rate].get_function_values(
                  solution, water_at_fq_points);

                for (const unsigned int f_q_point :
                     fe_face_values.quadrature_point_indices())
                  {
                    const double             &w = water_at_fq_points[f_q_point];
                    const Tensor<1, spacedim> d =
                      downhill_direction_from_elevation_gradient(
                        elevation_grad_at_fq_points[f_q_point]);
                    const Tensor<1, spacedim> &N =
                      fe_face_values.normal_vector(f_q_point);
                    const double JxW = fe_face_values.JxW(f_q_point);

                    output_from_water_rate += w * (d * N) * JxW;
                  }
              }
        }

    input_from_rain_rate =
      Utilities::MPI::sum(input_from_rain_rate, mpi_communicator);
    output_from_water_rate =
      Utilities::MPI::sum(output_from_water_rate, mpi_communicator);

    const double error_abs =
      std::abs(input_from_rain_rate - output_from_water_rate);
    const double error_rel = std::abs(error_abs / input_from_rain_rate);
    pcout << "Conservation check (water)" << std::endl
          << "   Input: " << input_from_rain_rate << std::endl
          << "   Output: " << output_from_water_rate << std::endl
          << "   Error (abs): " << error_abs << std::endl
          << "   Error (rel): " << error_rel << std::endl;

    AssertThrow(error_rel < 1e-9,
                ExcMessage("Conservation of water rate not satisfied."));
  }


  template <typename NumberType>
  void StreamPowerErosionProblem::compute_local_residual(
    const FEValues<dim, spacedim> &fe_values,
    const std::vector<NumberType> &local_solution_water_at_q_points,
    const std::vector<Tensor<1, spacedim, NumberType>>
                                  &local_gradient_elevation_at_q_points,
    const std::vector<NumberType> &div_Ih_d_wh_at_q_points,
    const std::vector<double>     &local_solution_dot_elevation_at_q_points,
    const std::vector<double>     &rain_fall_rate_rhs_values,
    const double                   cell_diameter,
    std::vector<NumberType>       &cell_residual,
    const bool                     debug) const
  {
    const FEValuesExtractors::Scalar elevation(0);
    const FEValuesExtractors::Scalar water_flow_rate(1);

    const unsigned int H_dof = 0;
    const unsigned int w_dof = 1;

    for (const unsigned int q : fe_values.quadrature_point_indices())
      {
        const NumberType &H_dot = local_solution_dot_elevation_at_q_points[q];
        const NumberType &w     = local_solution_water_at_q_points[q];
        const Tensor<1, spacedim, NumberType> &grad_H =
          local_gradient_elevation_at_q_points[q];
        const NumberType &div_Ih_d_wh = div_Ih_d_wh_at_q_points[q];
        const double      p           = rain_fall_rate_rhs_values[q];
        const double     &JxW         = fe_values.JxW(q);

        const NumberType S = slope_from_elevation_gradient(
          local_gradient_elevation_at_q_points[q]);
        const Tensor<1, spacedim, NumberType> d =
          -local_gradient_elevation_at_q_points[q] / S;

        constexpr double m  = ModelParameters::stream_power_exponent_m;
        constexpr double n  = ModelParameters::stream_power_exponent_n;
        constexpr double k  = ModelParameters::stream_power_coefficient_k;
        constexpr double Kd = ModelParameters::diffusion_coefficient_Kd;
        constexpr double c  = ModelParameters::stabilization_c;

        for (const unsigned int i : fe_values.dof_indices())
          {
            const unsigned int i_group = fe.system_to_component_index(i).first;

            if (i_group == H_dof)
              {
                const double Nx_i = fe_values[elevation].value(i, q);
                const Tensor<1, spacedim> grad_Nx_i =
                  fe_values[elevation].gradient(i, q);

                // The sign adopted here needs to be consistent with
                // the extra operation in the linearisation process
                // where we add a scaled mass-matrix term (accounting
                // for the linearisation of the elevation rate).
                cell_residual[i] +=
                  (Nx_i *
                     (H_dot + k * std::pow(std::max(w, NumberType(0.0)), m) *
                                std::pow(S, n)) +
                   grad_Nx_i * (Kd * grad_H)) *
                  JxW;

                if (debug)
                  std::cout << "RES(ELEVATION)"
                            << "\n  i (H): " << i
                            << "\n  res H: " << cell_residual[i]
                            << "\n  H_dot: " << H_dot << "\n  w: " << w
                            << "\n  S: " << S
                            << "\n  std::pow(w, m): " << std::pow(w, m)
                            << "\n  std::pow(S, n): " << std::pow(S, n)
                            << "\n  grad_H: " << grad_H << std::endl;
              }
            else if (i_group == w_dof)
              {
                const double Nx_i = fe_values[water_flow_rate].value(i, q);
                const Tensor<1, spacedim> grad_Nx_i =
                  fe_values[water_flow_rate].gradient(i, q);
                const auto stabNx_i =
                  (Nx_i + c * cell_diameter * d * grad_Nx_i);

                cell_residual[i] += (stabNx_i * (p - div_Ih_d_wh)) * JxW;

                if (debug)
                  std::cout << "RES(WATER)"
                            << "\n  i (w): " << i
                            << "\n  res W: " << cell_residual[i]
                            << "\n  stabNx_i: " << stabNx_i << "\n  p: " << p
                            << "\n  div_Ih_d_wh: " << div_Ih_d_wh << std::endl;
              }
            else
              {
                AssertThrow(i_group <= w_dof, ExcMessage("Unknown DoF group"));
              }
          }
      }
  }



  void StreamPowerErosionProblem::assemble_residual(
    const VectorType &locally_relevant_solution,
    const VectorType &locally_relevant_solution_dot,
    VectorType       &residual)
  {
    TimerOutput::Scope t(computing_timer, "linear system: assemble residual");

    residual = 0;

    const QGauss<dim> quadrature_formula(fe.degree + 1);

    FEValues<dim, spacedim> fe_values(fe,
                                      quadrature_formula,
                                      update_values | update_gradients |
                                        update_quadrature_points |
                                        update_JxW_values);
    FEValues<dim, spacedim> fe_values_at_node_points(
      fe,
      Quadrature<dim>(
        fe.get_unit_support_points()), // could be made more efficient
      update_gradients);

    const unsigned int dofs_per_cell = fe.n_dofs_per_cell();
    const unsigned int n_q_points    = quadrature_formula.size();

    const RainFallRate<spacedim> rain_fall_rate_rhs;
    std::vector<double>          rain_fall_rate_rhs_values(n_q_points);

    std::vector<double> elevation_dot_at_q_points(
      fe_values.n_quadrature_points);
    std::vector<double> water_at_q_points(fe_values.n_quadrature_points);
    std::vector<Tensor<1, spacedim>> elevation_grad_at_q_points(
      fe_values.n_quadrature_points);

    std::vector<Tensor<1, spacedim>> elevation_grad_at_node_points(
      dofs_per_cell);
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

          compute_div_Ih_d_wh_at_q_points_from_nodal_points(
            fe_values,
            local_dof_values,
            elevation_grad_at_node_points,
            div_Ih_d_wh_at_q_points);

          compute_local_residual(fe_values,
                                 water_at_q_points,
                                 elevation_grad_at_q_points,
                                 div_Ih_d_wh_at_q_points,
                                 elevation_dot_at_q_points,
                                 rain_fall_rate_rhs_values,
                                 cell->diameter(),
                                 cell_residual);

          constexpr bool debug_vec = false;

          if (debug_vec)
            {
              pcout << "Vector" << std::endl;
              for (unsigned int i = 0; i < cell_residual.size(); ++i)
                pcout << cell_residual[i] << std::endl;
            }

          for (const unsigned int i : fe_values.dof_indices())
            {
              (void)i;
              AssertIsFinite(cell_residual[i]);
            }

          constraints.distribute_local_to_global(cell_residual,
                                                 local_dof_indices,
                                                 residual);
        }

    residual.compress(VectorOperation::add);
  }



  void StreamPowerErosionProblem::assemble_jacobian(
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

    FEValues<dim, spacedim> fe_values(fe,
                                      quadrature_formula,
                                      update_values | update_gradients |
                                        update_quadrature_points |
                                        update_JxW_values);
    FEValues<dim, spacedim> fe_values_at_node_points(
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

    const RainFallRate<spacedim> rain_fall_rate_rhs;
    std::vector<double>          rain_fall_rate_rhs_values(n_q_points);

    std::vector<double> elevation_dot_at_q_points(
      fe_values.n_quadrature_points);
    std::vector<ADNumberType> water_at_q_points(fe_values.n_quadrature_points);
    std::vector<Tensor<1, spacedim, ADNumberType>> elevation_grad_at_q_points(
      fe_values.n_quadrature_points);

    std::vector<Tensor<1, spacedim, ADNumberType>>
                              elevation_grad_at_node_points(dofs_per_cell);
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

          compute_div_Ih_d_wh_at_q_points_from_nodal_points(
            fe_values,
            dof_values_ad,
            elevation_grad_at_node_points,
            div_Ih_d_wh_at_q_points);

          constexpr bool debug_res  = false;
          constexpr bool debug_mtrx = false;

          std::vector<ADNumberType> residual_ad(n_dependent_variables,
                                                ADNumberType(0.0));
          compute_local_residual(fe_values,
                                 water_at_q_points,
                                 elevation_grad_at_q_points,
                                 div_Ih_d_wh_at_q_points,
                                 elevation_dot_at_q_points,
                                 rain_fall_rate_rhs_values,
                                 cell->diameter(),
                                 residual_ad,
                                 debug_res);

          ad_helper.register_residual_vector(residual_ad);
          ad_helper.compute_linearization(cell_matrix);

          if (debug_mtrx)
            {
              pcout << std::endl << "Matrix (before)" << std::endl;
              cell_matrix.print_formatted(
                std::cout, 3, true, 0, "0.0", 1.0, 1e-12);
            }

          // Assemble the local contribution to the Jacobian that accounts
          // for the time integration scheme adopted by SUNDIALS:
          // J = K + alpha M
          for (const unsigned int q : fe_values.quadrature_point_indices())
            for (const unsigned int i : fe_values.dof_indices())
              for (const unsigned int j : fe_values.dof_indices())
                cell_matrix(i, j) += alpha * fe_values[elevation].value(i, q) *
                                     fe_values[elevation].value(j, q) *
                                     fe_values.JxW(q);

          if (debug_mtrx)
            {
              pcout << std::endl << "Matrix (after)" << std::endl;
              cell_matrix.print_formatted(
                std::cout, 3, true, 0, "0.0", 1.0, 1e-12);
            }

          for (const unsigned int i : fe_values.dof_indices())
            for (const unsigned int j : fe_values.dof_indices())
              {
                (void)i;
                (void)j;
                AssertIsFinite(cell_matrix(i, j));
              }

          constraints.distribute_local_to_global(cell_matrix,
                                                 local_dof_indices,
                                                 jacobian);
        }

    jacobian.compress(VectorOperation::add);
  }



  void StreamPowerErosionProblem::solve_with_jacobian(const VectorType &rhs,
                                                      VectorType       &dst,
                                                      const double tolerance)
  {
    (void)tolerance; // TODO. What tolerance is this? We set two in the IDA
                     // solver settings, and neither of them seem to be related
                     // in any obvious manner to the value that this takes!
    TimerOutput::Scope t(computing_timer, "solve");

    ReductionControl        solver_control(system_matrix.m(),
                                    1.e-12 * rhs.l2_norm(),
                                    1.e-6);
    SolverGMRES<VectorType> solver(solver_control);

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

    solver.solve(system_matrix, distributed_solution, rhs, preconditioner);

    pcout << "   Solved in " << solver_control.last_step() << " iterations."
          << std::endl;

    constraints.distribute(distributed_solution);

    dst = distributed_solution;
  }



  void StreamPowerErosionProblem::refine_grid()
  {
    TimerOutput::Scope t(computing_timer, "refine");

    triangulation.refine_global();
  }


  class DownhillFlowPostprocessor
    : public DataPostprocessorVector<StreamPowerErosionProblem::spacedim>
  {
  public:
    static constexpr int dim      = StreamPowerErosionProblem::dim;
    static constexpr int spacedim = StreamPowerErosionProblem::spacedim;


    DownhillFlowPostprocessor()
      : DataPostprocessorVector<spacedim>("downhill_direction",
                                          update_values | update_gradients)
    {}

    virtual void evaluate_vector_field(
      const DataPostprocessorInputs::Vector<spacedim> &input_data,
      std::vector<Vector<double>> &computed_quantities) const override
    {
      AssertDimension(input_data.solution_gradients.size(),
                      computed_quantities.size());

      for (unsigned int p = 0; p < input_data.solution_gradients.size(); ++p)
        {
          AssertDimension(computed_quantities[p].size(), spacedim);

          const double catchment_area = 1.; // input_data.solution_values[p][1];
          for (unsigned int d = 0; d < dim; ++d)
            computed_quantities[p][d] =
              -catchment_area *
              input_data.solution_gradients[p][/* vector component= */ 0][d] /
              input_data.solution_gradients[p][0].norm();
        }
    }
  };


  void StreamPowerErosionProblem::output_results(
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

    DownhillFlowPostprocessor downhill_flow_postprocessor;
    DataOut<dim, spacedim>    data_out;
    data_out.attach_dof_handler(dof_handler);
    data_out.add_data_vector(locally_relevant_solution,
                             solution_names,
                             DataOut<dim, spacedim>::type_dof_data,
                             data_component_interpretation);
    data_out.add_data_vector(
      locally_relevant_solution_dot,
      solution_dot_names,
      DataOut<dim, spacedim>::type_dof_data,
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



  void StreamPowerErosionProblem::run()
  {
    make_grid();
    setup_system();
    interpolate_initial_elevation();
    compute_initial_constraints();

    // Get rid of these three lines eventually (?):
    output_results(/* time= */ 0,
                   locally_relevant_solution,
                   locally_relevant_solution_dot,
                   /*cycle*/ 0);
    assemble_initial_waterflow_system();
    solve_initial_waterflow_system();
    check_conservation_for_waterflow_system(locally_relevant_solution);
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
      /* const InitialConditionCorrection ic_type */
      SUNDIALS::IDA<VectorType>::AdditionalData::use_y_diff, // ???
      /* const InitialConditionCorrection reset_type */
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
        pcout << "Assemble residual @ " << time << "... " << std::flush;
        this->assemble_residual(locally_relevant_solution,
                                locally_relevant_solution_dot,
                                residual);
        pcout << "done." << std::endl;
      };

    time_integrator.setup_jacobian =
      [this](const double      time,
             const VectorType &locally_relevant_solution,
             const VectorType &locally_relevant_solution_dot,
             const double      alpha) {
        pcout << "Assemble Jacobian @ " << time << "... " << std::flush;
        this->assemble_jacobian(locally_relevant_solution,
                                locally_relevant_solution_dot,
                                alpha,
                                system_matrix);
        pcout << "done." << std::endl;
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
        check_conservation_for_waterflow_system(locally_relevant_solution);
        this->output_results(time,
                             locally_relevant_solution,
                             locally_relevant_solution_dot,
                             step_number);
      };

    time_integrator.differential_components = [this]() {
      IndexSet indices(dof_handler.n_dofs());

      // We will solve for the initial conditions for the water flow, so
      // we mark these DoFs as having the correct initial value.
      // We delegate to IDA the responsibility solve for consistent initial
      // conditions for the elevation, as well as the rates for both fields.
      const FEValuesExtractors::Scalar water_flow_rate(1);
      const BlockMask water_mask = fe.block_mask(water_flow_rate);
      indices.add_indices(DoFTools::extract_dofs(dof_handler, water_mask));

      // const FEValuesExtractors::Scalar elevation(0);
      // const BlockMask elevation_mask = fe.block_mask(water_flow_rate);
      // indices.add_indices(DoFTools::extract_dofs(dof_handler,
      // elevation_mask));

      return indices;
    };

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

      StreamPowerErosionProblem problem;
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
