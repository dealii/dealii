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
 */

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function.h>
#include <deal.II/base/timer.h>
#include <deal.II/base/function_lib.h>

#include <deal.II/lac/generic_linear_algebra.h>

/* #define FORCE_USE_OF_TRILINOS */

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

#include <cmath>
#include <fstream>
#include <iostream>

namespace Step55
{
  using namespace dealii;



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
    return 1;
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
    void assemble_system();
    void solve();
    void refine_grid();
    void output_results(const unsigned int cycle);

    MPI_Comm mpi_communicator;

    const FESystem<dim>                       fe;
    parallel::distributed::Triangulation<dim> triangulation;
    DoFHandler<dim>                           dof_handler;

    std::vector<IndexSet> owned_partitioning;
    std::vector<IndexSet> relevant_partitioning;

    AffineConstraints<double> constraints;

    LA::MPI::BlockSparseMatrix system_matrix;
    LA::MPI::BlockVector       locally_relevant_solution;
    LA::MPI::BlockVector       system_rhs;

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
    GridGenerator::subdivided_hyper_rectangle(triangulation,
                                              {7, 4},
                                              Point<2>(-109., 37.),
                                              Point<2>(-102., 41.));
    triangulation.refine_global(4);
  }


  template <int dim>
  void StreamPowerErosionProblem<dim>::setup_system()
  {
    TimerOutput::Scope t(computing_timer, "setup");

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
      constraints.reinit(locally_owned_dofs, locally_relevant_dofs);

      // TODO: Figure out boundary values
      //      const FEValuesExtractors::Scalar waterflow_rate(1);
      //      DoFTools::make_hanging_node_constraints(dof_handler, constraints);
      //      VectorTools::interpolate_boundary_values(dof_handler,
      //                                               0,
      //                                               ExactSolution<dim>(),
      //                                               constraints,
      //                                               fe.component_mask(waterflow_rate));

      constraints.close();
    }

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
    system_rhs.reinit(owned_partitioning, mpi_communicator);
  }


  template <int dim>
  void StreamPowerErosionProblem<dim>::interpolate_initial_elevation()
  {
    TimerOutput::Scope t(computing_timer, "interpolating initial conditions");

    const ColoradoTopography colorado_topography(mpi_communicator);
    const VectorFunctionFromScalarFunctionObject<dim> initial_values(
      [&](const Point<dim> &p) { return colorado_topography.value(p); }, 0, 2);

    LA::MPI::BlockVector interpolated;
    interpolated.reinit(owned_partitioning, MPI_COMM_WORLD);
    VectorTools::interpolate(dof_handler, initial_values, interpolated);

    locally_relevant_solution = interpolated;
  }


  template <int dim>
  void StreamPowerErosionProblem<dim>::assemble_system()
  {
    TimerOutput::Scope t(computing_timer, "assembly");

    system_matrix = 0;
    system_rhs    = 0;

    const QGauss<dim> quadrature_formula(fe.degree + 1);

    FEValues<dim> fe_values(fe,
                            quadrature_formula,
                            update_values | update_gradients |
                              update_quadrature_points | update_JxW_values);

    const unsigned int dofs_per_cell = fe.n_dofs_per_cell();
    const unsigned int n_q_points    = quadrature_formula.size();

    FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);
    Vector<double>     cell_rhs(dofs_per_cell);

    const RainFallRate<dim> rain_fall_rate_rhs;
    std::vector<double>     rain_fall_rate_rhs_values(n_q_points);

    std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);
    const FEValuesExtractors::Scalar     elevation(0);
    const FEValuesExtractors::Scalar     rain_fall_rate(1);

    for (const auto &cell : dof_handler.active_cell_iterators())
      if (cell->is_locally_owned())
        {
          cell_matrix = 0;
          cell_rhs    = 0;

          fe_values.reinit(cell);
          rain_fall_rate_rhs.value_list(fe_values.get_quadrature_points(),
                                        rain_fall_rate_rhs_values);
          for (unsigned int q = 0; q < n_q_points; ++q)
            {
              for (unsigned int i = 0; i < dofs_per_cell; ++i)
                {
                  for (unsigned int j = 0; j < dofs_per_cell; ++j)
                    {
                      cell_matrix(i, j) += (0) * fe_values.JxW(q);
                    }

                  cell_rhs(i) += fe_values[rain_fall_rate].value(i, q) *
                                 rain_fall_rate_rhs_values[q] *
                                 fe_values.JxW(q);
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
  }



  template <int dim>
  void StreamPowerErosionProblem<dim>::solve()
  {
    return;

    TimerOutput::Scope t(computing_timer, "solve");

    SolverControl solver_control(system_matrix.m(),
                                 1e-6 * system_rhs.l2_norm());

    SolverGMRES<LA::MPI::BlockVector> solver(solver_control);

    LA::MPI::BlockVector distributed_solution(owned_partitioning,
                                              mpi_communicator);

    constraints.set_zero(distributed_solution);

    solver.solve(system_matrix,
                 distributed_solution,
                 system_rhs,
                 PreconditionIdentity());

    pcout << "   Solved in " << solver_control.last_step() << " iterations."
          << std::endl;

    constraints.distribute(distributed_solution);

    locally_relevant_solution = distributed_solution;
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
      : DataPostprocessorVector<dim>("direction",
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
              -catchment_area * input_data.solution_gradients[p][0][d] /
              input_data.solution_gradients[p][0].norm();
        }
    }
  };


  template <int dim>
  void StreamPowerErosionProblem<dim>::output_results(const unsigned int cycle)
  {
    TimerOutput::Scope t(computing_timer, "output");

    const std::vector<std::string> solution_names = {"elevation",
                                                     "water_flow_rate"};
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
    data_out.add_data_vector(locally_relevant_solution,
                             downhill_flow_postprocessor);

    Vector<float> subdomain(triangulation.n_active_cells());
    for (unsigned int i = 0; i < subdomain.size(); ++i)
      subdomain(i) = triangulation.locally_owned_subdomain();
    data_out.add_data_vector(subdomain, "subdomain");

    data_out.build_patches();

    data_out.write_vtu_with_pvtu_record(
      "./", "solution", cycle, mpi_communicator, 2);
  }



  template <int dim>
  void StreamPowerErosionProblem<dim>::run()
  {
    make_grid();
    setup_system();
    interpolate_initial_elevation();

    assemble_system();
    solve();

    if (Utilities::MPI::n_mpi_processes(mpi_communicator) <= 32)
      output_results(/*cycle*/ 0);

    computing_timer.print_summary();
    computing_timer.reset();
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
