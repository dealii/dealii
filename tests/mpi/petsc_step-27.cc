/* ---------------------------------------------------------------------
 *
 * Copyright (C) 2006 - 2018 by the deal.II authors
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
 */



// parallelized version of step-27 with PETSc


#include <deal.II/lac/generic_linear_algebra.h>
namespace LA
{
  using namespace dealii::LinearAlgebraPETSc;
} // namespace LA

#include <deal.II/base/conditional_ostream.h>
#include <deal.II/base/function.h>
#include <deal.II/base/logstream.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/timer.h>
#include <deal.II/base/utilities.h>

#include <deal.II/distributed/grid_refinement.h>
#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_series.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/hp/dof_handler.h>
#include <deal.II/hp/fe_values.h>
#include <deal.II/hp/refinement.h>

#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/sparsity_tools.h>
#include <deal.II/lac/vector.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/error_estimator.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/vector_tools.h>

#include <complex>
#include <fstream>
#include <iostream>

#include "../tests.h"


namespace Step27
{
  template <int dim>
  class LaplaceProblem
  {
  public:
    LaplaceProblem();
    ~LaplaceProblem();

    void
    run();

  private:
    void
    setup_system();
    void
    assemble_system();
    void
    solve();
    void
    create_coarse_grid();
    void
    estimate_smoothness(Vector<float> &smoothness_indicators);
    void
    postprocess();
    std::pair<bool, unsigned int>
    predicate(const TableIndices<dim> &indices);

    MPI_Comm mpi_communicator;

    parallel::distributed::Triangulation<dim> triangulation;

    hp::DoFHandler<dim>      dof_handler;
    hp::FECollection<dim>    fe_collection;
    hp::QCollection<dim>     quadrature_collection;
    hp::QCollection<dim - 1> face_quadrature_collection;

    hp::QCollection<dim>                    fourier_q_collection;
    std::shared_ptr<FESeries::Fourier<dim>> fourier;
    std::vector<double>                     ln_k;
    Table<dim, std::complex<double>>        fourier_coefficients;

    AffineConstraints<double> constraints;

    IndexSet locally_owned_dofs;
    IndexSet locally_relevant_dofs;

    LA::MPI::SparseMatrix system_matrix;

    LA::MPI::Vector solution;
    LA::MPI::Vector system_rhs;

    const unsigned int max_degree;

    ConditionalOStream pcout;
  };



  template <int dim>
  class RightHandSide : public Function<dim>
  {
  public:
    RightHandSide()
      : Function<dim>()
    {}

    virtual double
    value(const Point<dim> &p, const unsigned int component) const override;
  };


  template <int dim>
  double
  RightHandSide<dim>::value(const Point<dim> &p,
                            const unsigned int /*component*/) const
  {
    double product = 1;
    for (unsigned int d = 0; d < dim; ++d)
      product *= (p[d] + 1);
    return product;
  }



  template <int dim, typename T>
  void
  resize(Table<dim, T> &coeff, const unsigned int N)
  {
    TableIndices<dim> size;
    for (unsigned int d = 0; d < dim; d++)
      size[d] = N;
    coeff.reinit(size);
  }



  template <int dim>
  LaplaceProblem<dim>::LaplaceProblem()
    : mpi_communicator(MPI_COMM_WORLD)
    , triangulation(mpi_communicator)
    , dof_handler(triangulation)
    , max_degree(dim <= 2 ? 7 : 5)
    , pcout(std::cout,
            (Utilities::MPI::this_mpi_process(mpi_communicator) == 0))
  {
    for (unsigned int degree = 2; degree <= max_degree; ++degree)
      {
        fe_collection.push_back(FE_Q<dim>(degree));
        quadrature_collection.push_back(QGauss<dim>(degree + 1));
        face_quadrature_collection.push_back(QGauss<dim - 1>(degree + 1));
      }

    const unsigned int N = max_degree;

    QGauss<1>      base_quadrature(2);
    QIterated<dim> quadrature(base_quadrature, N);
    for (unsigned int i = 0; i < fe_collection.size(); i++)
      fourier_q_collection.push_back(quadrature);

    fourier = std::make_shared<FESeries::Fourier<dim>>(N,
                                                       fe_collection,
                                                       fourier_q_collection);

    resize(fourier_coefficients, N);
  }



  template <int dim>
  LaplaceProblem<dim>::~LaplaceProblem()
  {
    dof_handler.clear();
  }



  template <int dim>
  void
  LaplaceProblem<dim>::setup_system()
  {
    dof_handler.distribute_dofs(fe_collection);

    locally_owned_dofs = dof_handler.locally_owned_dofs();
    DoFTools::extract_locally_relevant_dofs(dof_handler, locally_relevant_dofs);

    solution.reinit(locally_owned_dofs,
                    locally_relevant_dofs,
                    mpi_communicator);
    system_rhs.reinit(locally_owned_dofs, mpi_communicator);

    constraints.clear();
    constraints.reinit(locally_relevant_dofs);
    DoFTools::make_hanging_node_constraints(dof_handler, constraints);
    VectorTools::interpolate_boundary_values(dof_handler,
                                             0,
                                             Functions::ZeroFunction<dim>(),
                                             constraints);
#ifdef DEBUG
    // We have not dealt with chains of constraints on ghost cells yet.
    // Thus, we are content with verifying their consistency for now.
    IndexSet locally_active_dofs;
    DoFTools::extract_locally_active_dofs(dof_handler, locally_active_dofs);
    AssertThrow(constraints.is_consistent_in_parallel(
                  dof_handler.locally_owned_dofs_per_processor(),
                  locally_active_dofs,
                  mpi_communicator,
                  /*verbose=*/true),
                ExcMessage(
                  "AffineConstraints object contains inconsistencies!"));
#endif
    constraints.close();

    DynamicSparsityPattern dsp(locally_relevant_dofs);
    DoFTools::make_sparsity_pattern(dof_handler, dsp, constraints, false);
    SparsityTools::distribute_sparsity_pattern(dsp,
                                               locally_owned_dofs,
                                               mpi_communicator,
                                               locally_relevant_dofs);

    system_matrix.reinit(locally_owned_dofs,
                         locally_owned_dofs,
                         dsp,
                         mpi_communicator);
  }



  template <int dim>
  void
  LaplaceProblem<dim>::assemble_system()
  {
    hp::FEValues<dim> hp_fe_values(fe_collection,
                                   quadrature_collection,
                                   update_values | update_gradients |
                                     update_quadrature_points |
                                     update_JxW_values);

    const RightHandSide<dim> rhs_function;

    FullMatrix<double> cell_matrix;
    Vector<double>     cell_rhs;

    std::vector<types::global_dof_index> local_dof_indices;

    for (const auto &cell : dof_handler.active_cell_iterators())
      if (cell->is_locally_owned())
        {
          const unsigned int dofs_per_cell = cell->get_fe().dofs_per_cell;

          cell_matrix.reinit(dofs_per_cell, dofs_per_cell);
          cell_matrix = 0;

          cell_rhs.reinit(dofs_per_cell);
          cell_rhs = 0;

          hp_fe_values.reinit(cell);

          const FEValues<dim> &fe_values = hp_fe_values.get_present_fe_values();

          std::vector<double> rhs_values(fe_values.n_quadrature_points);
          rhs_function.value_list(fe_values.get_quadrature_points(),
                                  rhs_values);

          for (unsigned int q_point = 0;
               q_point < fe_values.n_quadrature_points;
               ++q_point)
            for (unsigned int i = 0; i < dofs_per_cell; ++i)
              {
                for (unsigned int j = 0; j < dofs_per_cell; ++j)
                  cell_matrix(i, j) +=
                    (fe_values.shape_grad(i, q_point) * // grad phi_i(x_q)
                     fe_values.shape_grad(j, q_point) * // grad phi_j(x_q)
                     fe_values.JxW(q_point));           // dx

                cell_rhs(i) +=
                  (fe_values.shape_value(i, q_point) * // phi_i(x_q)
                   rhs_values[q_point] *               // f(x_q)
                   fe_values.JxW(q_point));            // dx
              }

          local_dof_indices.resize(dofs_per_cell);
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
  void
  LaplaceProblem<dim>::solve()
  {
    LA::MPI::Vector completely_distributed_solution(locally_owned_dofs,
                                                    mpi_communicator);

    SolverControl solver_control(system_rhs.size(),
                                 1e-8 * system_rhs.l2_norm());
    //                           ^~~~
    // Loss of precision with a factor of 1e-12 with Trilinos
    LA::SolverCG cg(solver_control, mpi_communicator);

    LA::MPI::PreconditionAMG                 preconditioner;
    LA::MPI::PreconditionAMG::AdditionalData data;
    data.symmetric_operator = true;
    preconditioner.initialize(system_matrix, data);

    check_solver_within_range(cg.solve(system_matrix,
                                       completely_distributed_solution,
                                       system_rhs,
                                       preconditioner),
                              solver_control.last_step(),
                              5,
                              40);

    pcout << "   Solved in " << solver_control.last_step() << " iterations."
          << std::endl;

    constraints.distribute(completely_distributed_solution);

    solution = completely_distributed_solution;
  }



  template <int dim>
  void
  LaplaceProblem<dim>::postprocess()
  {
    Vector<float> estimated_error_per_cell(triangulation.n_active_cells());
    KellyErrorEstimator<dim>::estimate(
      dof_handler,
      face_quadrature_collection,
      std::map<types::boundary_id, const Function<dim> *>(),
      solution,
      estimated_error_per_cell);

    Vector<float> smoothness_indicators(triangulation.n_active_cells());
    estimate_smoothness(smoothness_indicators);

    parallel::distributed::GridRefinement::refine_and_coarsen_fixed_number(
      triangulation, estimated_error_per_cell, 0.3, 0.03);

    hp::Refinement::p_adaptivity_from_relative_threshold(dof_handler,
                                                         smoothness_indicators,
                                                         0.5,
                                                         0.);
    hp::Refinement::choose_p_over_h(dof_handler);

    triangulation.execute_coarsening_and_refinement();
  }



  template <>
  void
  LaplaceProblem<2>::create_coarse_grid()
  {
    const unsigned int dim = 2;

    const std::vector<Point<2>> vertices = {
      {-1.0, -1.0}, {-0.5, -1.0}, {+0.0, -1.0}, {+0.5, -1.0}, {+1.0, -1.0}, //
      {-1.0, -0.5}, {-0.5, -0.5}, {+0.0, -0.5}, {+0.5, -0.5}, {+1.0, -0.5}, //
      {-1.0, +0.0}, {-0.5, +0.0}, {+0.5, +0.0}, {+1.0, +0.0},               //
      {-1.0, +0.5}, {-0.5, +0.5}, {+0.0, +0.5}, {+0.5, +0.5}, {+1.0, +0.5}, //
      {-1.0, +1.0}, {-0.5, +1.0}, {+0.0, +1.0}, {+0.5, +1.0}, {+1.0, +1.0}};

    const std::vector<std::array<int, GeometryInfo<dim>::vertices_per_cell>>
      cell_vertices = {{{0, 1, 5, 6}},
                       {{1, 2, 6, 7}},
                       {{2, 3, 7, 8}},
                       {{3, 4, 8, 9}},
                       {{5, 6, 10, 11}},
                       {{8, 9, 12, 13}},
                       {{10, 11, 14, 15}},
                       {{12, 13, 17, 18}},
                       {{14, 15, 19, 20}},
                       {{15, 16, 20, 21}},
                       {{16, 17, 21, 22}},
                       {{17, 18, 22, 23}}};

    const unsigned int n_cells = cell_vertices.size();

    std::vector<CellData<dim>> cells(n_cells, CellData<dim>());
    for (unsigned int i = 0; i < n_cells; ++i)
      {
        for (const unsigned int j : GeometryInfo<dim>::vertex_indices())
          cells[i].vertices[j] = cell_vertices[i][j];
        cells[i].material_id = 0;
      }

    triangulation.create_triangulation(vertices, cells, SubCellData());
    triangulation.refine_global(3);
  }



  template <int dim>
  void
  LaplaceProblem<dim>::run()
  {
    for (unsigned int cycle = 0; cycle < 5; ++cycle)
      {
        pcout << "Cycle " << cycle << ':' << std::endl;

        if (cycle == 0)
          create_coarse_grid();

        setup_system();

        pcout << "   Number of active cells      : "
              << triangulation.n_global_active_cells() << std::endl
              << "   Number of degrees of freedom: " << dof_handler.n_dofs()
              << std::endl
              << "   Number of constraints       : "
              << Utilities::MPI::sum(constraints.n_constraints(),
                                     mpi_communicator)
              << std::endl;

        assemble_system();
        solve();
        postprocess();
      }
  }



  template <int dim>
  std::pair<bool, unsigned int>
  LaplaceProblem<dim>::predicate(const TableIndices<dim> &ind)
  {
    unsigned int v = 0;
    for (unsigned int i = 0; i < dim; i++)
      v += ind[i] * ind[i];
    if (v > 0 && v < max_degree * max_degree)
      return std::make_pair(true, v);
    else
      return std::make_pair(false, v);
  }



  template <int dim>
  void
  LaplaceProblem<dim>::estimate_smoothness(Vector<float> &smoothness_indicators)
  {
    Vector<double> local_dof_values;

    for (const auto &cell : dof_handler.active_cell_iterators())
      if (cell->is_locally_owned())
        {
          local_dof_values.reinit(cell->get_fe().dofs_per_cell);
          cell->get_dof_values(solution, local_dof_values);

          fourier->calculate(local_dof_values,
                             cell->active_fe_index(),
                             fourier_coefficients);

          std::pair<std::vector<unsigned int>, std::vector<double>> res =
            FESeries::process_coefficients<dim>(
              fourier_coefficients,
              std::bind(&LaplaceProblem<dim>::predicate,
                        this,
                        std::placeholders::_1),
              VectorTools::Linfty_norm);

          Assert(res.first.size() == res.second.size(), ExcInternalError());

          if (ln_k.size() == 0)
            {
              ln_k.resize(res.first.size(), 0);
              for (unsigned int f = 0; f < ln_k.size(); f++)
                ln_k[f] =
                  std::log(2.0 * numbers::PI * std::sqrt(1. * res.first[f]));
            }

          for (double &residual_element : res.second)
            residual_element = std::log(residual_element);

          std::pair<double, double> fit =
            FESeries::linear_regression(ln_k, res.second);

          smoothness_indicators(cell->active_cell_index()) =
            -fit.first - 1. * dim / 2;
        }
  }
} // namespace Step27



int
main(int argc, char *argv[])
{
  try
    {
      using namespace Step27;

      Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

      LaplaceProblem<2> laplace_problem;
      laplace_problem.run();
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
