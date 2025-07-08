// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2023 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// Tests the NonlinearSolverSelector class using an example based on the
// step-77 tutorial program. This test checks the compatibility of the
// class with MPI.

#include <deal.II/base/function.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/utilities.h>

// Included from step-40
#include <deal.II/lac/generic_linear_algebra.h>

namespace LA
{
#if defined(DEAL_II_WITH_PETSC) && \
  (!defined(DEAL_II_TRILINOS_WITH_NOX) || defined(FORCE_USE_OF_PETSC))
  using namespace dealii::LinearAlgebraPETSc;
#  define USE_PETSC_LA
#elif defined(DEAL_II_TRILINOS_WITH_NOX)
  using namespace dealii::LinearAlgebraTrilinos;
#else
#  error DEAL_II_WITH_PETSC or DEAL_II_TRILINOS_WITH_NOX required
#endif
} // namespace LA

#include <deal.II/base/conditional_ostream.h>

#include <deal.II/distributed/grid_refinement.h>
#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_refinement.h>
#include <deal.II/grid/tria.h>

#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/sparse_direct.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/sparsity_tools.h>
#include <deal.II/lac/vector.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/error_estimator.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/nonlinear.h>
#include <deal.II/numerics/solution_transfer.h>
#include <deal.II/numerics/vector_tools.h>

#include <fstream>
#include <iostream>

#include "../tests.h"

namespace MPI_nonlinear_solver_selector_test
{

  using NLSolve = NonlinearSolverSelector<LA::MPI::Vector>;

#ifndef SOLVER
#  define SOLVER NLSolve::AdditionalData::kinsol
#endif

  template <int dim>
  class MinimalSurfaceProblem
  {
  public:
    MinimalSurfaceProblem();
    void
    run();

  private:
    void
    setup_system(const bool initial_step);
    void
    solve(const LA::MPI::Vector &rhs,
          LA::MPI::Vector       &solution,
          const double           tolerance);
    void
    compute_and_factorize_jacobian(const LA::MPI::Vector &evaluation_point);
    void
    compute_residual(const LA::MPI::Vector &evaluation_point,
                     LA::MPI::Vector       &residual);

    MPI_Comm mpi_communicator;

    parallel::distributed::Triangulation<dim> triangulation;

    DoFHandler<dim> dof_handler;
    FE_Q<dim>       fe;

    IndexSet locally_owned_dofs;
    IndexSet locally_relevant_dofs;

    AffineConstraints<double> nonzero_constraints;
    AffineConstraints<double> zero_constraints;

    LA::MPI::SparseMatrix jacobian_matrix;

    LA::MPI::Vector current_solution;
  };


  template <int dim>
  MinimalSurfaceProblem<dim>::MinimalSurfaceProblem()
    : mpi_communicator(MPI_COMM_WORLD)
    , triangulation(mpi_communicator,
                    typename Triangulation<dim>::MeshSmoothing(
                      Triangulation<dim>::smoothing_on_refinement |
                      Triangulation<dim>::smoothing_on_coarsening))
    , dof_handler(triangulation)
    , fe(1)
  {}



  template <int dim>
  class BoundaryValues : public Function<dim>
  {
  public:
    virtual double
    value(const Point<dim> &p, const unsigned int component = 0) const override;
  };


  template <int dim>
  double
  BoundaryValues<dim>::value(const Point<dim> &p,
                             const unsigned int /*component*/) const
  {
    return std::sin(2 * numbers::PI * (p[0] + p[1]));
  };


  template <int dim>
  void
  MinimalSurfaceProblem<dim>::setup_system(const bool initial_step)
  {
    if (initial_step)
      {
        dof_handler.distribute_dofs(fe);

        locally_owned_dofs = dof_handler.locally_owned_dofs();
        locally_relevant_dofs =
          DoFTools::extract_locally_relevant_dofs(dof_handler);

        current_solution.reinit(locally_owned_dofs, mpi_communicator);

        {
          nonzero_constraints.clear();
          nonzero_constraints.reinit(locally_owned_dofs, locally_relevant_dofs);
          DoFTools::make_hanging_node_constraints(dof_handler,
                                                  nonzero_constraints);

          nonzero_constraints.close();

          nonzero_constraints.distribute(current_solution);

          std::map<types::global_dof_index, double> boundary_values;
          VectorTools::interpolate_boundary_values(dof_handler,
                                                   0,
                                                   BoundaryValues<dim>(),
                                                   boundary_values);

          for (const auto &boundary_value : boundary_values)
            current_solution(boundary_value.first) = boundary_value.second;

          current_solution.compress(VectorOperation::insert);
        }

        {
          zero_constraints.clear();
          zero_constraints.reinit(locally_owned_dofs, locally_relevant_dofs);
          DoFTools::make_hanging_node_constraints(dof_handler,
                                                  zero_constraints);
          VectorTools::interpolate_boundary_values(
            dof_handler, 0, Functions::ZeroFunction<dim>(), zero_constraints);
        }
        zero_constraints.close();
      }

    DynamicSparsityPattern dsp(locally_relevant_dofs);

    DoFTools::make_sparsity_pattern(dof_handler, dsp, zero_constraints, false);

    SparsityTools::distribute_sparsity_pattern(dsp,
                                               dof_handler.locally_owned_dofs(),
                                               mpi_communicator,
                                               locally_relevant_dofs);

    jacobian_matrix.reinit(locally_owned_dofs,
                           locally_owned_dofs,
                           dsp,
                           mpi_communicator);
  }



  template <int dim>
  void
  MinimalSurfaceProblem<dim>::compute_and_factorize_jacobian(
    const LA::MPI::Vector &evaluation_point)
  {
    LA::MPI::Vector evaluation_point_1;
    evaluation_point_1.reinit(locally_owned_dofs,
                              locally_relevant_dofs,
                              mpi_communicator);
    evaluation_point_1 = evaluation_point;

    {
      deallog << "  Computing Jacobian matrix" << std::endl;

      const QGauss<dim> quadrature_formula(fe.degree + 1);

      jacobian_matrix = 0;

      FEValues<dim> fe_values(fe,
                              quadrature_formula,
                              update_gradients | update_quadrature_points |
                                update_JxW_values);

      const unsigned int dofs_per_cell = fe.n_dofs_per_cell();
      const unsigned int n_q_points    = quadrature_formula.size();

      FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);

      std::vector<Tensor<1, dim>> evaluation_point_gradients(n_q_points);

      std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

      for (const auto &cell : dof_handler.active_cell_iterators())
        {
          if (cell->is_locally_owned())
            {
              cell_matrix = 0.;

              fe_values.reinit(cell);

              fe_values.get_function_gradients(evaluation_point_1,
                                               evaluation_point_gradients);

              for (unsigned int q = 0; q < n_q_points; ++q)
                {
                  const double coeff =
                    1.0 / std::sqrt(1 + evaluation_point_gradients[q] *
                                          evaluation_point_gradients[q]);

                  for (unsigned int i = 0; i < dofs_per_cell; ++i)
                    {
                      for (unsigned int j = 0; j < dofs_per_cell; ++j)
                        {
                          cell_matrix(i, j) +=
                            (((fe_values.shape_grad(i, q) // ((\nabla \phi_i
                               * coeff                    //   * a_n
                               *
                               fe_values.shape_grad(j, q)) //   * \nabla \phi_j)
                              -                            //  -
                              (fe_values.shape_grad(i, q)  //  (\nabla \phi_i
                               * coeff * coeff * coeff     //   * a_n^3
                               *
                               (fe_values.shape_grad(j, q) //   * (\nabla \phi_j
                                *
                                evaluation_point_gradients[q]) //      * \nabla
                                                               //      u_n)
                               * evaluation_point_gradients[q])) //   * \nabla
                                                                 //   u_n)))
                             * fe_values.JxW(q));                // * dx
                        }
                    }
                }

              cell->get_dof_indices(local_dof_indices);

              zero_constraints.distribute_local_to_global(cell_matrix,
                                                          local_dof_indices,
                                                          jacobian_matrix);
            }
        }
    }

    jacobian_matrix.compress(VectorOperation::add);

    deallog << "  Factorizing Jacobian matrix" << std::endl;
  }



  template <int dim>
  void
  MinimalSurfaceProblem<dim>::compute_residual(
    const LA::MPI::Vector &evaluation_point,
    LA::MPI::Vector       &residual)
  {
    deallog << "  Computing residual vector..." << std::flush;
    residual = 0.;

    LA::MPI::Vector evaluation_point_1;
    evaluation_point_1.reinit(locally_owned_dofs,
                              locally_relevant_dofs,
                              mpi_communicator);
    evaluation_point_1 = evaluation_point;

    const QGauss<dim> quadrature_formula(fe.degree + 1);
    FEValues<dim>     fe_values(fe,
                            quadrature_formula,
                            update_gradients | update_quadrature_points |
                              update_JxW_values);

    const unsigned int dofs_per_cell = fe.n_dofs_per_cell();
    const unsigned int n_q_points    = quadrature_formula.size();

    Vector<double>              cell_residual(dofs_per_cell);
    std::vector<Tensor<1, dim>> evaluation_point_gradients(n_q_points);

    std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

    for (const auto &cell : dof_handler.active_cell_iterators())
      {
        if (cell->is_locally_owned())
          {
            cell_residual = 0.;
            fe_values.reinit(cell);

            fe_values.get_function_gradients(evaluation_point_1,
                                             evaluation_point_gradients);


            for (unsigned int q = 0; q < n_q_points; ++q)
              {
                const double coeff =
                  1.0 / std::sqrt(1 + evaluation_point_gradients[q] *
                                        evaluation_point_gradients[q]);

                for (unsigned int i = 0; i < dofs_per_cell; ++i)
                  cell_residual(i) =
                    (fe_values.shape_grad(i, q)      // \nabla \phi_i
                     * coeff                         // * a_n
                     * evaluation_point_gradients[q] // * \nabla u_n
                     * fe_values.JxW(q));            // * dx
              }

            cell->get_dof_indices(local_dof_indices);

            zero_constraints.distribute_local_to_global(cell_residual,
                                                        local_dof_indices,
                                                        residual);
          }
      }

    residual.compress(VectorOperation::add);
    zero_constraints.set_zero(residual);

    deallog << " norm=" << residual.l2_norm() << std::endl;
  }



  template <int dim>
  void
  MinimalSurfaceProblem<dim>::solve(const LA::MPI::Vector &rhs,
                                    LA::MPI::Vector       &solution,
                                    const double /*tolerance*/)
  {
    deallog << "  Solving linear system" << std::endl;

    SolverControl solver_control(dof_handler.n_dofs(), 1e-12);

#ifdef USE_PETSC_LA
    LA::SolverCG solver(solver_control, mpi_communicator);
#else
    LA::SolverCG solver(solver_control);
#endif

    LA::MPI::PreconditionAMG preconditioner;

    LA::MPI::PreconditionAMG::AdditionalData data;

#ifdef USE_PETSC_LA
    data.symmetric_operator = true;
#else
    /* Trilinos defaults are good */
#endif
    preconditioner.initialize(jacobian_matrix, data);

    check_solver_within_range(
      solver.solve(jacobian_matrix, solution, rhs, preconditioner),
      solver_control.last_step(),
      0,
      10);

    zero_constraints.distribute(solution);
  }

  template <int dim>
  void
  MinimalSurfaceProblem<dim>::run()
  {
    deallog << "Running with "
#ifdef USE_PETSC_LA
            << "PETSc"
#else
            << "Trilinos"
#endif
            << " on " << Utilities::MPI::n_mpi_processes(mpi_communicator)
            << " MPI rank(s)..." << std::endl;

    GridGenerator::hyper_ball(triangulation);
    triangulation.refine_global(4);

    const bool initial_step = true;

    setup_system(initial_step);

    const double target_tolerance = 1e-3;
    deallog << "  Target_tolerance: " << target_tolerance << std::endl
            << std::endl;

    typename NLSolve::AdditionalData additional_data;
    additional_data.function_tolerance = target_tolerance;
    additional_data.solver_type        = SOLVER;

    NLSolve nonlinear_solver(additional_data, mpi_communicator);

    nonlinear_solver.reinit_vector = [&](LA::MPI::Vector &x) {
      x.reinit(locally_owned_dofs, mpi_communicator);
    };

    nonlinear_solver.residual = [&](const LA::MPI::Vector &evaluation_point,
                                    LA::MPI::Vector       &residual) {
      compute_residual(evaluation_point, residual);
    };

    nonlinear_solver.setup_jacobian = [&](const LA::MPI::Vector &current_u) {
      compute_and_factorize_jacobian(current_u);
    };

    nonlinear_solver.solve_with_jacobian = [&](const LA::MPI::Vector &rhs,
                                               LA::MPI::Vector       &dst,
                                               const double tolerance) {
      solve(rhs, dst, tolerance);
    };

    nonlinear_solver.solve(current_solution);

    deallog << std::endl;
  }
} // namespace MPI_nonlinear_solver_selector_test


int
main(int argc, char *argv[])
{
  initlog();

  using namespace MPI_nonlinear_solver_selector_test;

  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

  MinimalSurfaceProblem<2> laplace_problem_2d;
  laplace_problem_2d.run();
}
