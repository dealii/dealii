/* ------------------------------------------------------------------------
 *
 * SPDX-License-Identifier: LGPL-2.1-or-later
 * Copyright (C) 2023 - 2024 by the deal.II authors
 *
 * This file is part of the deal.II library.
 *
 * Part of the source code is dual licensed under Apache-2.0 WITH
 * LLVM-exception OR LGPL-2.1-or-later. Detailed license information
 * governing the source code and code contributions can be found in
 * LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
 *
 * ------------------------------------------------------------------------
 */

// A variation of step-77 that uses PETSc's SNES library as a
// nonlinear solver. Similar to tests/sundials/step-77.cc.

#include <deal.II/base/function.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/timer.h>
#include <deal.II/base/utilities.h>

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
#include <deal.II/lac/petsc_precondition.h>
#include <deal.II/lac/petsc_snes.h>
#include <deal.II/lac/petsc_sparse_matrix.h>
#include <deal.II/lac/petsc_vector.h>
#include <deal.II/lac/vector.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/error_estimator.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/solution_transfer.h>
#include <deal.II/numerics/vector_tools.h>

#include "../tests.h"

// The following classes are used in parallel distributed computations and
// have all already been introduced in step-40:
#include <deal.II/base/conditional_ostream.h>

#include <deal.II/distributed/grid_refinement.h>
#include <deal.II/distributed/tria.h>

#include <deal.II/grid/grid_tools.h>

#include <iostream>


namespace Step77
{
  // Before writing the main class to solve the problem, we define
  // shortcuts for the types we are going to use within this tutorial.
  using VectorType         = PETScWrappers::MPI::Vector;
  using MatrixType         = PETScWrappers::MPI::SparseMatrix;
  using PreconditionerType = PETScWrappers::PreconditionLU;
  using NonlinearSolver =
    PETScWrappers::NonlinearSolver<VectorType, MatrixType>;

  // @sect3{The <code>MinimalSurfaceProblem</code> class template}

  // The main class of this program is essentially a copy of the one
  // in step-15 and step-77. This class does, however, support parallel
  // computations, and splits the setup of the tangent linear system solve
  // and the computation of the residual into separate functions for the reasons
  // outlined in the introduction. Non-homogeneous boundary conditions are
  // handled using AffineConstraints.
  //
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
    solve(const VectorType &rhs, VectorType &solution);
    void
    refine_mesh();
    void
    compute_jacobian_and_initialize_preconditioner(
      const VectorType &evaluation_point);
    void
    compute_jacobian(const VectorType &evaluation_point);
    void
    compute_residual(const VectorType &evaluation_point, VectorType &residual);

    MPI_Comm mpi_communicator;

    parallel::distributed::Triangulation<dim> triangulation;

    DoFHandler<dim> dof_handler;
    FE_Q<dim>       fe;

    AffineConstraints<double> zero_constraints;
    AffineConstraints<double> bc_constraints;

    MatrixType         jacobian_matrix;
    PreconditionerType jacobian_matrix_factorization;

    VectorType current_solution;
    VectorType scratch_vector;
    VectorType locally_relevant_solution;
  };



  // @sect3{Boundary condition}

  // The class implementing boundary values is a copy from step-15:
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
  }


  // @sect3{The <code>MinimalSurfaceProblem</code> class implementation}

  // @sect4{Constructor and set up functions}

  // The following few functions are also essentially copies of what
  // step-15 and step-77 already do, and so there is little to discuss.
  // The only difference is in using PETSc vectors and matrices and in
  // the handling of boundary conditions.
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
  void
  MinimalSurfaceProblem<dim>::setup_system(const bool initial_step)
  {
    if (initial_step)
      {
        dof_handler.distribute_dofs(fe);
      }

    const IndexSet &locally_owned_dofs = dof_handler.locally_owned_dofs();
    const IndexSet  locally_relevant_dofs =
      DoFTools::extract_locally_relevant_dofs(dof_handler);

    // Specifically, we need two types of AffineConstraints.
    // One to handle homogeneous boundary conditions for the update step.
    zero_constraints.clear();
    zero_constraints.reinit(locally_owned_dofs, locally_relevant_dofs);
    DoFTools::make_hanging_node_constraints(dof_handler, zero_constraints);
    VectorTools::interpolate_boundary_values(dof_handler,
                                             0,
                                             Functions::ZeroFunction<dim>(),
                                             zero_constraints);
    zero_constraints.close();

    // And another one to handle non-homogeneous boundary conditions
    // when computing the residual function.
    bc_constraints.clear();
    bc_constraints.reinit(locally_owned_dofs, locally_relevant_dofs);
    DoFTools::make_hanging_node_constraints(dof_handler, bc_constraints);
    VectorTools::interpolate_boundary_values(dof_handler,
                                             0,
                                             BoundaryValues<dim>(),
                                             bc_constraints);
    bc_constraints.close();


    DynamicSparsityPattern dsp(locally_relevant_dofs);
    DoFTools::make_sparsity_pattern(dof_handler, dsp);

    zero_constraints.condense(dsp);

    SparsityTools::distribute_sparsity_pattern(dsp,
                                               locally_owned_dofs,
                                               mpi_communicator,
                                               locally_relevant_dofs);

    if (initial_step)
      current_solution.reinit(locally_owned_dofs, mpi_communicator);
    scratch_vector.reinit(locally_owned_dofs, mpi_communicator);

    jacobian_matrix.reinit(locally_owned_dofs,
                           locally_owned_dofs,
                           dsp,
                           mpi_communicator);

    locally_relevant_solution.reinit(locally_owned_dofs,
                                     locally_relevant_dofs,
                                     mpi_communicator);
  }



  // @sect4{Computing the residual vector}

  // The following function is similar to that in step-77, except that it
  // supports parallel computations.
  // The residual is assembled using homogeneous boundary conditions using
  // the AffineConstraints class since we always solve for the update step.
  // However, the `locally_relevant_solution` vector needs to satisfy
  // non-homogeneous boundary conditions and resolve hanging node constraints.
  // For doing that, we need a scratch vector since ghosted vectors are
  // read-only.
  template <int dim>
  void
  MinimalSurfaceProblem<dim>::compute_residual(
    const VectorType &evaluation_point,
    VectorType       &residual)
  {
    deallog << "  Computing residual vector " << std::endl;
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

    scratch_vector = evaluation_point;
    bc_constraints.distribute(scratch_vector);
    locally_relevant_solution = scratch_vector;

    residual = 0;

    for (const auto &cell : dof_handler.active_cell_iterators())
      {
        if (!cell->is_locally_owned())
          continue;

        cell_residual = 0;
        fe_values.reinit(cell);

        fe_values.get_function_gradients(locally_relevant_solution,
                                         evaluation_point_gradients);

        for (unsigned int q = 0; q < n_q_points; ++q)
          {
            const double coeff =
              1.0 / std::sqrt(1 + evaluation_point_gradients[q] *
                                    evaluation_point_gradients[q]);

            for (unsigned int i = 0; i < dofs_per_cell; ++i)
              cell_residual(i) +=
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
    residual.compress(VectorOperation::add);
  }


  // @sect4{Assembling and factorizing the Jacobian matrix}

  // The only difference with step-77, is that here we do not
  // factor the Jacobian matrix, but we only associate it with a
  // PETSc preconditioner. An explicit call to the setup of the
  // preconditioner is not needed since PETSc will do it for us
  // right before using it for the first time.
  // Hardcoding the factorization at Jacobian setup time has the
  // big disadvantage that if for some reason we want to change the
  // preconditioner type at command line, we will waste computational
  // resources by constructing the factors that will be then thrown away.
  template <int dim>
  void
  MinimalSurfaceProblem<dim>::compute_jacobian_and_initialize_preconditioner(
    const VectorType &evaluation_point)
  {
    compute_jacobian(evaluation_point);
    jacobian_matrix_factorization.initialize(jacobian_matrix);
  }

  // The following function is similar to that in step-77,
  // except that it supports parallel assembly.
  // The Jacobian is assembled using homogeneous boundary conditions
  // since we always solve for the update step.
  // Here we don't need to reevaluate the `locally_relevant_solution`
  // vector since SNES guaranties that the Jacobian callback is called
  // always after a residual callback.
  template <int dim>
  void
  MinimalSurfaceProblem<dim>::compute_jacobian(
    const VectorType &evaluation_point)
  {
    deallog << "  Computing Jacobian matrix" << std::endl;
    const QGauss<dim> quadrature_formula(fe.degree + 1);

    (void)evaluation_point;

    jacobian_matrix = 0;

    FEValues<dim> fe_values(fe,
                            quadrature_formula,
                            update_gradients | update_quadrature_points |
                              update_JxW_values);

    const unsigned int dofs_per_cell = fe.n_dofs_per_cell();
    const unsigned int n_q_points    = quadrature_formula.size();

    FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);

    Tensor<2, dim> identity;
    for (unsigned int i = 0; i < dim; i++)
      identity[i][i] = 1.0;

    std::vector<Tensor<1, dim>> evaluation_point_gradients(n_q_points);

    std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

    for (const auto &cell : dof_handler.active_cell_iterators())
      {
        if (!cell->is_locally_owned())
          continue;

        cell_matrix = 0;

        fe_values.reinit(cell);

        fe_values.get_function_gradients(locally_relevant_solution,
                                         evaluation_point_gradients);

        for (unsigned int q = 0; q < n_q_points; ++q)
          {
            const double coeff =
              1.0 / std::sqrt(1 + evaluation_point_gradients[q] *
                                    evaluation_point_gradients[q]);
            auto B =
              fe_values.JxW(q) * coeff *
              (identity - coeff * coeff *
                            outer_product(evaluation_point_gradients[q],
                                          evaluation_point_gradients[q]));
            for (unsigned int i = 0; i < dofs_per_cell; ++i)
              {
                for (unsigned int j = 0; j < dofs_per_cell; ++j)
                  cell_matrix(i, j) +=
                    fe_values.shape_grad(i, q) * B * fe_values.shape_grad(j, q);
              }
          }

        cell->get_dof_indices(local_dof_indices);

        zero_constraints.distribute_local_to_global(cell_matrix,
                                                    local_dof_indices,
                                                    jacobian_matrix);
      }
    jacobian_matrix.compress(VectorOperation::add);
  }



  // @sect4{Solving linear systems with the Jacobian matrix}

  // Again, this is basically a verbatim copy of the function in step-77.
  // The actual factorization of the Jacobian matrix will be performed
  // the first time `vmult` is called.
  template <int dim>
  void
  MinimalSurfaceProblem<dim>::solve(const VectorType &rhs, VectorType &solution)
  {
    jacobian_matrix_factorization.vmult(solution, rhs);
  }



  // @sect4{Refining the mesh, setting boundary values, and generating graphical
  // output}

  // The following three functions are again simply copies of the ones in
  // step-15 with the exception of resizing PETSc vectors and using parallel
  // functions:
  template <int dim>
  void
  MinimalSurfaceProblem<dim>::refine_mesh()
  {
    Vector<float> estimated_error_per_cell(triangulation.n_active_cells());

    KellyErrorEstimator<dim>::estimate(
      dof_handler,
      QGauss<dim - 1>(fe.degree + 1),
      std::map<types::boundary_id, const Function<dim> *>(),
      locally_relevant_solution,
      estimated_error_per_cell);

    parallel::distributed::GridRefinement::refine_and_coarsen_fixed_number(
      triangulation, estimated_error_per_cell, 0.3, 0.03);

    triangulation.prepare_coarsening_and_refinement();

    SolutionTransfer<dim, PETScWrappers::MPI::Vector> solution_transfer(
      dof_handler);

    PETScWrappers::MPI::Vector current_solution_tmp(locally_relevant_solution);
    solution_transfer.prepare_for_coarsening_and_refinement(
      current_solution_tmp);

    triangulation.execute_coarsening_and_refinement();

    dof_handler.distribute_dofs(fe);

    const IndexSet locally_owned_dofs = dof_handler.locally_owned_dofs();
    current_solution.reinit(locally_owned_dofs, mpi_communicator);
    solution_transfer.interpolate(current_solution);

    setup_system(/*initial_step=*/false);

    bc_constraints.distribute(current_solution);
  }



  // @sect4{The run() function and the overall logic of the program}

  // Again, this is basically a verbatim copy of the function in step-77.
  // The only differences are in how we setup the nonlinear solver and in
  // the way we handle non-homogeneous boundary conditions.
  template <int dim>
  void
  MinimalSurfaceProblem<dim>::run()
  {
    GridGenerator::hyper_ball(triangulation);
    triangulation.refine_global(2);

    setup_system(/*initial_step=*/true);

    // Here we make sure the initial guess satisfies the boundary conditions.
    bc_constraints.distribute(current_solution);

    for (unsigned int refinement_cycle = 0; refinement_cycle < 6;
         ++refinement_cycle)
      {
        deallog << "Mesh refinement step " << refinement_cycle << std::endl;

        if (refinement_cycle != 0)
          refine_mesh();

        const double target_tolerance = 1e-3 * std::pow(0.1, refinement_cycle);
        deallog << "  Target_tolerance: " << target_tolerance << std::endl
                << std::endl;

        // This is where we create the nonlinear solver
        // and feed it with an object that encodes a number of additional
        // specifics (of which we only change the nonlinear tolerance we want to
        // reach; but you might want to look into what other members of the
        // PETScWrappers::NonlinearSolverData class has and play with them).
        //
        // When using the PETSc nonlinear solver, we have two possibilities,
        // both of them are coded below for this example.
        //  - In the case with `user_control` set to true
        //    there is complete control of the linear system solution process
        //    using the `setup_jacobian` and `solve_with_jacobian` callbacks.
        //  - When `user_control` is set to false, this tutorials follows
        //    an a-la-PETSc style and only assembles the Jacobian when asked.
        //    PETSc will handle the linear system solves.
        //
        //  For additional details on these solutions processes, see
        //  PETScWrappers::NonlinearSolver.
        //
        //  When using SNES, we can also check the
        //  accuracy of our Jacobian matrix with command line options
        //  **-snes_test_jacobian -snes_test_jacobian_view**.
        //  Note that in our case the test will report a non-negligible error
        //  in Frobenius norm; however, the only nonzero rows in the
        //  differences between our Jacobian and the finite-difference Jacobian
        //  computed by PETSc will be the ones associated with boundary dofs.
        //  These differences are harmless since these dofs
        //  correspond to isolated linear equations with zero right-hand side.
        {
          PETScWrappers::NonlinearSolverData additional_data;
          additional_data.absolute_tolerance = target_tolerance;

          NonlinearSolver nonlinear_solver(additional_data);

          // Here we inform the nonlinear_solver about how to sample the
          // residual of our nonlinear equations.
          nonlinear_solver.residual = [&](const VectorType &evaluation_point,
                                          VectorType       &residual) {
            compute_residual(evaluation_point, residual);
            return 0;
          };

          bool user_control = true;
          if (user_control)
            {
              // Then we tell PETSc what to do when a
              // new Jacobian is requested. Here we do as in step-77.
              nonlinear_solver.setup_jacobian =
                [&](const VectorType &current_u) {
                  compute_jacobian_and_initialize_preconditioner(current_u);
                  return 0;
                };

              // We also need to tell PETSc how we solve the Jacobian
              // system.
              nonlinear_solver.solve_with_jacobian = [&](const VectorType &rhs,
                                                         VectorType &dst) {
                solve(rhs, dst);

                return 0;
              };

              // When using the `user_control` approach, we do not
              // need to specify the Jacobian matrix to the solver.
              // We do it here because we want to be able to test
              // the correctness of the Jacobian.
              nonlinear_solver.set_matrix(jacobian_matrix);
            }
          else
            {
              // When using the PETSc style interface, we specify the matrix we
              // want to use to construct the preconditioner and the routine to
              // resample it when requested
              // In the `jacobian` callback below, we make sure that PETSc is
              // returning to us the correct matrix.
              nonlinear_solver.set_matrix(jacobian_matrix);

              nonlinear_solver.jacobian =
                [&](const VectorType &current_u, MatrixType &, MatrixType &P) {
                  Assert(P == jacobian_matrix, ExcInternalError());
                  compute_jacobian(current_u);
                  (void)P;
                  return 0;
                };
            }

          // Solver diagnostics can be performed by using a monitoring routine
          // that will be called at each Newton step. Here PETSc will give us
          // the current solution (unused here), the current step, and the
          // value of the norm of the function.
          nonlinear_solver.monitor =
            [&](const VectorType &, unsigned int step, double gnorm) {
              deallog << step << " norm=" << gnorm << std::endl;
              return 0;
            };

          // We are now set up to solve the nonlinear system
          nonlinear_solver.solve(current_solution);

          // Differently from step-77, we apply non-homogeneous boundary
          // conditions only once, after the algebraic solve is done.
          // Note that this call is only needed since this example uses hanging
          // nodes constraints.
          bc_constraints.distribute(current_solution);
        }

        deallog << std::endl;
      }
  }
} // namespace Step77


int
main(int argc, char **argv)
{
  initlog();
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

  using namespace Step77;

  MinimalSurfaceProblem<2> laplace_problem_2d;
  laplace_problem_2d.run();
}
