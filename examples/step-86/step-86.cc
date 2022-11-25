/* ---------------------------------------------------------------------
 *
 * Copyright (C) 2021 by the deal.II authors
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
 * Author: Stefano Zampini, King Abdullah University of Science and Technology,
 2022
 * Based on step-15 by Sven Wetterauer, University of Heidelberg, 2012,
 * and on step-77 by Wolfgang Bangerth, Colorado State University, 2021.
 */


// @sect3{Include files}

// This program starts out like most others with well known include
// files. Compared to the step-15 and step-77 programs from which most
// of what we do here is copied (including this comment), the only difference
// is the include of the header files from which we import the needed PETSc
// classes:

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function.h>
#include <deal.II/base/timer.h>
#include <deal.II/base/utilities.h>

#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>

#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_refinement.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/fe_q.h>

#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/error_estimator.h>
#include <deal.II/numerics/solution_transfer.h>

#include <deal.II/lac/petsc_snes.h>
#include <deal.II/lac/petsc_precondition.h>
#include <deal.II/lac/petsc_sparse_matrix.h>
#include <deal.II/lac/petsc_vector.h>

#include <fstream>
#include <iostream>


namespace Step86
{
  using namespace dealii;

  // Shortcuts for the PETSc types we will use throughout this tutorial
  using VectorType         = PETScWrappers::MPI::Vector;
  using MatrixType         = PETScWrappers::MPI::SparseMatrix;
  using PreconditionerType = PETScWrappers::PreconditionLU;
  using NonlinearSolver =
    PETScWrappers::NonlinearSolver<VectorType, MatrixType>;

  // @sect3{The <code>MinimalSurfaceProblem</code> class template}

  // The main class of this program is essentially a copy of the one
  // in step-15 and step-77. This class does, however, split the computation of
  // the Jacobian (system) matrix (and its factorization using a direct solver)
  // and residual into separate functions for the reasons outlined in the
  // introduction.
  //
  template <int dim>
  class MinimalSurfaceProblem
  {
  public:
    MinimalSurfaceProblem();
    void run();

  private:
    void setup_system(const bool initial_step);
    void solve(const VectorType &rhs, VectorType &solution);
    void refine_mesh();
    void output_results(const unsigned int refinement_cycle);
    void set_boundary_values();
    void set_boundary_values(VectorType &evaluation_point);
    void compute_and_factorize_jacobian(const VectorType &evaluation_point);
    void compute_jacobian(const VectorType &evaluation_point);
    void compute_residual(const VectorType &evaluation_point,
                          VectorType &      residual);

    Triangulation<dim> triangulation;

    DoFHandler<dim> dof_handler;
    FE_Q<dim>       fe;

    AffineConstraints<double> hanging_node_constraints;

    MatrixType         jacobian_matrix;
    PreconditionerType jacobian_matrix_factorization;

    VectorType current_solution;
    VectorType work;

    TimerOutput computing_timer;
  };



  // @sect3{Boundary condition}

  // The classes implementing boundary values are a copy from step-15:
  template <int dim>
  class BoundaryValues : public Function<dim>
  {
  public:
    virtual double value(const Point<dim> & p,
                         const unsigned int component = 0) const override;
  };


  template <int dim>
  double BoundaryValues<dim>::value(const Point<dim> &p,
                                    const unsigned int /*component*/) const
  {
    return std::sin(2 * numbers::PI * (p[0] + p[1]));
  }


  // @sect3{The <code>MinimalSurfaceProblem</code> class implementation}

  // @sect4{Constructor and set up functions}

  // The following few functions are also essentially copies of what
  // step-15 and step-77 already do, and so there is little to discuss.
  // The only difference is in using PETSc vectors and matrices.
  template <int dim>
  MinimalSurfaceProblem<dim>::MinimalSurfaceProblem()
    : dof_handler(triangulation)
    , fe(1)
    , computing_timer(std::cout, TimerOutput::never, TimerOutput::wall_times)
  {}



  template <int dim>
  void MinimalSurfaceProblem<dim>::setup_system(const bool initial_step)
  {
    TimerOutput::Scope t(computing_timer, "set up");

    if (initial_step)
      {
        dof_handler.distribute_dofs(fe);
        current_solution.reinit(MPI_COMM_SELF,
                                dof_handler.n_dofs(),
                                dof_handler.n_dofs());

        hanging_node_constraints.clear();
        DoFTools::make_hanging_node_constraints(dof_handler,
                                                hanging_node_constraints);
        hanging_node_constraints.close();
      }

    DynamicSparsityPattern dsp(dof_handler.n_dofs());
    DoFTools::make_sparsity_pattern(dof_handler, dsp);

    hanging_node_constraints.condense(dsp);

    const std::vector<IndexSet> locally_owned_dofs_per_proc =
      DoFTools::locally_owned_dofs_per_subdomain(dof_handler);
    const IndexSet locally_owned_dofs = locally_owned_dofs_per_proc[0];

    jacobian_matrix.reinit(locally_owned_dofs,
                           locally_owned_dofs,
                           dsp,
                           MPI_COMM_SELF);
  }



  // @sect4{Assembling and factorizing the Jacobian matrix}

  // The only difference with step-77, is that here we do not
  // factor the jacobian matrix, but we only associate it with a
  // PETSc preconditioner. An explicit call to the setup of the
  // preconditioner is not needed since PETSc will do it for us
  // right before using it for the first time.
  // Hardcoding the factorization at jacobian setup time has the
  // big disadvantage that if for some reason we want to change the
  // preconditioner type at command line, we will waste computational
  // resources by constructing the factors that will be then thrown away.
  template <int dim>
  void MinimalSurfaceProblem<dim>::compute_and_factorize_jacobian(
    const VectorType &evaluation_point)
  {
    compute_jacobian(evaluation_point);
    jacobian_matrix_factorization.initialize(jacobian_matrix);
  }

  template <int dim>
  void MinimalSurfaceProblem<dim>::compute_jacobian(
    const VectorType &evaluation_point)
  {
    TimerOutput::Scope t(computing_timer, "assembling the Jacobian");

    std::cout << "  Computing Jacobian matrix" << std::endl;
    const QGauss<dim> quadrature_formula(fe.degree + 1);

    // Enforce boundary values
    // This should not be needed during our Newton solve but we do it
    // anyway to mirror what is done in the compute_residual function
    work = evaluation_point;
    set_boundary_values(work);

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
        cell_matrix = 0;

        fe_values.reinit(cell);

        fe_values.get_function_gradients(work, evaluation_point_gradients);

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
        hanging_node_constraints.distribute_local_to_global(cell_matrix,
                                                            local_dof_indices,
                                                            jacobian_matrix);
      }
    jacobian_matrix.compress(VectorOperation::add);

    std::map<types::global_dof_index, double> boundary_values;
    VectorTools::interpolate_boundary_values(dof_handler,
                                             0,
                                             Functions::ZeroFunction<dim>(),
                                             boundary_values);
    MatrixTools::apply_boundary_values(boundary_values,
                                       jacobian_matrix,
                                       work, // dummy vectors
                                       work);
  }


  // @sect4{Computing the residual vector}

  // The following function is the same as in step-77, except that we use
  // PETSc vectors for input and output.
  template <int dim>
  void MinimalSurfaceProblem<dim>::compute_residual(
    const VectorType &evaluation_point,
    VectorType &      residual)
  {
    TimerOutput::Scope t(computing_timer, "assembling the residual");

    std::cout << "  Computing residual vector " << std::endl;
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

    // This should not be needed during our Newton solve but we do it
    // anyway, since this can be called when testing the Jacobian
    // exactness
    work = evaluation_point;
    set_boundary_values(work);

    residual = 0;

    for (const auto &cell : dof_handler.active_cell_iterators())
      {
        cell_residual = 0;
        fe_values.reinit(cell);

        fe_values.get_function_gradients(work, evaluation_point_gradients);

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
        hanging_node_constraints.distribute_local_to_global(cell_residual,
                                                            local_dof_indices,
                                                            residual);
      }
    residual.compress(VectorOperation::add);

    // We had strongly imposed boundary conditions
    // thus we zero the residual corresponding to boundary points
    for (const types::global_dof_index i :
         DoFTools::extract_boundary_dofs(dof_handler))
      residual(i) = 0;

    for (const types::global_dof_index i :
         DoFTools::extract_hanging_node_dofs(dof_handler))
      residual(i) = 0;
    residual.compress(VectorOperation::insert);
  }



  // @sect4{Solving linear systems with the Jacobian matrix}

  // Again, this is basically a verbatim copy of the function in step-77.
  // Note that PETSc nonlinear solver can also handle the solution of linear
  // system without the need for users to do it. Here we do it ourselves.
  template <int dim>
  void MinimalSurfaceProblem<dim>::solve(const VectorType &rhs,
                                         VectorType &      solution)
  {
    TimerOutput::Scope t(computing_timer, "linear system solve");
    jacobian_matrix_factorization.vmult(solution, rhs);
  }



  // @sect4{Refining the mesh, setting boundary values, and generating graphical output}

  // The following three functions are again simply copies of the ones in
  // step-15 with the exception of resizing PETSc vectors:
  template <int dim>
  void MinimalSurfaceProblem<dim>::refine_mesh()
  {
    Vector<float> estimated_error_per_cell(triangulation.n_active_cells());

    KellyErrorEstimator<dim>::estimate(
      dof_handler,
      QGauss<dim - 1>(fe.degree + 1),
      std::map<types::boundary_id, const Function<dim> *>(),
      current_solution,
      estimated_error_per_cell);

    GridRefinement::refine_and_coarsen_fixed_number(triangulation,
                                                    estimated_error_per_cell,
                                                    0.3,
                                                    0.03);

    triangulation.prepare_coarsening_and_refinement();

    SolutionTransfer<dim> solution_transfer(dof_handler);
    Vector<double>        current_solution_tmp(current_solution);
    solution_transfer.prepare_for_coarsening_and_refinement(
      current_solution_tmp);

    triangulation.execute_coarsening_and_refinement();

    dof_handler.distribute_dofs(fe);

    Vector<double> tmp(dof_handler.n_dofs());
    solution_transfer.interpolate(current_solution_tmp, tmp);

    current_solution.reinit(MPI_COMM_SELF,
                            dof_handler.n_dofs(),
                            dof_handler.n_dofs());
    current_solution = tmp;

    hanging_node_constraints.clear();

    DoFTools::make_hanging_node_constraints(dof_handler,
                                            hanging_node_constraints);
    hanging_node_constraints.close();

    hanging_node_constraints.distribute(current_solution);

    set_boundary_values();

    setup_system(/*initial_step=*/false);
  }



  template <int dim>
  void MinimalSurfaceProblem<dim>::set_boundary_values()
  {
    set_boundary_values(current_solution);
  }

  template <int dim>
  void
  MinimalSurfaceProblem<dim>::set_boundary_values(VectorType &evaluation_point)
  {
    std::map<types::global_dof_index, double> boundary_values;
    VectorTools::interpolate_boundary_values(dof_handler,
                                             0,
                                             BoundaryValues<dim>(),
                                             boundary_values);
    for (const auto &boundary_value : boundary_values)
      evaluation_point(boundary_value.first) = boundary_value.second;
    evaluation_point.compress(VectorOperation::insert);

    hanging_node_constraints.distribute(evaluation_point);
  }



  template <int dim>
  void MinimalSurfaceProblem<dim>::output_results(
    const unsigned int refinement_cycle)
  {
    TimerOutput::Scope t(computing_timer, "graphical output");

    DataOut<dim> data_out;

    data_out.attach_dof_handler(dof_handler);
    data_out.add_data_vector(current_solution, "solution");
    data_out.build_patches();

    const std::string filename =
      "solution-" + Utilities::int_to_string(refinement_cycle, 2) + ".vtu";
    std::ofstream output(filename);
    data_out.write_vtu(output);
  }



  // @sect4{The run() function and the overall logic of the program}

  // Again, this is basically a verbatim copy of the function in step-77.
  // The only difference is in how we setup the nonlinear solver.
  template <int dim>
  void MinimalSurfaceProblem<dim>::run()
  {
    GridGenerator::hyper_ball(triangulation);
    triangulation.refine_global(2);

    setup_system(/*initial_step=*/true);
    set_boundary_values();

    for (unsigned int refinement_cycle = 0; refinement_cycle < 6;
         ++refinement_cycle)
      {
        computing_timer.reset();
        std::cout << "Mesh refinement step " << refinement_cycle << std::endl;

        if (refinement_cycle != 0)
          refine_mesh();

        const double target_tolerance = 1e-3 * std::pow(0.1, refinement_cycle);
        std::cout << "  Target_tolerance: " << target_tolerance << std::endl
                  << std::endl;

        // This is where we create the nonlinear solver
        // and feed it with an object that encodes a number of additional
        // specifics (of which we only change the nonlinear tolerance we want to
        // reach; but you might want to look into what other members of the
        // PETScWrappers::NonlinearSolverData class has and play with them).
        // When using the PETSc nonlinear solver, we have two possibilites,
        // both of them are coded below for this example.
        //  - In the case with "user_control" set to true
        //    there is complete control of the linear system solution process
        //    using the setup_jacobian and solve_for_jacobian_system functions.
        //  - When "user_control" is set to false, this tutorials follows
        //    an a-la-PETSc style and only assembles the Jacobian when asked.
        //    PETSc will handle the linear system solves.
        //
        //  When using SNES, we can also check the
        //  quality of our Jacobian matrix with command line options
        //  "-snes_test_jacobian -snes_test_jacobian_view"
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

          // First we inform the nonlinear_solver about how to sample the
          // residual of our nonlinear equations
          nonlinear_solver.residual = [&](const VectorType &evaluation_point,
                                          VectorType &      residual) {
            compute_residual(evaluation_point, residual);
            return 0;
          };

          bool user_control = true;
          if (user_control)
            {
              // In this case we need to tell PETSc what to do when a
              // new Jacobian is requested. Here we do as in step-77
              // As noted above, we actually don't need to factorize
              // the matrix ourselves, since PETSc can do it for us.
              nonlinear_solver.setup_jacobian =
                [&](const VectorType &current_u) {
                  compute_and_factorize_jacobian(current_u);
                  return 0;
                };

              // We also need to tell PETSc how we solve the jacobian
              // system. In this case, this call is used within an
              // internally created preconditioner.
              // By default, the deal.II interface uses a Jacobian-free
              // Newton-Krylov solver to allow for approximate solvers.
              // Note that, by default, the Krylov solver selected
              // only runs the preconditioner.
              nonlinear_solver.solve_for_jacobian_system =
                [&](const VectorType &rhs, VectorType &dst) {
                  this->solve(rhs, dst);

                  return 0;
                };

              // When using the user_control approach, we do not
              // need to specify the Jacobian matrix to the solver.
              // We do it here because we want to be able to test
              // the correctness of the Jacobian.
              nonlinear_solver.reinit_matrices(jacobian_matrix);
            }
          else
            {
              // In this case we need to pass the matrix we will use
              // to store the Jacobian, since PETSc must return it to us
              // for the jacobian callback below. This could be also done
              // passing the matrix as a second argument of
              // nonlinear_solve.solve()
              nonlinear_solver.reinit_matrices(jacobian_matrix);

              // Last is the routine to resample the Jacobian matrix when
              // requested
              nonlinear_solver.jacobian =
                [&](const VectorType &current_u, MatrixType &A, MatrixType &P) {
                  // In this case P == jacobian_matrix;
                  compute_jacobian(current_u);
                  // Prevent from warnings
                  (void)A;
                  (void)P;
                  return 0;
                };
            }

          // We can use a monitoring routine that will be called at each
          // Newton step. Here PETSc will give us the current solution, the
          // current step, and the value of the norm of the function
          nonlinear_solver.monitor =
            [&](const VectorType &current_u, unsigned int step, double gnorm) {
              (void)current_u;
              std::cout << step << " norm=" << gnorm << std::endl;
              return 0;
            };

          // Solve the nonlinear system
          nonlinear_solver.solve(current_solution);

          // Differently from step-77, we distribute constraints after the
          // algebraic solve is done.
          hanging_node_constraints.distribute(current_solution);
        }

        // The rest is then just house-keeping: Writing data to a file for
        // visualizing, and showing a summary of the timing collected so that we
        // can interpret how long each operation has taken, how often it was
        // executed, etc:
        output_results(refinement_cycle);

        computing_timer.print_summary();

        std::cout << std::endl;
      }
  }
} // namespace Step86


int main(int argc, char **argv)
{
  try
    {
      using namespace Step86;

      Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

      MinimalSurfaceProblem<2> laplace_problem_2d;
      laplace_problem_2d.run();
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
