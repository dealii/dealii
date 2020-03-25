// ---------------------------------------------------------------------
//
// Copyright (C) 2015 - 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE.md at
// the top level directory of deal.II.
//
// ---------------------------------------------------------------------

// Set up a laplace problem with inhomogeneity 1 on the boundary, and run a
// cycle with local refinement to have hanging nodes.
//
// This test compares the solution vector obtained by the classical method
// of applying the constraints directly to the solution obtained by solving
// the system:
//   (C^T * A * C + Id_c) * x = C^t * (b - A * k)
// with the help of constrainted_linear_operator and
// constrained_right_hand_side

#include <deal.II/base/quadrature_lib.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_refinement.h>
#include <deal.II/grid/tria.h>

#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/constrained_linear_operator.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/linear_operator.h>
#include <deal.II/lac/packaged_operation.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/vector.h>

#include <deal.II/numerics/error_estimator.h>
#include <deal.II/numerics/vector_tools.h>

#include <iostream>

#include "../tests.h"


template <int dim>
class Step6
{
public:
  Step6();

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
  refine_grid();

  Triangulation<dim>        triangulation;
  FE_Q<dim>                 fe;
  DoFHandler<dim>           dof_handler;
  AffineConstraints<double> constraints;

  SparsityPattern      sparsity_pattern;
  SparseMatrix<double> system_matrix;
  SparseMatrix<double> system_matrix_lo;

  Vector<double> solution;
  Vector<double> solution_lo;
  Vector<double> system_rhs;
  Vector<double> system_rhs_lo;
};


template <int dim>
Step6<dim>::Step6()
  : fe(2)
  , dof_handler(triangulation)
{}


template <int dim>
void
Step6<dim>::setup_system()
{
  dof_handler.distribute_dofs(fe);

  solution.reinit(dof_handler.n_dofs());
  solution_lo.reinit(dof_handler.n_dofs());

  system_rhs.reinit(dof_handler.n_dofs());
  system_rhs_lo.reinit(dof_handler.n_dofs());

  constraints.clear();
  DoFTools::make_hanging_node_constraints(dof_handler, constraints);
  VectorTools::interpolate_boundary_values(dof_handler,
                                           0,
                                           Functions::ConstantFunction<dim>(1.),
                                           constraints);
  constraints.close();

  DynamicSparsityPattern dsp(dof_handler.n_dofs());
  DoFTools::make_sparsity_pattern(dof_handler, dsp, constraints);
  sparsity_pattern.copy_from(dsp);

  system_matrix.reinit(sparsity_pattern);
  system_matrix_lo.reinit(sparsity_pattern);
}


template <int dim>
void
Step6<dim>::assemble_system()
{
  const QGauss<dim> quadrature_formula(3);

  FEValues<dim> fe_values(fe,
                          quadrature_formula,
                          update_values | update_gradients |
                            update_quadrature_points | update_JxW_values);

  const unsigned int dofs_per_cell = fe.dofs_per_cell;
  const unsigned int n_q_points    = quadrature_formula.size();

  FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);
  Vector<double>     cell_rhs(dofs_per_cell);

  std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

  typename DoFHandler<dim>::active_cell_iterator cell =
                                                   dof_handler.begin_active(),
                                                 endc = dof_handler.end();
  for (; cell != endc; ++cell)
    {
      cell_matrix = 0;
      cell_rhs    = 0;

      fe_values.reinit(cell);

      for (unsigned int q_index = 0; q_index < n_q_points; ++q_index)
        for (unsigned int i = 0; i < dofs_per_cell; ++i)
          {
            for (unsigned int j = 0; j < dofs_per_cell; ++j)
              cell_matrix(i, j) +=
                (fe_values.shape_grad(i, q_index) *
                 fe_values.shape_grad(j, q_index) * fe_values.JxW(q_index));

            cell_rhs(i) += (fe_values.shape_value(i, q_index) * 1.0 *
                            fe_values.JxW(q_index));
          }

      cell->get_dof_indices(local_dof_indices);
      constraints.distribute_local_to_global(
        cell_matrix, cell_rhs, local_dof_indices, system_matrix, system_rhs);
      cell->distribute_local_to_global(cell_matrix, system_matrix_lo);
      cell->distribute_local_to_global(cell_rhs, system_rhs_lo);
    }
}


template <int dim>
void
Step6<dim>::solve()
{
  SolverControl solver_control(1000, 1e-8);
  SolverCG<>    solver(solver_control);

  PreconditionSSOR<> preconditioner;
  preconditioner.initialize(system_matrix, 1.2);

  check_solver_within_range(
    solver.solve(system_matrix, solution, system_rhs, preconditioner),
    solver_control.last_step(),
    10,
    60);
  constraints.distribute(solution);

  const auto A   = linear_operator(system_matrix_lo);
  const auto M   = constrained_linear_operator(constraints, A);
  const auto rhs = constrained_right_hand_side(constraints, A, system_rhs_lo);

  check_solver_within_range(solver.solve(M, solution_lo, rhs, preconditioner),
                            solver_control.last_step(),
                            10,
                            60);
  constraints.distribute(solution_lo);
}


template <int dim>
void
Step6<dim>::refine_grid()
{
  Vector<float> estimated_error_per_cell(triangulation.n_active_cells());

  KellyErrorEstimator<dim>::estimate(
    dof_handler,
    QGauss<dim - 1>(3),
    std::map<types::boundary_id, const Function<dim> *>(),
    solution,
    estimated_error_per_cell);

  GridRefinement::refine_and_coarsen_fixed_number(triangulation,
                                                  estimated_error_per_cell,
                                                  0.3,
                                                  0.03);

  triangulation.execute_coarsening_and_refinement();
}


template <int dim>
void
Step6<dim>::run()
{
  for (unsigned int cycle = 0; cycle < 3; ++cycle)
    {
      if (cycle == 0)
        {
          GridGenerator::hyper_cube(triangulation);
          triangulation.refine_global(3);
        }
      else
        refine_grid();

      setup_system();
      assemble_system();
      solve();

      Vector<double> diff = solution - solution_lo;
      if (diff.l2_norm() < 1e-10)
        deallog << "OK" << std::endl;
      else
        deallog << "ERROR! Output does not match!" << std::endl;
    }
}


int
main()
{
  initlog();
  deallog << std::setprecision(10);

  Step6<2> laplace_problem_2d;
  laplace_problem_2d.run();

  return 0;
}
