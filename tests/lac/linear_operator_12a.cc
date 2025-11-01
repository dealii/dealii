// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2017 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

// tests Trilinos direct solvers on a 2D Poisson equation for linear elements
// Note: This test is a modified version of tests/trilinos/direct_solver_2.cc

#include <deal.II/base/conditional_ostream.h>
#include <deal.II/base/function.h>
#include <deal.II/base/index_set.h>
#include <deal.II/base/utilities.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>

#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/packaged_operation.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/sparsity_tools.h>
#include <deal.II/lac/trilinos_linear_operator.h>
#include <deal.II/lac/trilinos_precondition.h>
#include <deal.II/lac/trilinos_solver.h>
#include <deal.II/lac/trilinos_sparse_matrix.h>
#include <deal.II/lac/trilinos_vector.h>
#include <deal.II/lac/vector.h>

#include <deal.II/numerics/vector_tools.h>

#include "../tests.h"



template <int dim>
class Step4
{
public:
  Step4();
  void
  run();

private:
  void
  make_grid();
  void
  setup_system();
  void
  assemble_system();
  void
  solve();

  parallel::distributed::Triangulation<dim> triangulation;
  FE_Q<dim>                                 fe;
  DoFHandler<dim>                           dof_handler;

  AffineConstraints<double> constraints;
  SparsityPattern           sparsity_pattern;

  TrilinosWrappers::SparseMatrix system_matrix;

  TrilinosWrappers::MPI::Vector solution;
  TrilinosWrappers::MPI::Vector system_rhs;
};


template <int dim>
class RightHandSide : public Function<dim>
{
public:
  RightHandSide()
    : Function<dim>()
  {}

  virtual double
  value(const Point<dim> &p, const unsigned int component = 0) const;
};


template <int dim>
class RightHandSideTwo : public Function<dim>
{
public:
  RightHandSideTwo()
    : Function<dim>()
  {}

  virtual double
  value(const Point<dim> &p, const unsigned int component = 0) const;
};



template <int dim>
class BoundaryValues : public Function<dim>
{
public:
  BoundaryValues()
    : Function<dim>()
  {}

  virtual double
  value(const Point<dim> &p, const unsigned int component = 0) const;
};



template <int dim>
double
RightHandSide<dim>::value(const Point<dim> &p,
                          const unsigned int /*component*/) const
{
  double return_value = 0;
  for (unsigned int i = 0; i < dim; ++i)
    return_value += 2 * std::pow(p[i], 2);

  return return_value;
}


template <int dim>
double
RightHandSideTwo<dim>::value(const Point<dim> &p,
                             const unsigned int /*component*/) const
{
  double return_value = 0;
  for (unsigned int i = 0; i < dim; ++i)
    return_value += 4 * std::pow(p(i), 4);

  return return_value;
}


template <int dim>
double
BoundaryValues<dim>::value(const Point<dim> &p,
                           const unsigned int /*component*/) const
{
  return p.square();
}



template <int dim>
Step4<dim>::Step4()
  : triangulation(MPI_COMM_WORLD,
                  typename Triangulation<dim>::MeshSmoothing(
                    Triangulation<dim>::smoothing_on_refinement |
                    Triangulation<dim>::smoothing_on_coarsening))
  , fe(1)
  , dof_handler(triangulation)
{}


template <int dim>
void
Step4<dim>::make_grid()
{
  GridGenerator::hyper_cube(triangulation, -1, 1);
  triangulation.refine_global(6);
}



template <int dim>
void
Step4<dim>::setup_system()
{
  dof_handler.distribute_dofs(fe);

  const IndexSet &locally_owned_dofs = dof_handler.locally_owned_dofs();
  const IndexSet  locally_relevant_dofs =
    DoFTools::extract_locally_relevant_dofs(dof_handler);

  constraints.clear();
  constraints.reinit(locally_owned_dofs, locally_relevant_dofs);
  VectorTools::interpolate_boundary_values(dof_handler,
                                           0,
                                           BoundaryValues<dim>(),
                                           constraints);
  constraints.close();

  DynamicSparsityPattern dsp(locally_relevant_dofs);
  DoFTools::make_sparsity_pattern(dof_handler, dsp, constraints, false);
  SparsityTools::distribute_sparsity_pattern(dsp,
                                             locally_owned_dofs,
                                             MPI_COMM_WORLD,
                                             locally_relevant_dofs);

  system_matrix.reinit(locally_owned_dofs,
                       locally_owned_dofs,
                       dsp,
                       MPI_COMM_WORLD);

  solution.reinit(locally_relevant_dofs, MPI_COMM_WORLD);

  system_rhs.reinit(locally_owned_dofs,
                    locally_relevant_dofs,
                    MPI_COMM_WORLD,
                    true);
}


template <int dim>
void
Step4<dim>::assemble_system()
{
  QGauss<dim> quadrature_formula(fe.degree + 1);

  const RightHandSide<dim> right_hand_side;

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
      if (cell->is_locally_owned())
        {
          fe_values.reinit(cell);
          cell_matrix = 0;
          cell_rhs    = 0;

          for (unsigned int q_point = 0; q_point < n_q_points; ++q_point)
            for (unsigned int i = 0; i < dofs_per_cell; ++i)
              {
                for (unsigned int j = 0; j < dofs_per_cell; ++j)
                  cell_matrix(i, j) +=
                    (fe_values.shape_grad(i, q_point) *
                     fe_values.shape_grad(j, q_point) * fe_values.JxW(q_point));

                cell_rhs(i) +=
                  (fe_values.shape_value(i, q_point) *
                   right_hand_side.value(fe_values.quadrature_point(q_point)) *
                   fe_values.JxW(q_point));
              }

          cell->get_dof_indices(local_dof_indices);
          constraints.distribute_local_to_global(cell_matrix,
                                                 cell_rhs,
                                                 local_dof_indices,
                                                 system_matrix,
                                                 system_rhs);
        }
    }
  system_matrix.compress(VectorOperation::add);
  system_rhs.compress(VectorOperation::add);
}



template <int dim>
void
Step4<dim>::solve()
{
  using VectorType = TrilinosWrappers::MPI::Vector;

  // Compute 'reference' solution with direct solver
  VectorType temp_solution(system_rhs);
  {
    deallog.push("DirectKLU");
    temp_solution = 0;
    TrilinosWrappers::SolverDirect::AdditionalData data;
    data.solver_type = "Amesos_Klu";
    SolverControl                  solver_control(1000, 1e-10);
    TrilinosWrappers::SolverDirect solver(solver_control, data);
    solver.solve(system_matrix, temp_solution, system_rhs);
    constraints.distribute(temp_solution);
    solution = temp_solution;
    deallog.pop();
  }

  VectorType output(system_rhs);
  {
    deallog.push("deal_II_CG_SSOR");
    output = 0;
    SolverControl                           solver_control(1000, 1e-12);
    SolverCG<TrilinosWrappers::MPI::Vector> solver(solver_control);
    TrilinosWrappers::PreconditionSSOR      preconditioner;
    preconditioner.initialize(system_matrix);
    solver.solve(system_matrix, output, system_rhs, preconditioner);
    constraints.distribute(output);
    output -= temp_solution;
    const double local_error = output.l2_norm();
    const double global_error =
      std::sqrt(Utilities::MPI::sum(local_error * local_error, MPI_COMM_WORLD));
    deallog << "Norm of error in standard solve: " << global_error << std::endl;
    deallog.pop();
  }

  {
    deallog.push("LinearOperator_deal_II_CG_SSOR");
    output = 0;
    SolverControl                           solver_control(1000, 1e-12);
    SolverCG<TrilinosWrappers::MPI::Vector> solver(solver_control);
    TrilinosWrappers::PreconditionSSOR      preconditioner;
    preconditioner.initialize(system_matrix);
    const auto lo_A     = linear_operator<VectorType>(system_matrix);
    const auto lo_A_inv = inverse_operator(lo_A, solver, preconditioner);
    output              = lo_A_inv * system_rhs;
    constraints.distribute(output);
    output -= temp_solution;
    const double local_error = output.l2_norm();
    const double global_error =
      std::sqrt(Utilities::MPI::sum(local_error * local_error, MPI_COMM_WORLD));
    deallog << "Norm of error in LinearOperator solve: " << global_error
            << std::endl;
    deallog.pop();
  }
}



template <int dim>
void
Step4<dim>::run()
{
  make_grid();
  setup_system();
  assemble_system();
  solve();
}


int
main(int argc, char **argv)
{
  initlog();

  Utilities::MPI::MPI_InitFinalize mpi_initialization(
    argc, argv, testing_max_num_threads());

  try
    {
      Step4<2> test;
      test.run();
    }
  catch (const std::exception &exc)
    {
      deallog << std::endl
              << std::endl
              << "----------------------------------------------------"
              << std::endl;
      deallog << "Exception on processing: " << std::endl
              << exc.what() << std::endl
              << "Aborting!" << std::endl
              << "----------------------------------------------------"
              << std::endl;

      return 1;
    }
  catch (...)
    {
      deallog << std::endl
              << std::endl
              << "----------------------------------------------------"
              << std::endl;
      deallog << "Unknown exception!" << std::endl
              << "Aborting!" << std::endl
              << "----------------------------------------------------"
              << std::endl;
      return 1;
    };
}
