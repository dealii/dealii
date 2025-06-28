/* ------------------------------------------------------------------------
 *
 * SPDX-License-Identifier: LGPL-2.1-or-later
 * Copyright (C) 2013 - 2025 by the deal.II authors
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

#include <deal.II/base/function.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/table_handler.h>
#include <deal.II/base/utilities.h>

#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/trilinos_precondition.h>
#include <deal.II/lac/trilinos_solver.h>
#include <deal.II/lac/trilinos_sparse_matrix.h>
#include <deal.II/lac/trilinos_vector.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/vector_tools.h>

#include <iostream>

using namespace dealii;

// Test that deal.II is working with Trilinos by solving the Laplace's
// problem in 2d.
class LaplaceProblem
{
public:
  LaplaceProblem();
  void
  run();

private:
  void
  setup_system();
  void
  assemble_system();
  void
  solve();

  Triangulation<2> triangulation;
  FE_Q<2>          fe;
  DoFHandler<2>    dof_handler;

  TrilinosWrappers::SparseMatrix A;
  TrilinosWrappers::MPI::Vector  b, x;
  AffineConstraints<double>      constraints;

  TableHandler output_table;
};

LaplaceProblem::LaplaceProblem()
  : fe(1)
  , dof_handler(triangulation)
{}

void
LaplaceProblem::setup_system()
{
  dof_handler.distribute_dofs(fe);

  constraints.clear();
  DoFTools::make_zero_boundary_constraints(dof_handler, constraints);
  constraints.close();

  DynamicSparsityPattern dsp(dof_handler.n_dofs());
  DoFTools::make_sparsity_pattern(dof_handler, dsp, constraints, false);

  A.reinit(dsp);
  b.reinit(complete_index_set(dof_handler.n_dofs()));
  x.reinit(complete_index_set(dof_handler.n_dofs()));

  // some output
  output_table.add_value("cells", triangulation.n_active_cells());
  output_table.add_value("dofs", dof_handler.n_dofs());
}

void
LaplaceProblem::assemble_system()
{
  QGauss<2> quadrature_formula(2);

  FEValues<2> fe_values(fe,
                        quadrature_formula,
                        update_values | update_gradients |
                          update_quadrature_points | update_JxW_values);

  const unsigned int dofs_per_cell = fe.dofs_per_cell;
  const unsigned int n_q_points    = quadrature_formula.size();

  FullMatrix<double> cell_A(dofs_per_cell, dofs_per_cell);
  Vector<double>     cell_b(dofs_per_cell);

  std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

  DoFHandler<2>::active_cell_iterator cell = dof_handler.begin_active(),
                                      endc = dof_handler.end();

  for (; cell != endc; ++cell)
    {
      fe_values.reinit(cell);
      cell_A = 0;
      cell_b = 0;

      for (unsigned int q_point = 0; q_point < n_q_points; ++q_point)
        for (unsigned int i = 0; i < dofs_per_cell; ++i)
          {
            for (unsigned int j = 0; j < dofs_per_cell; ++j)
              {
                cell_A(i, j) += fe_values.shape_grad(i, q_point) *
                                fe_values.shape_grad(j, q_point) *
                                fe_values.JxW(q_point);
              }

            cell_b(i) +=
              fe_values.shape_value(i, q_point) * fe_values.JxW(q_point);
          }

      cell->get_dof_indices(local_dof_indices);

      constraints.distribute_local_to_global(cell_A, local_dof_indices, A);
      constraints.distribute_local_to_global(cell_b, local_dof_indices, b);
    }

  A.compress(VectorOperation::add);
  b.compress(VectorOperation::add);
}

void
LaplaceProblem::solve()
{
  SolverControl                             solver_control(1e03, 1e-03);
  TrilinosWrappers::SolverCG                cg_solver(solver_control);
  TrilinosWrappers::PreconditionBlockJacobi preconditioner;
  preconditioner.initialize(A);
  cg_solver.solve(A, x, b, preconditioner);

  TrilinosWrappers::MPI::Vector res(x);
  A.residual(res, x, b);
  AssertThrow(res.l2_norm() < 1e-3, ExcInternalError());
}

void
LaplaceProblem::run()
{
  GridGenerator::hyper_cube(triangulation, -1, 1);

  for (unsigned int c = 0; c < 5; ++c)
    {
      triangulation.refine_global(1);
      setup_system();
      assemble_system();
      solve();
    }

  // finalise output
  output_table.write_text(std::cout);
  deallog << std::endl;
}


int
main(int argc, char **argv)
{
  try
    {
      Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
      {
        LaplaceProblem problem;
        problem.run();
        deallog << "OK" << std::endl;
      }
    }

  catch (const std::exception &exc)
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
