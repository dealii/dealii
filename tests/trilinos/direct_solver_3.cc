// ---------------------------------------------------------------------
//
// Copyright (C) 2017 - 2018 by the deal.II authors
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

// Trilinos direct solvers on a 2D Poisson equation for linear elements with
// deal.II's parallel vector

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
#include <deal.II/lac/la_parallel_vector.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/sparsity_tools.h>
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

  LinearAlgebra::distributed::Vector<double> solution;
  LinearAlgebra::distributed::Vector<double> system_rhs;
  LinearAlgebra::distributed::Vector<double> system_rhs_two;
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
BoundaryValues<dim>::value(const Point<dim> &p,
                           const unsigned int /*component*/) const
{
  return -0.5 / dim * p.norm_square();
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

  constraints.clear();
  std::map<unsigned int, double> boundary_values;
  VectorTools::interpolate_boundary_values(dof_handler,
                                           0,
                                           BoundaryValues<dim>(),
                                           constraints);
  constraints.close();

  IndexSet locally_owned_dofs = dof_handler.locally_owned_dofs();
  IndexSet locally_relevant_dofs;

  DoFTools::extract_locally_relevant_dofs(dof_handler, locally_relevant_dofs);


  DynamicSparsityPattern dsp(dof_handler.n_dofs());
  DoFTools::make_sparsity_pattern(dof_handler, dsp, constraints, false);
  SparsityTools::distribute_sparsity_pattern(dsp,
                                             locally_owned_dofs,
                                             MPI_COMM_WORLD,
                                             locally_relevant_dofs);

  system_matrix.reinit(locally_owned_dofs,
                       locally_owned_dofs,
                       dsp,
                       MPI_COMM_WORLD);

  solution.reinit(locally_owned_dofs, locally_relevant_dofs, MPI_COMM_WORLD);

  system_rhs.reinit(locally_owned_dofs, locally_relevant_dofs, MPI_COMM_WORLD);

  system_rhs_two.reinit(locally_owned_dofs,
                        locally_relevant_dofs,
                        MPI_COMM_WORLD);
}


template <int dim>
void
Step4<dim>::assemble_system()
{
  QGauss<dim> quadrature_formula(fe.degree + 1);

  FEValues<dim> fe_values(fe,
                          quadrature_formula,
                          update_values | update_gradients |
                            update_quadrature_points | update_JxW_values);

  const unsigned int dofs_per_cell = fe.dofs_per_cell;
  const unsigned int n_q_points    = quadrature_formula.size();

  FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);
  Vector<double>     cell_rhs(dofs_per_cell);
  Vector<double>     cell_rhs_two(dofs_per_cell);

  std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

  typename DoFHandler<dim>::active_cell_iterator cell =
                                                   dof_handler.begin_active(),
                                                 endc = dof_handler.end();

  for (; cell != endc; ++cell)
    {
      if (cell->is_locally_owned())
        {
          fe_values.reinit(cell);
          cell_matrix  = 0;
          cell_rhs     = 0;
          cell_rhs_two = 0;

          for (unsigned int q_point = 0; q_point < n_q_points; ++q_point)
            for (unsigned int i = 0; i < dofs_per_cell; ++i)
              {
                for (unsigned int j = 0; j < dofs_per_cell; ++j)
                  cell_matrix(i, j) +=
                    (fe_values.shape_grad(i, q_point) *
                     fe_values.shape_grad(j, q_point) * fe_values.JxW(q_point));

                cell_rhs(i) += (fe_values.shape_value(i, q_point) * 1.0 *
                                fe_values.JxW(q_point));

                cell_rhs_two(i) += (fe_values.shape_value(i, q_point) * 2.0 *
                                    fe_values.JxW(q_point));
              }

          cell->get_dof_indices(local_dof_indices);
          constraints.distribute_local_to_global(cell_matrix,
                                                 local_dof_indices,
                                                 system_matrix);

          constraints.distribute_local_to_global(cell_rhs,
                                                 local_dof_indices,
                                                 system_rhs,
                                                 cell_matrix);
          constraints.distribute_local_to_global(cell_rhs_two,
                                                 local_dof_indices,
                                                 system_rhs_two,
                                                 cell_matrix);
        }
    }
  system_matrix.compress(VectorOperation::add);
  system_rhs.compress(VectorOperation::add);
  system_rhs_two.compress(VectorOperation::add);
}



template <int dim>
void
Step4<dim>::solve()
{
  SolverControl solver_control(100, 1e-12);
  // factorize matrix for direct solver
  solution = 0;

  deallog.push("DirectKLU");
  TrilinosWrappers::SolverDirect::AdditionalData data;
  data.solver_type = "Amesos_Klu";
  TrilinosWrappers::SolverDirect direct_solver(solver_control, data);
  direct_solver.initialize(system_matrix);

  // do solve 1
  direct_solver.solve(solution, system_rhs);
  deallog << "Vector norm: " << solution.l2_norm();

  constraints.distribute(solution);
  deallog << " and " << solution.l2_norm() << std::endl;
  solution.update_ghost_values();

  Vector<double> cellwise_error(triangulation.n_active_cells());
  VectorTools::integrate_difference(dof_handler,
                                    solution,
                                    BoundaryValues<dim>(),
                                    cellwise_error,
                                    QGauss<dim>(3),
                                    VectorTools::L2_norm);
  deallog << " L2 error: "
          << VectorTools::compute_global_error(triangulation,
                                               cellwise_error,
                                               VectorTools::L2_norm)
          << std::endl;

  // do solve 2 without refactorizing
  solution = 0;
  direct_solver.solve(solution, system_rhs_two);
  deallog << "Vector norm: " << solution.l2_norm();

  constraints.distribute(solution);
  deallog << " and " << solution.l2_norm() << std::endl;

  deallog.pop();
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
  Utilities::MPI::MPI_InitFinalize mpi_initialization(
    argc, argv, testing_max_num_threads());
  mpi_initlog();
  deallog << std::setprecision(10);

  Step4<2> test;
  test.run();
}
