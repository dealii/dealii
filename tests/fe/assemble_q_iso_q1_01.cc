// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2013 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// Assemble a 1d,2d,3d Poisson problem with FE_Q_iso_Q1 elements to
// check that the sparsity pattern is computed correctly.  The problem
// is taken from step-4

#include <deal.II/base/function.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_q_iso_q1.h>
#include <deal.II/fe/fe_values.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>

#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/sparse_matrix.h>
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
  assemble_preconditioner();

  Triangulation<dim> triangulation;

  FE_Q_iso_Q1<dim> fe_precondition;
  DoFHandler<dim>  dof_handler_precondition;

  AffineConstraints<double> constraints;

  SparsityPattern      prec_pattern;
  SparseMatrix<double> preconditioner_matrix;

  Vector<double> system_rhs;
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
RightHandSide<dim>::value(const Point<dim> &p,
                          const unsigned int /*component*/) const
{
  double return_value = 0;
  for (unsigned int i = 0; i < dim; ++i)
    return_value += 4 * std::pow(p[i], 4);

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
  : fe_precondition(3)
  , dof_handler_precondition(triangulation)
{}


template <int dim>
void
Step4<dim>::make_grid()
{
  GridGenerator::hyper_cube(triangulation, -1, 1);
  triangulation.refine_global(1);
}



template <int dim>
void
Step4<dim>::setup_system()
{
  dof_handler_precondition.distribute_dofs(fe_precondition);

  constraints.clear();
  std::map<unsigned int, double> boundary_values;
  VectorTools::interpolate_boundary_values(dof_handler_precondition,
                                           0,
                                           BoundaryValues<dim>(),
                                           constraints);
  constraints.close();

  {
    DynamicSparsityPattern c_sparsity(dof_handler_precondition.n_dofs());
    DoFTools::make_sparsity_pattern(dof_handler_precondition,
                                    c_sparsity,
                                    constraints,
                                    false);
    prec_pattern.copy_from(c_sparsity);
    preconditioner_matrix.reinit(prec_pattern);
  }

  system_rhs.reinit(dof_handler_precondition.n_dofs());
}



template <int dim>
void
Step4<dim>::assemble_preconditioner()
{
  QIterated<dim> quadrature_formula(QGauss<1>(2), 3);

  const RightHandSide<dim> right_hand_side;

  FEValues<dim> fe_values(fe_precondition,
                          quadrature_formula,
                          update_values | update_gradients |
                            update_quadrature_points | update_JxW_values);

  const unsigned int dofs_per_cell = fe_precondition.dofs_per_cell;
  const unsigned int n_q_points    = quadrature_formula.size();

  FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);
  Vector<double>     cell_rhs(dofs_per_cell);

  std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

  typename DoFHandler<dim>::active_cell_iterator
    cell = dof_handler_precondition.begin_active(),
    endc = dof_handler_precondition.end();

  for (; cell != endc; ++cell)
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
          }

      cell->get_dof_indices(local_dof_indices);
      constraints.distribute_local_to_global(cell_matrix,
                                             local_dof_indices,
                                             preconditioner_matrix);
    }
  preconditioner_matrix.compress(VectorOperation::add);
}



template <int dim>
void
Step4<dim>::run()
{
  deallog.push("dim " + std::to_string(dim));

  make_grid();

  setup_system();
  assemble_preconditioner();
  deallog << "nnz: " << preconditioner_matrix.n_nonzero_elements() << std::endl;

  deallog.pop();
}


int
main(int argc, char **argv)
{
  initlog(true);

  {
    Step4<1> test;
    test.run();
  }
  {
    Step4<2> test;
    test.run();
  }
  {
    Step4<3> test;
    test.run();
  }
}
