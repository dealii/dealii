// ---------------------------------------------------------------------
//
// Copyright (C) 2009 - 2013 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------



// this function tests the correctness of the implementation of
// inhomogeneous constraints on a nonsymmetric matrix that comes from a
// discretization of the first derivative

#include "../tests.h"

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function.h>
#include <deal.II/base/logstream.h>
#include <deal.II/base/utilities.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_refinement.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/lac/constraint_matrix.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/error_estimator.h>
#include <deal.II/lac/compressed_simple_sparsity_pattern.h>

#include <fstream>
#include <iostream>
#include <complex>

std::ofstream logfile("output");

using namespace dealii;

template <int dim>
class AdvectionProblem
{
public:
  AdvectionProblem ();
  ~AdvectionProblem ();

  void run ();

private:
  void setup_system ();
  void test_equality ();
  void assemble_reference ();
  void assemble_test_1 ();
  void assemble_test_2 ();

  Triangulation<dim>   triangulation;

  DoFHandler<dim>      dof_handler;
  FE_Q<dim>            fe;

  ConstraintMatrix     hanging_nodes_only;
  ConstraintMatrix     test_all_constraints;

  SparsityPattern      sparsity_pattern;
  SparseMatrix<double> reference_matrix;
  SparseMatrix<double> test_matrix;

  Vector<double>       reference_rhs;
  Vector<double>       test_rhs;
};



template <int dim>
class RightHandSide : public Function<dim>
{
public:
  RightHandSide () : Function<dim> () {}

  virtual double value (const Point<dim>   &p,
                        const unsigned int  component) const;
};


template <int dim>
double
RightHandSide<dim>::value (const Point<dim>   &p,
                           const unsigned int  /*component*/) const
{
  double product = 1;
  for (unsigned int d=0; d<dim; ++d)
    product *= (p[d]+1);
  return product;
}


template <int dim>
AdvectionProblem<dim>::AdvectionProblem ()
  :
  dof_handler (triangulation),
  fe (2)
{}


template <int dim>
AdvectionProblem<dim>::~AdvectionProblem ()
{
  dof_handler.clear ();
}


template <int dim>
void AdvectionProblem<dim>::setup_system ()
{
  dof_handler.distribute_dofs (fe);

  reference_rhs.reinit (dof_handler.n_dofs());
  test_rhs.reinit (dof_handler.n_dofs());

  hanging_nodes_only.clear ();
  test_all_constraints.clear ();

  // add boundary conditions as
  // inhomogeneous constraints here. just
  // take the right hand side function as
  // boundary function
  {
    std::map<types::global_dof_index,double> boundary_values;
    VectorTools::interpolate_boundary_values (dof_handler,
                                              0,
                                              RightHandSide<dim>(),
                                              boundary_values);
    std::map<types::global_dof_index,double>::const_iterator boundary_value =
      boundary_values.begin();
    for ( ; boundary_value !=boundary_values.end(); ++boundary_value)
      {
        test_all_constraints.add_line(boundary_value->first);
        test_all_constraints.set_inhomogeneity (boundary_value->first,
                                                boundary_value->second);
      }
  }
  DoFTools::make_hanging_node_constraints (dof_handler,
                                           hanging_nodes_only);
  DoFTools::make_hanging_node_constraints (dof_handler,
                                           test_all_constraints);
  hanging_nodes_only.close ();
  test_all_constraints.close ();

  CompressedSimpleSparsityPattern csp (dof_handler.n_dofs(),
                                       dof_handler.n_dofs());
  DoFTools::make_sparsity_pattern (dof_handler, csp,
                                   hanging_nodes_only, true);
  sparsity_pattern.copy_from (csp);

  reference_matrix.reinit (sparsity_pattern);
  test_matrix.reinit (sparsity_pattern);
}



// test whether we are equal with the
// standard matrix and right hand side
template <int dim>
void AdvectionProblem<dim>::test_equality ()
{
  // need to manually go through the
  // matrix, since we can have different
  // entries in constrained lines.
  for (unsigned int i=0; i<reference_matrix.m(); ++i)
    {
      SparseMatrix<double>::const_iterator reference = reference_matrix.begin(i);
      SparseMatrix<double>::iterator test = test_matrix.begin(i);
      if (test_all_constraints.is_constrained(i) == false)
        {
          for ( ; test != test_matrix.end(i); ++test, ++reference)
            test->value() -= reference->value();
        }
      else
        for ( ; test != test_matrix.end(i); ++test)
          test->value() = 0;
    }

  deallog << "  Matrix difference norm: "
          << test_matrix.frobenius_norm() << std::endl;
  Assert (test_matrix.frobenius_norm() < 1e-13, ExcInternalError());

  // same here -- Dirichlet lines will have
  // nonzero rhs, whereas we will have zero
  // rhs when using inhomogeneous
  // constraints.
  for (unsigned int i=0; i<reference_matrix.m(); ++i)
    if (test_all_constraints.is_constrained(i) == false)
      test_rhs(i) -= reference_rhs(i);
    else
      test_rhs(i) = 0;

  deallog << "  RHS difference norm: "
          << test_rhs.l2_norm() << std::endl;

  Assert (test_rhs.l2_norm() < 1e-14, ExcInternalError());
}




template <int dim>
void AdvectionProblem<dim>::assemble_reference ()
{
  reference_matrix = 0;
  reference_rhs = 0;

  QGauss<dim> quadrature_formula (3);
  FEValues<dim> fe_values (fe, quadrature_formula,
                           update_values    |  update_gradients |
                           update_quadrature_points  |  update_JxW_values);

  const RightHandSide<dim> rhs_function;
  const unsigned int dofs_per_cell = fe.dofs_per_cell;
  const unsigned int n_q_points = quadrature_formula.size();

  FullMatrix<double>   cell_matrix (dofs_per_cell, dofs_per_cell);
  Vector<double>       cell_rhs (dofs_per_cell);

  std::vector<types::global_dof_index> local_dof_indices (dofs_per_cell);
  std::vector<double>  rhs_values (n_q_points);

  typename DoFHandler<dim>::active_cell_iterator
  cell = dof_handler.begin_active(),
  endc = dof_handler.end();
  for (; cell!=endc; ++cell)
    {
      cell_matrix = 0;
      cell_rhs = 0;
      fe_values.reinit (cell);

      rhs_function.value_list (fe_values.get_quadrature_points(),
                               rhs_values);

      Tensor<1,dim> advection_direction;
      advection_direction[0] = 1;
      advection_direction[1] = 1;
      advection_direction[dim-1] = -1;

      for (unsigned int q_point=0; q_point<n_q_points; ++q_point)
        for (unsigned int i=0; i<dofs_per_cell; ++i)
          {
            for (unsigned int j=0; j<dofs_per_cell; ++j)
              cell_matrix(i,j) += (fe_values.shape_value(i,q_point) *
                                   advection_direction *
                                   fe_values.shape_grad(j,q_point) *
                                   fe_values.JxW(q_point));

            cell_rhs(i) += (fe_values.shape_value(i,q_point) *
                            rhs_values[q_point] *
                            fe_values.JxW(q_point));
          }

      local_dof_indices.resize (dofs_per_cell);
      cell->get_dof_indices (local_dof_indices);

      reference_matrix.add(local_dof_indices, cell_matrix);
      for (unsigned int i=0; i<dofs_per_cell; ++i)
        reference_rhs(local_dof_indices[i]) += cell_rhs(i);
    }

  hanging_nodes_only.condense (reference_matrix, reference_rhs);

  // use some other rhs vector as dummy for
  // application of Dirichlet conditions
  std::map<types::global_dof_index,double> boundary_values;
  VectorTools::interpolate_boundary_values (dof_handler,
                                            0,
                                            RightHandSide<dim>(),
                                            boundary_values);
  MatrixTools::apply_boundary_values (boundary_values,
                                      reference_matrix,
                                      test_rhs,
                                      reference_rhs);
}



template <int dim>
void AdvectionProblem<dim>::assemble_test_1 ()
{
  test_matrix = 0;
  test_rhs = 0;


  QGauss<dim> quadrature_formula (3);
  FEValues<dim> fe_values (fe, quadrature_formula,
                           update_values    |  update_gradients |
                           update_quadrature_points  |  update_JxW_values);

  const RightHandSide<dim> rhs_function;
  const unsigned int dofs_per_cell = fe.dofs_per_cell;
  const unsigned int n_q_points = quadrature_formula.size();

  FullMatrix<double>   cell_matrix (dofs_per_cell, dofs_per_cell);
  Vector<double>       cell_rhs (dofs_per_cell);

  std::vector<types::global_dof_index> local_dof_indices (dofs_per_cell);
  std::vector<double>  rhs_values (n_q_points);

  typename DoFHandler<dim>::active_cell_iterator
  cell = dof_handler.begin_active(),
  endc = dof_handler.end();
  for (; cell!=endc; ++cell)
    {
      cell_matrix = 0;
      cell_rhs = 0;
      fe_values.reinit (cell);

      rhs_function.value_list (fe_values.get_quadrature_points(),
                               rhs_values);

      Tensor<1,dim> advection_direction;
      advection_direction[0] = 1;
      advection_direction[1] = 1;
      advection_direction[dim-1] = -1;

      for (unsigned int q_point=0; q_point<n_q_points; ++q_point)
        for (unsigned int i=0; i<dofs_per_cell; ++i)
          {
            for (unsigned int j=0; j<dofs_per_cell; ++j)
              cell_matrix(i,j) += (fe_values.shape_value(i,q_point) *
                                   advection_direction *
                                   fe_values.shape_grad(j,q_point) *
                                   fe_values.JxW(q_point));

            cell_rhs(i) += (fe_values.shape_value(i,q_point) *
                            rhs_values[q_point] *
                            fe_values.JxW(q_point));
          }

      local_dof_indices.resize (dofs_per_cell);
      cell->get_dof_indices (local_dof_indices);

      test_matrix.add(local_dof_indices, cell_matrix);
      for (unsigned int i=0; i<dofs_per_cell; ++i)
        test_rhs(local_dof_indices[i]) += cell_rhs(i);

    }

  test_all_constraints.condense (test_matrix, test_rhs);

  test_equality();
}



template <int dim>
void AdvectionProblem<dim>::assemble_test_2 ()
{
  test_matrix = 0;
  test_rhs = 0;

  QGauss<dim> quadrature_formula (3);
  FEValues<dim> fe_values (fe, quadrature_formula,
                           update_values    |  update_gradients |
                           update_quadrature_points  |  update_JxW_values);

  const RightHandSide<dim> rhs_function;
  const unsigned int dofs_per_cell = fe.dofs_per_cell;
  const unsigned int n_q_points = quadrature_formula.size();

  FullMatrix<double>   cell_matrix (dofs_per_cell, dofs_per_cell);
  Vector<double>       cell_rhs (dofs_per_cell);

  std::vector<types::global_dof_index> local_dof_indices (dofs_per_cell);
  std::vector<double>  rhs_values (n_q_points);

  typename DoFHandler<dim>::active_cell_iterator
  cell = dof_handler.begin_active(),
  endc = dof_handler.end();
  for (; cell!=endc; ++cell)
    {
      cell_matrix = 0;
      cell_rhs = 0;
      fe_values.reinit (cell);

      rhs_function.value_list (fe_values.get_quadrature_points(),
                               rhs_values);

      Tensor<1,dim> advection_direction;
      advection_direction[0] = 1;
      advection_direction[1] = 1;
      advection_direction[dim-1] = -1;

      for (unsigned int q_point=0; q_point<n_q_points; ++q_point)
        for (unsigned int i=0; i<dofs_per_cell; ++i)
          {
            for (unsigned int j=0; j<dofs_per_cell; ++j)
              cell_matrix(i,j) += (fe_values.shape_value(i,q_point) *
                                   advection_direction *
                                   fe_values.shape_grad(j,q_point) *
                                   fe_values.JxW(q_point));

            cell_rhs(i) += (fe_values.shape_value(i,q_point) *
                            rhs_values[q_point] *
                            fe_values.JxW(q_point));
          }

      local_dof_indices.resize (dofs_per_cell);
      cell->get_dof_indices (local_dof_indices);

      test_all_constraints.distribute_local_to_global (cell_matrix,
                                                       cell_rhs,
                                                       local_dof_indices,
                                                       test_matrix,
                                                       test_rhs);
    }
  test_equality();
}


template <int dim>
void AdvectionProblem<dim>::run ()
{
  GridGenerator::hyper_ball (triangulation);
  triangulation.refine_global (3-dim);

  // manually refine the first two cells
  // and then one of these cells once again
  // to create some hanging nodes
  {
    typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active();
    cell->set_refine_flag();
  }
  triangulation.execute_coarsening_and_refinement();
  {
    // find the last cell and mark it
    // for refinement
    for (typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active();
         cell != dof_handler.end(); ++cell)
      if (++typename DoFHandler<dim>::active_cell_iterator(cell) ==
          dof_handler.end())
        cell->set_refine_flag();
  }
  triangulation.execute_coarsening_and_refinement();

  setup_system ();

  deallog << std::endl << std::endl
          << "  Number of active cells:       "
          << triangulation.n_active_cells()
          << std::endl
          << "  Number of degrees of freedom: "
          << dof_handler.n_dofs()
          << std::endl
          << "  Number of constraints       : "
          << hanging_nodes_only.n_constraints()
          << std::endl;

  assemble_reference ();
  assemble_test_1 ();
  assemble_test_2 ();
}



int main ()
{
  deallog << std::setprecision (2);
  logfile << std::setprecision (2);
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-12);

  {
    AdvectionProblem<2> advection_problem;
    advection_problem.run ();
  }
  {
    AdvectionProblem<3> advection_problem;
    advection_problem.run ();
  }
}
