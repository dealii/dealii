// ---------------------------------------------------------------------
//
// Copyright (C) 2005 - 2013 by the deal.II authors
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



// a test adapted from little programs by Markus Buerg that makes sure
// we compute the same solution for hp and non-hp cases


#include "../tests.h"
#include <deal.II/base/logstream.h>
#include <deal.II/grid/tria.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/grid/grid_generator.h>

#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/dofs/dof_accessor.h>

#include <deal.II/fe/fe_q.h>

#include <deal.II/dofs/dof_tools.h>

#include <deal.II/hp/fe_values.h>
#include <deal.II/base/quadrature_lib.h>

#include <deal.II/base/function.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/matrix_tools.h>

#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/precondition.h>

#include <deal.II/numerics/data_out.h>
#include <fstream>
#include <iostream>

std::ofstream logfile("output");


template <int dim>
class ExactSolution: public Function<dim>
{
public:
  ExactSolution ();
  virtual double value (const Point<dim> &p, const unsigned int component = 0) const;
};

template <int dim>
double ExactSolution<dim>::value (const Point<dim> &p, const unsigned int) const
{
  return p (0) * p (0);
}

template <int dim>
ExactSolution<dim>::ExactSolution (): Function<dim> ()
{
}



namespace with_hp
{
  template <int dim>
  class LaplaceProblem
  {
  public:
    LaplaceProblem<dim> ();

    void run (Vector<double> &sol);

  private:
    const ExactSolution<dim> exact_solution;
    void make_grid_and_dofs ();
    void assemble_system ();
    void solve ();

    Triangulation<dim>     triangulation;
    hp::FECollection<dim>      fe;
    hp::DoFHandler<dim>        dof_handler;
    hp::QCollection<dim>   quadrature;
    SparsityPattern      sparsity_pattern;
    SparseMatrix<double> system_matrix;

    Vector<double>       solution;
    Vector<double>       system_rhs;
  };

  template <int dim>
  LaplaceProblem<dim>::LaplaceProblem () :
    dof_handler (triangulation)
  {
    fe.push_back (FE_Q<dim> (1));
    quadrature.push_back (QGauss<dim> (2));
  }

  template <int dim>
  void LaplaceProblem<dim>::make_grid_and_dofs ()
  {
    GridGenerator::hyper_cube (triangulation);
    triangulation.refine_global (2);
    deallog << "Number of active cells: "
            << triangulation.n_active_cells()
            << std::endl;
    deallog << "Total number of cells: "
            << triangulation.n_cells()
            << std::endl;

    dof_handler.distribute_dofs (fe);
    deallog << "Number of degrees of freedom: "
            << dof_handler.n_dofs()
            << std::endl;

    sparsity_pattern.reinit (dof_handler.n_dofs(),
                             dof_handler.n_dofs(),
                             dof_handler.max_couplings_between_dofs());
    DoFTools::make_sparsity_pattern (dof_handler, sparsity_pattern);
    sparsity_pattern.compress();

    system_matrix.reinit (sparsity_pattern);

    solution.reinit (dof_handler.n_dofs());
    system_rhs.reinit (dof_handler.n_dofs());
  }

  template <int dim>
  void LaplaceProblem<dim>::assemble_system ()
  {
    hp::FEValues<dim> fe_values (fe, quadrature,
                                 update_values | update_gradients | update_JxW_values);

    unsigned int   dofs_per_cell;
    unsigned int   n_q_points;

    FullMatrix<double>   cell_matrix;
    Vector<double>       cell_rhs;

    std::vector<types::global_dof_index> local_dof_indices;

    typename hp::DoFHandler<dim>::active_cell_iterator
    cell = dof_handler.begin_active(),
    endc = dof_handler.end();
    for (; cell!=endc; ++cell)
      {
        fe_values.reinit (cell);
        dofs_per_cell = cell->get_fe ().dofs_per_cell;
        cell_matrix.reinit (dofs_per_cell, dofs_per_cell);
        cell_matrix = 0;
        cell_rhs.reinit (dofs_per_cell);
        cell_rhs = 0;
        n_q_points = fe_values.get_present_fe_values ().n_quadrature_points;

        for (unsigned int i=0; i<dofs_per_cell; ++i)
          for (unsigned int j=0; j<dofs_per_cell; ++j)
            for (unsigned int q_point=0; q_point<n_q_points; ++q_point)
              cell_matrix(i,j) += (fe_values.get_present_fe_values ().shape_grad (i, q_point) *
                                   fe_values.get_present_fe_values ().shape_grad (j, q_point) *
                                   fe_values.get_present_fe_values ().JxW (q_point));

        for (unsigned int i=0; i<dofs_per_cell; ++i)
          for (unsigned int q_point=0; q_point<n_q_points; ++q_point)
            cell_rhs(i) += (-2.0 * fe_values.get_present_fe_values ().shape_value (i, q_point) *
                            fe_values.get_present_fe_values ().JxW (q_point));

        local_dof_indices.resize (dofs_per_cell);
        cell->get_dof_indices (local_dof_indices);

        for (unsigned int i=0; i<dofs_per_cell; ++i)
          for (unsigned int j=0; j<dofs_per_cell; ++j)
            system_matrix.add (local_dof_indices[i],
                               local_dof_indices[j],
                               cell_matrix(i,j));

        for (unsigned int i=0; i<dofs_per_cell; ++i)
          system_rhs(local_dof_indices[i]) += cell_rhs(i);
      }


    std::map<types::global_dof_index,double> boundary_values;
    VectorTools::interpolate_boundary_values (dof_handler,
                                              0,
                                              exact_solution,
                                              boundary_values);
    VectorTools::interpolate_boundary_values (dof_handler,
                                              1,
                                              exact_solution,
                                              boundary_values);
    MatrixTools::apply_boundary_values (boundary_values,
                                        system_matrix,
                                        solution,
                                        system_rhs);
  }

  template <int dim>
  void LaplaceProblem<dim>::solve ()
  {
    SolverControl           solver_control (1000, 1e-12);
    SolverCG<>              cg (solver_control);

    cg.solve (system_matrix, solution, system_rhs,
              PreconditionIdentity());
  }

  template <int dim>
  void LaplaceProblem<dim>::run (Vector<double> &sol)
  {
    make_grid_and_dofs ();
    assemble_system ();
    solve ();

    sol = solution;
  }
}


namespace without_hp
{
  template <int dim>
  class LaplaceProblem
  {
  public:
    LaplaceProblem<dim> ();

    void run (Vector<double> &sol);

  private:
    const ExactSolution<dim> exact_solution;
    void make_grid_and_dofs ();
    void assemble_system ();
    void solve ();

    Triangulation<dim>     triangulation;
    FE_Q<dim>              fe;
    DoFHandler<dim>        dof_handler;

    SparsityPattern      sparsity_pattern;
    SparseMatrix<double> system_matrix;

    Vector<double>       solution;
    Vector<double>       system_rhs;
  };

  template <int dim>
  LaplaceProblem<dim>::LaplaceProblem () :
    fe (1),
    dof_handler (triangulation)
  {}

  template <int dim>
  void LaplaceProblem<dim>::make_grid_and_dofs ()
  {
    GridGenerator::hyper_cube (triangulation);
    triangulation.refine_global (2);
    deallog << "Number of active cells: "
            << triangulation.n_active_cells()
            << std::endl;
    deallog << "Total number of cells: "
            << triangulation.n_cells()
            << std::endl;

    dof_handler.distribute_dofs (fe);
    deallog << "Number of degrees of freedom: "
            << dof_handler.n_dofs()
            << std::endl;

    sparsity_pattern.reinit (dof_handler.n_dofs(),
                             dof_handler.n_dofs(),
                             dof_handler.max_couplings_between_dofs());
    DoFTools::make_sparsity_pattern (dof_handler, sparsity_pattern);
    sparsity_pattern.compress();

    system_matrix.reinit (sparsity_pattern);

    solution.reinit (dof_handler.n_dofs());
    system_rhs.reinit (dof_handler.n_dofs());
  }

  template <int dim>
  void LaplaceProblem<dim>::assemble_system ()
  {
    QGauss<dim>  quadrature_formula(2);
    FEValues<dim> fe_values (fe, quadrature_formula,
                             update_values | update_gradients | update_JxW_values);

    const unsigned int   dofs_per_cell = fe.dofs_per_cell;
    const unsigned int   n_q_points    = quadrature_formula.size();

    FullMatrix<double>   cell_matrix (dofs_per_cell, dofs_per_cell);
    Vector<double>       cell_rhs (dofs_per_cell);

    std::vector<types::global_dof_index> local_dof_indices (dofs_per_cell);

    typename DoFHandler<dim>::active_cell_iterator
    cell = dof_handler.begin_active(),
    endc = dof_handler.end();
    for (; cell!=endc; ++cell)
      {
        fe_values.reinit (cell);

        cell_matrix = 0;
        cell_rhs = 0;

        for (unsigned int i=0; i<dofs_per_cell; ++i)
          for (unsigned int j=0; j<dofs_per_cell; ++j)
            for (unsigned int q_point=0; q_point<n_q_points; ++q_point)
              cell_matrix(i,j) += (fe_values.shape_grad (i, q_point) *
                                   fe_values.shape_grad (j, q_point) *
                                   fe_values.JxW (q_point));

        for (unsigned int i=0; i<dofs_per_cell; ++i)
          for (unsigned int q_point=0; q_point<n_q_points; ++q_point)
            cell_rhs(i) += (-2.0 * fe_values.shape_value (i, q_point) *
                            fe_values.JxW (q_point));

        cell->get_dof_indices (local_dof_indices);

        for (unsigned int i=0; i<dofs_per_cell; ++i)
          for (unsigned int j=0; j<dofs_per_cell; ++j)
            system_matrix.add (local_dof_indices[i],
                               local_dof_indices[j],
                               cell_matrix(i,j));

        for (unsigned int i=0; i<dofs_per_cell; ++i)
          system_rhs(local_dof_indices[i]) += cell_rhs(i);
      }


    std::map<types::global_dof_index,double> boundary_values;
    VectorTools::interpolate_boundary_values (dof_handler,
                                              0,
                                              exact_solution,
                                              boundary_values);
    VectorTools::interpolate_boundary_values (dof_handler,
                                              1,
                                              exact_solution,
                                              boundary_values);
    MatrixTools::apply_boundary_values (boundary_values,
                                        system_matrix,
                                        solution,
                                        system_rhs);
  }

  template <int dim>
  void LaplaceProblem<dim>::solve ()
  {
    SolverControl           solver_control (1000, 1e-12);
    SolverCG<>              cg (solver_control);

    cg.solve (system_matrix, solution, system_rhs,
              PreconditionIdentity());
  }

  template <int dim>
  void LaplaceProblem<dim>::run (Vector<double> &sol)
  {
    make_grid_and_dofs ();
    assemble_system ();
    solve ();

    sol = solution;
  }
}


template <int dim>
void test ()
{
  deallog << "dim=" << dim << std::endl;

  Vector<double> sol1, sol2;
  with_hp::LaplaceProblem<dim>().run (sol1);
  without_hp::LaplaceProblem<dim>().run (sol2);

  deallog << sol1.l2_norm() << ' ' << sol2.l2_norm() << std::endl;

  sol1 -= sol2;
  Assert (sol1.l2_norm () <= 1e-8 * sol2.l2_norm(),
          ExcInternalError());
}



int main ()
{
  logfile.precision(2);
  deallog << std::setprecision(2);

  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  test<1> ();
  test<2> ();
  test<3> ();

  return 0;
}
