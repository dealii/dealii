// ---------------------------------------------------------------------
//
// Copyright (C) 2015 by the deal.II authors
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



// solves a 2D Poisson equation for FE_Q elements using FE_Q_iso_Q1 elements
// for the MueLu preconditioner. The problem is taken from step-4

#include "../tests.h"
#include <deal.II/lac/trilinos_sparse_matrix.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_q_iso_q1.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/constraint_matrix.h>
#include <deal.II/lac/compressed_sparsity_pattern.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/lac/trilinos_precondition.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/base/function.h>
#include <deal.II/grid/tria.h>

#include <fstream>
#include <iomanip>


template <int dim>
class Step4
{
public:
  Step4 ();
  void run ();

private:
  void make_grid ();
  void setup_system();
  void assemble_system ();
  void assemble_preconditioner ();
  void solve ();

  Triangulation<dim>   triangulation;
  FE_Q<dim>            fe;
  DoFHandler<dim>      dof_handler;

  FE_Q_iso_Q1<dim>     fe_precondition;
  DoFHandler<dim>      dof_handler_precondition;

  ConstraintMatrix     constraints;

  TrilinosWrappers::SparseMatrix system_matrix;
  TrilinosWrappers::SparseMatrix preconditioner_matrix;

  Vector<double>       solution;
  Vector<double>       system_rhs;
};


template <int dim>
class RightHandSide : public Function<dim>
{
public:
  RightHandSide () : Function<dim>() {}

  virtual double value (const Point<dim>   &p,
                        const unsigned int  component = 0) const;
};



template <int dim>
class BoundaryValues : public Function<dim>
{
public:
  BoundaryValues () : Function<dim>() {}

  virtual double value (const Point<dim>   &p,
                        const unsigned int  component = 0) const;
};




template <int dim>
double RightHandSide<dim>::value (const Point<dim> &p,
                                  const unsigned int /*component*/) const
{
  double return_value = 0;
  for (unsigned int i=0; i<dim; ++i)
    return_value += 4*std::pow(p(i), 4);

  return return_value;
}



template <int dim>
double BoundaryValues<dim>::value (const Point<dim> &p,
                                   const unsigned int /*component*/) const
{
  return p.square();
}



template <int dim>
Step4<dim>::Step4 ()
  :
  fe (3),
  dof_handler (triangulation),
  fe_precondition (3),
  dof_handler_precondition(triangulation)
{}


template <int dim>
void Step4<dim>::make_grid ()
{
  GridGenerator::hyper_cube (triangulation, -1, 1);
  triangulation.refine_global (4);
}



template <int dim>
void Step4<dim>::setup_system ()
{
  dof_handler.distribute_dofs (fe);
  dof_handler_precondition.distribute_dofs (fe_precondition);

  constraints.clear();
  std::map<unsigned int,double> boundary_values;
  VectorTools::interpolate_boundary_values (dof_handler,
                                            0,
                                            BoundaryValues<dim>(),
                                            constraints);
  constraints.close();

  CompressedSparsityPattern c_sparsity(dof_handler.n_dofs());
  DoFTools::make_sparsity_pattern (dof_handler, c_sparsity, constraints, false);
  system_matrix.reinit (c_sparsity);
  preconditioner_matrix.reinit(c_sparsity);

  solution.reinit (dof_handler.n_dofs());
  system_rhs.reinit (dof_handler.n_dofs());
}


template <int dim>
void Step4<dim>::assemble_system ()
{
  QGauss<dim>  quadrature_formula(fe.degree+1);

  const RightHandSide<dim> right_hand_side;

  FEValues<dim> fe_values (fe, quadrature_formula,
                           update_values   | update_gradients |
                           update_quadrature_points | update_JxW_values);

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

      for (unsigned int q_point=0; q_point<n_q_points; ++q_point)
        for (unsigned int i=0; i<dofs_per_cell; ++i)
          {
            for (unsigned int j=0; j<dofs_per_cell; ++j)
              cell_matrix(i,j) += (fe_values.shape_grad (i, q_point) *
                                   fe_values.shape_grad (j, q_point) *
                                   fe_values.JxW (q_point));

            cell_rhs(i) += (fe_values.shape_value (i, q_point) *
                            right_hand_side.value (fe_values.quadrature_point (q_point)) *
                            fe_values.JxW (q_point));
          }

      cell->get_dof_indices (local_dof_indices);
      constraints.distribute_local_to_global(cell_matrix, cell_rhs,
                                             local_dof_indices,
                                             system_matrix, system_rhs);
    }
  system_matrix.compress(VectorOperation::add);
}




template <int dim>
void Step4<dim>::assemble_preconditioner ()
{
  QIterated<dim>  quadrature_formula(QGauss<1>(2), fe.degree);

  const RightHandSide<dim> right_hand_side;

  FEValues<dim> fe_values (fe_precondition, quadrature_formula,
                           update_values   | update_gradients |
                           update_quadrature_points | update_JxW_values);

  const unsigned int   dofs_per_cell = fe_precondition.dofs_per_cell;
  const unsigned int   n_q_points    = quadrature_formula.size();

  FullMatrix<double>   cell_matrix (dofs_per_cell, dofs_per_cell);
  Vector<double>       cell_rhs (dofs_per_cell);

  std::vector<types::global_dof_index> local_dof_indices (dofs_per_cell);

  typename DoFHandler<dim>::active_cell_iterator
  cell = dof_handler_precondition.begin_active(),
  endc = dof_handler_precondition.end();

  for (; cell!=endc; ++cell)
    {
      fe_values.reinit (cell);
      cell_matrix = 0;
      cell_rhs = 0;

      for (unsigned int q_point=0; q_point<n_q_points; ++q_point)
        for (unsigned int i=0; i<dofs_per_cell; ++i)
          {
            for (unsigned int j=0; j<dofs_per_cell; ++j)
              cell_matrix(i,j) += (fe_values.shape_grad (i, q_point) *
                                   fe_values.shape_grad (j, q_point) *
                                   fe_values.JxW (q_point));
          }

      cell->get_dof_indices (local_dof_indices);
      constraints.distribute_local_to_global(cell_matrix,
                                             local_dof_indices,
                                             preconditioner_matrix);
    }
  preconditioner_matrix.compress(VectorOperation::add);
}


template <int dim>
void Step4<dim>::solve ()
{

  // variant 1: solve with MueLu
  deallog.push(Utilities::int_to_string(dof_handler.n_dofs(),5));
  deallog.push("MueLu_Q");
  {
    solution = 0;
    SolverControl           solver_control (1000, 1e-8);
    SolverCG<>              solver (solver_control);
    TrilinosWrappers::PreconditionAMGMueLu preconditioner;
    preconditioner.initialize(system_matrix);
    check_solver_within_range(
      solver.solve (system_matrix, solution, system_rhs,
                    preconditioner),
      solver_control.last_step(), 25, 34);
  }
  deallog.pop();

  deallog.push("MueLu_Q_iso_Q1");
  {
    solution = 0;
    SolverControl           solver_control (1000, 1e-8);
    SolverCG<>              solver (solver_control);
    TrilinosWrappers::PreconditionAMGMueLu preconditioner;
    preconditioner.initialize(preconditioner_matrix);
    check_solver_within_range(
      solver.solve (system_matrix, solution, system_rhs,
                    preconditioner),
      solver_control.last_step(), 24, 40);
  }
  deallog.pop();
  deallog << std::endl;
  deallog.pop();
}



template <int dim>
void Step4<dim>::run()
{
  for (unsigned int cycle = 0; cycle < 2; ++cycle)
    {
      if (cycle == 0)
        make_grid();
      else
        triangulation.refine_global(1);

      setup_system();
      assemble_system();
      assemble_preconditioner();
      solve();
    }
}


int main (int argc, char **argv)
{
  std::ofstream logfile("output");
  deallog.attach(logfile);
  deallog.threshold_double(1.e-10);

  Utilities::MPI::MPI_InitFinalize mpi_initialization (argc, argv, testing_max_num_threads());

  try
    {
      Step4<2> test;
      test.run();
    }
  catch (std::exception &exc)
    {
      deallog << std::endl << std::endl
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
      deallog << std::endl << std::endl
              << "----------------------------------------------------"
              << std::endl;
      deallog << "Unknown exception!" << std::endl
              << "Aborting!" << std::endl
              << "----------------------------------------------------"
              << std::endl;
      return 1;
    };
}
