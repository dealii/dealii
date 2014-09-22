/* ---------------------------------------------------------------------
 *
 * Copyright (C) 2013 by the deal.II authors
 *
 * This file is part of the deal.II library.
 *
 * The deal.II library is free software; you can use it, redistribute
 * it, and/or modify it under the terms of the GNU Lesser General
 * Public License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 * The full text of the license can be found in the file LICENSE at
 * the top level of the deal.II distribution.
 *
 * ---------------------------------------------------------------------
 */

#include <deal.II/base/logstream.h>
#include <deal.II/base/table_handler.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function.h>
#include <deal.II/base/utilities.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/petsc_sparse_matrix.h>
#include <deal.II/lac/petsc_vector.h>
#include <deal.II/lac/slepc_solver.h>

#include <fstream>
#include <iostream>

using namespace dealii;

// Test that deal.II is working with SLEPc by solving the Laplace's
// eigenspectrum problem in 2d.
class LaplaceEigenspectrumProblem
{
public:
  LaplaceEigenspectrumProblem ();
  void run ();
  
private:
  void setup_system ();
  void assemble_system ();
  void solve ();
  
  Triangulation<2> triangulation;
  FE_Q<2>          fe;
  DoFHandler<2>    dof_handler;
  
  PETScWrappers::SparseMatrix        A, B;
  std::vector<PETScWrappers::Vector> x;
  std::vector<double>                lambda;  
  ConstraintMatrix                   constraints;

  TableHandler output_table;
};

LaplaceEigenspectrumProblem::LaplaceEigenspectrumProblem ()
  :
  fe (1),
  dof_handler (triangulation)
{}

void LaplaceEigenspectrumProblem::setup_system ()
{
  dof_handler.distribute_dofs (fe);
  
  constraints.clear ();
  DoFTools::make_zero_boundary_constraints (dof_handler, constraints);
  constraints.close ();
  
  A.reinit (dof_handler.n_dofs(), dof_handler.n_dofs(),
	    dof_handler.max_couplings_between_dofs());
  B.reinit (dof_handler.n_dofs(), dof_handler.n_dofs(),
	    dof_handler.max_couplings_between_dofs());
  
  x.resize (1);
  x[0].reinit (dof_handler.n_dofs ());
  lambda.resize (1);
  lambda[0] = 0.;

  // some output
  output_table.add_value ("cells", triangulation.n_active_cells());
  output_table.add_value ("dofs",  dof_handler.n_dofs());
}

void LaplaceEigenspectrumProblem::assemble_system ()
{
  QGauss<2> quadrature_formula(2);
  
  FEValues<2> fe_values (fe, quadrature_formula,
			 update_values            | 
			 update_gradients         |
			 update_quadrature_points | 
			 update_JxW_values);
  
  const unsigned int dofs_per_cell = fe.dofs_per_cell;
  const unsigned int n_q_points    = quadrature_formula.size();
  
  FullMatrix<double> cell_A (dofs_per_cell, dofs_per_cell);
  FullMatrix<double> cell_B (dofs_per_cell, dofs_per_cell);
  
  std::vector<types::global_dof_index> local_dof_indices (dofs_per_cell);
  
  typename DoFHandler<2>::active_cell_iterator
    cell = dof_handler.begin_active (),
    endc = dof_handler.end ();
  
  for (; cell!=endc; ++cell)
    {
      fe_values.reinit (cell);
      cell_A = 0;
      cell_B = 0;
      
      for (unsigned int q_point=0; q_point<n_q_points; ++q_point)
	for (unsigned int i=0; i<dofs_per_cell; ++i)
	  for (unsigned int j=0; j<dofs_per_cell; ++j)
	    {
	      cell_A (i, j)
		+= 
		fe_values.shape_grad (i, q_point) *
		fe_values.shape_grad (j, q_point)
		* 
		fe_values.JxW (q_point);
	      
	      cell_B (i, j)
                += 
		fe_values.shape_value (i, q_point) *
		fe_values.shape_value (j, q_point)
		* 
		fe_values.JxW (q_point);
	    }
      
      cell->get_dof_indices (local_dof_indices);
      
      constraints.distribute_local_to_global (cell_A, local_dof_indices, A);
      constraints.distribute_local_to_global (cell_B, local_dof_indices, B);
    }
  
  A.compress (VectorOperation::add);
  B.compress (VectorOperation::add);
}

void LaplaceEigenspectrumProblem::solve ()
{
  SolverControl solver_control (1000, 1e-03);
  SLEPcWrappers::SolverArnoldi eigensolver (solver_control);
  eigensolver.set_which_eigenpairs (EPS_SMALLEST_REAL);
  eigensolver.solve (A, B, lambda, x, x.size());

  // some output
  output_table.add_value ("lambda", lambda[0]);
  output_table.add_value ("error", std::fabs(2.-lambda[0]));
}

void LaplaceEigenspectrumProblem::run ()
{
  const double radius = dealii::numbers::PI/2.;
  GridGenerator::hyper_cube (triangulation, -radius, radius);

  // set the old eigenvalue to a silly number.
  double old_lambda = 1000;

  for (unsigned int c=0; c<5; ++c)
    {
      // obtain numerical result
      triangulation.refine_global (1);
      setup_system ();
      assemble_system ();
      solve ();

      // check energy convergence with previous result
      AssertThrow (lambda[0]<old_lambda, ExcMessage("solution is not converging"));
      old_lambda = lambda[0];
    }

  // push back analytic result
  output_table.add_value ("cells", "inf");
  output_table.add_value ("dofs",  "inf");
  output_table.add_value ("lambda", 2.);
  output_table.add_value ("error", "-");

  // finialise output
  output_table.write_text (std::cout);
  deallog << std::endl;
}


int main (int argc, char **argv)
{
  try
    {
      dealii::Utilities::MPI::MPI_InitFinalize mpi_initialization (argc, argv, 1);
      {
        deallog.depth_console (0);
        LaplaceEigenspectrumProblem problem;
        problem.run ();
      }
    }

  catch (std::exception &exc)
    {
      std::cerr << std::endl << std::endl
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
      std::cerr << std::endl << std::endl
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
