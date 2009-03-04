//----------------------------  step-11.cc  ---------------------------
//    $Id$
//    Version: $Name$ 
//
//    Copyright (C) 2005, 2006 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  step-11.cc  ---------------------------


// a un-hp-ified version of hp/step-11


#include "../tests.h"
#include <base/logstream.h>
#include <fstream>
std::ofstream logfile("step-11/output");

#include <base/quadrature_lib.h>
#include <base/function.h>
#include "../tests.h"
#include <base/logstream.h>
#include <base/table_handler.h>
#include <lac/vector.h>
#include <lac/sparse_matrix.h>
#include <lac/solver_cg.h>
#include <lac/precondition.h>
#include <grid/tria.h>
#include <grid/grid_generator.h>
#include <grid/tria_boundary_lib.h>
#include <grid/tria_accessor.h>
#include <grid/tria_iterator.h>
#include <dofs/dof_handler.h>
#include <lac/constraint_matrix.h>
#include <dofs/dof_accessor.h>
#include <dofs/dof_tools.h>
#include <fe/fe_q.h>
#include <fe/fe_values.h>
#include <fe/mapping_q.h>
#include <numerics/vectors.h>
#include <numerics/matrices.h>

#include <lac/compressed_sparsity_pattern.h>

#include <algorithm>
#include <iomanip>
#include <iomanip>
#include <cmath>


template <int dim>
class LaplaceProblem 
{
  public:
    LaplaceProblem (const unsigned int mapping_degree);
    void run ();
    
  private:
    void setup_system ();
    void assemble_and_solve ();
    void solve ();

    Triangulation<dim>   triangulation;
    FE_Q<dim>            fe;
    DoFHandler<dim>      dof_handler;
    MappingQ<dim>        mapping;

    SparsityPattern      sparsity_pattern;
    SparseMatrix<double> system_matrix;
    ConstraintMatrix     mean_value_constraints;

    Vector<double>       solution;
    Vector<double>       system_rhs;

    TableHandler         output_table;
};



template <int dim>
LaplaceProblem<dim>::LaplaceProblem (const unsigned int mapping_degree) :
                fe (1),
		dof_handler (triangulation),
		mapping (mapping_degree)
{
  deallog << "Using mapping with degree " << mapping_degree << ":"
	    << std::endl
	    << "============================"
	    << std::endl;
}



template <int dim>
void LaplaceProblem<dim>::setup_system ()
{
  dof_handler.distribute_dofs (fe);
  solution.reinit (dof_handler.n_dofs());
  system_rhs.reinit (dof_handler.n_dofs());

  std::vector<bool> boundary_dofs (dof_handler.n_dofs(), false);
  DoFTools::extract_boundary_dofs (dof_handler, std::vector<bool>(1,true),
				   boundary_dofs);

  const unsigned int first_boundary_dof
    = std::distance (boundary_dofs.begin(),
		     std::find (boundary_dofs.begin(),
				boundary_dofs.end(),
				true));

  mean_value_constraints.clear ();
  mean_value_constraints.add_line (first_boundary_dof);
  for (unsigned int i=first_boundary_dof+1; i<dof_handler.n_dofs(); ++i)
    if (boundary_dofs[i] == true)
      mean_value_constraints.add_entry (first_boundary_dof,
					i, -1);
  mean_value_constraints.close ();

  CompressedSparsityPattern csp (dof_handler.n_dofs(),
				 dof_handler.n_dofs());
  DoFTools::make_sparsity_pattern (dof_handler, csp);
  mean_value_constraints.condense (csp);

  sparsity_pattern.copy_from (csp);
  system_matrix.reinit (sparsity_pattern);
}



template <int dim>
void LaplaceProblem<dim>::assemble_and_solve () 
{

  const unsigned int gauss_degree
    = std::max (static_cast<unsigned int>(std::ceil(1.*(mapping.get_degree()+1)/2)),
		2U);
  MatrixTools::create_laplace_matrix (mapping, dof_handler,
				      QGauss<dim>(gauss_degree),
				      system_matrix);
  VectorTools::create_right_hand_side (mapping, dof_handler,
				       QGauss<dim>(gauss_degree),
				       ConstantFunction<dim>(-2),
				       system_rhs);
  Vector<double> tmp (system_rhs.size());
  VectorTools::create_boundary_right_hand_side (mapping, dof_handler,
						QGauss<dim-1>(gauss_degree),
						ConstantFunction<dim>(1),
						tmp);
  system_rhs += tmp;

  mean_value_constraints.condense (system_matrix);
  mean_value_constraints.condense (system_rhs);  

  solve ();
  mean_value_constraints.distribute (solution);

  Vector<float> norm_per_cell (triangulation.n_active_cells());
  VectorTools::integrate_difference (mapping, dof_handler,
				     solution,
				     ZeroFunction<dim>(),
				     norm_per_cell,
				     QGauss<dim>(gauss_degree+1),
				     VectorTools::H1_seminorm);
  const double norm = norm_per_cell.l2_norm();

  output_table.add_value ("cells", triangulation.n_active_cells());
  output_table.add_value ("|u|_1", norm);
  output_table.add_value ("error", std::fabs(norm-std::sqrt(3.14159265358/2)));
}



template <int dim>
void LaplaceProblem<dim>::solve () 
{
  SolverControl           solver_control (1000, 1e-12);
  SolverCG<>              cg (solver_control);

  PreconditionSSOR<> preconditioner;
  preconditioner.initialize(system_matrix, 1.2);

  cg.solve (system_matrix, solution, system_rhs,
	    preconditioner);
}



template <int dim>
void LaplaceProblem<dim>::run () 
{
  GridGenerator::hyper_ball (triangulation);
  static const HyperBallBoundary<dim> boundary;
  triangulation.set_boundary (0, boundary);
  
  for (unsigned int cycle=0; cycle<6; ++cycle, triangulation.refine_global(1))
    {
      setup_system ();
      assemble_and_solve ();
    };

  output_table.set_precision("|u|_1", 6);
  output_table.set_precision("error", 6);
  output_table.write_text (deallog.get_file_stream());
  deallog << std::endl;
}

    

int main () 
{
  try
    {
      deallog << std::setprecision(2);
  
      deallog.attach(logfile);
      deallog.depth_console(0);
      deallog.threshold_double(1.e-10);  

      for (unsigned int mapping_degree=1; mapping_degree<=3; ++mapping_degree)
	LaplaceProblem<2>(mapping_degree).run ();
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

  return 0;
}
