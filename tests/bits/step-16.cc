//----------------------------  step-16.cc  ---------------------------
//    $Id$
//    Version: $Name$ 
//
//    Copyright (C) 2005, 2006, 2008, 2010 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  step-16.cc  ---------------------------


// a un-hp-ified version of hp/step-16


#include "../tests.h"
#include <base/logstream.h>
#include <fstream>
std::ofstream logfile("step-16/output");

#include <base/quadrature_lib.h>
#include <base/function.h>
#include "../tests.h"
#include <base/logstream.h>
#include <lac/vector.h>
#include <lac/full_matrix.h>
#include <lac/sparse_matrix.h>
#include <lac/solver_cg.h>
#include <lac/precondition.h>
#include <grid/tria.h>
#include <grid/tria_accessor.h>
#include <grid/tria_iterator.h>
#include <grid/grid_generator.h>
#include <dofs/dof_handler.h>
#include <dofs/dof_accessor.h>
#include <dofs/dof_tools.h>
#include <fe/fe_q.h>
#include <fe/fe_values.h>
#include <numerics/vectors.h>
#include <numerics/matrices.h>
#include <numerics/data_out.h>
#include <multigrid/multigrid.h>
#include <multigrid/mg_dof_handler.h>
#include <multigrid/mg_dof_accessor.h>

#include <multigrid/mg_transfer.h>
#include <multigrid/mg_tools.h>

#include <multigrid/mg_coarse.h>
#include <multigrid/mg_smoother.h>
#include <multigrid/mg_matrix.h>

#include <fstream>
#include <sstream>


template <int dim>
class LaplaceProblem 
{
  public:
    LaplaceProblem ();
    void run ();
    
  private:
    void setup_system ();
    void assemble_system ();
    void assemble_multigrid ();
    void solve ();
    void output_results (const unsigned int cycle) const;

    Triangulation<dim>   triangulation;
    FE_Q<dim>            fe;
    MGDoFHandler<dim>      mg_dof_handler;

    SparsityPattern      sparsity_pattern;
    SparseMatrix<double> system_matrix;

    MGLevelObject<SparsityPattern> mg_sparsity;
    MGLevelObject<SparseMatrix<float> > mg_matrices;
    
    Vector<double>       solution;
    Vector<double>       system_rhs;
};


template <int dim>
LaplaceProblem<dim>::LaplaceProblem () :
                fe (1),
		mg_dof_handler (triangulation)
{}



template <int dim>
void LaplaceProblem<dim>::setup_system ()
{
  mg_dof_handler.distribute_dofs (fe);

  deallog << "   Number of degrees of freedom: "
	    << mg_dof_handler.n_dofs()
	    << std::endl;

  sparsity_pattern.reinit (mg_dof_handler.n_dofs(),
			   mg_dof_handler.n_dofs(),
			   mg_dof_handler.max_couplings_between_dofs());
  DoFTools::make_sparsity_pattern (static_cast<const DoFHandler<dim>&>(mg_dof_handler),
                                   sparsity_pattern);
  sparsity_pattern.compress();

  system_matrix.reinit (sparsity_pattern);

  solution.reinit (mg_dof_handler.n_dofs());
  system_rhs.reinit (mg_dof_handler.n_dofs());

  const unsigned int nlevels = triangulation.n_levels();
  mg_matrices.resize(0, nlevels-1);
  mg_sparsity.resize(0, nlevels-1);

  for (unsigned int level=0;level<nlevels;++level)
    {
      mg_sparsity[level].reinit (mg_dof_handler.n_dofs(level),
				 mg_dof_handler.n_dofs(level),
				 mg_dof_handler.max_couplings_between_dofs());
      MGTools::make_sparsity_pattern (mg_dof_handler, mg_sparsity[level], level);
      mg_sparsity[level].compress();
      mg_matrices[level].reinit(mg_sparsity[level]);
    }
}

template <int dim>
void LaplaceProblem<dim>::assemble_system () 
{  
  QGauss<dim>  quadrature_formula(2);

  FEValues<dim> fe_values (fe, quadrature_formula, 
			   update_values   | update_gradients |
                           update_q_points | update_JxW_values);

  const unsigned int   dofs_per_cell = fe.dofs_per_cell;
  const unsigned int   n_q_points    = quadrature_formula.size();

  FullMatrix<double>   cell_matrix (dofs_per_cell, dofs_per_cell);
  Vector<double>       cell_rhs (dofs_per_cell);

  std::vector<unsigned int> local_dof_indices (dofs_per_cell);

  typename DoFHandler<dim>::active_cell_iterator cell = mg_dof_handler.begin_active(),
						 endc = mg_dof_handler.end();
  for (; cell!=endc; ++cell)
    {
      cell_matrix = 0;
      cell_rhs = 0;

      fe_values.reinit (cell);
      
      for (unsigned int q_point=0; q_point<n_q_points; ++q_point)
	for (unsigned int i=0; i<dofs_per_cell; ++i)
	  {
	    for (unsigned int j=0; j<dofs_per_cell; ++j)
	      cell_matrix(i,j) += ((fe_values.shape_grad(i,q_point)
				    * fe_values.shape_grad(j,q_point))
				    + fe_values.shape_value(i,q_point)
				    * fe_values.shape_value(j,q_point))
				   * fe_values.JxW(q_point);

	    cell_rhs(i) += (fe_values.shape_value(i,q_point)
			    * 1.0 * fe_values.JxW(q_point));
	  };


      cell->get_dof_indices (local_dof_indices);
      for (unsigned int i=0; i<dofs_per_cell; ++i)
	{
	  for (unsigned int j=0; j<dofs_per_cell; ++j)
	    system_matrix.add (local_dof_indices[i],
			       local_dof_indices[j],
			       cell_matrix(i,j));
	  
	  system_rhs(local_dof_indices[i]) += cell_rhs(i);
	};
    };

}



template <int dim>
void LaplaceProblem<dim>::assemble_multigrid () 
{  
  QGauss<dim>  quadrature_formula(2);

  FEValues<dim> fe_values (fe, quadrature_formula, 
			   update_values   | update_gradients |
                           update_q_points | update_JxW_values);

  const unsigned int   dofs_per_cell = fe.dofs_per_cell;
  const unsigned int   n_q_points    = quadrature_formula.size();

  FullMatrix<double>   cell_matrix (dofs_per_cell, dofs_per_cell);

  std::vector<unsigned int> local_dof_indices (dofs_per_cell);

  typename MGDoFHandler<dim>::cell_iterator cell = mg_dof_handler.begin(),
					    endc = mg_dof_handler.end();
  for (; cell!=endc; ++cell)
    {
      const unsigned int level = cell->level();
      cell_matrix = 0;

      fe_values.reinit (cell);

      for (unsigned int q_point=0; q_point<n_q_points; ++q_point)
	for (unsigned int i=0; i<dofs_per_cell; ++i)
	  {
	    for (unsigned int j=0; j<dofs_per_cell; ++j)
	      cell_matrix(i,j) += ((fe_values.shape_grad(i,q_point)
				    * fe_values.shape_grad(j,q_point))
				   + fe_values.shape_value(i,q_point)
				   * fe_values.shape_value(j,q_point))
				  * fe_values.JxW(q_point);
	  };


      cell->get_mg_dof_indices (local_dof_indices);
      for (unsigned int i=0; i<dofs_per_cell; ++i)
	{
	  for (unsigned int j=0; j<dofs_per_cell; ++j)
	    mg_matrices[level].add (local_dof_indices[i],
                                    local_dof_indices[j],
                                    cell_matrix(i,j));
	};
    };

}



template <int dim>
void LaplaceProblem<dim>::solve () 
{
  GrowingVectorMemory<>   vector_memory;

  MGTransferPrebuilt<Vector<double> > mg_transfer;
  mg_transfer.build_matrices(mg_dof_handler);

  FullMatrix<float> coarse_matrix;
  coarse_matrix.copy_from (mg_matrices[0]);
  MGCoarseGridHouseholder<float, Vector<double> > mg_coarse;
  mg_coarse.initialize(coarse_matrix);

  typedef PreconditionSOR<SparseMatrix<float> > RELAXATION;
  MGSmootherRelaxation<SparseMatrix<float>, RELAXATION, Vector<double> >
    mg_smoother(vector_memory);

  RELAXATION::AdditionalData smoother_data;
  mg_smoother.initialize(mg_matrices, smoother_data);
  
  mg_smoother.set_steps(2);
  mg_smoother.set_symmetric(true);

  MGMatrix<SparseMatrix<float>, Vector<double> >
    mg_matrix(&mg_matrices);
  Multigrid<Vector<double> > mg(mg_dof_handler,
				mg_matrix,
				mg_coarse,
				mg_transfer,
				mg_smoother,
				mg_smoother);
  PreconditionMG<dim, Vector<double>,
    MGTransferPrebuilt<Vector<double> > >
    preconditioner(mg_dof_handler, mg, mg_transfer);
  
  SolverControl           solver_control (1000, 1e-12);
  SolverCG<>              cg (solver_control);

  
  cg.solve (system_matrix, solution, system_rhs,
	    preconditioner);

  deallog << "   " << solver_control.last_step()
	    << " CG iterations needed to obtain convergence."
	    << std::endl;
}



template <int dim>
void LaplaceProblem<dim>::output_results (const unsigned int cycle) const
{
  DataOut<dim> data_out;

  data_out.attach_dof_handler (mg_dof_handler);
  data_out.add_data_vector (solution, "solution");
  data_out.build_patches ();

  std::ostringstream filename;
  filename << "solution-"
	   << cycle
	   << ".gnuplot";

  data_out.write_gnuplot (deallog.get_file_stream());
}



template <int dim>
void LaplaceProblem<dim>::run () 
{
  for (unsigned int cycle=0; cycle<6; ++cycle)
    {
      deallog << "Cycle " << cycle << ':' << std::endl;

      if (cycle == 0)
	{
	  GridGenerator::hyper_cube(triangulation);
	  
	}
      else
	triangulation.refine_global (1);

      deallog << "   Number of active cells: "
		<< triangulation.n_active_cells()
		<< std::endl
		<< "   Total number of cells: "
		<< triangulation.n_cells()
		<< std::endl;

      setup_system ();
      assemble_system ();
      assemble_multigrid ();
      solve ();
      output_results (cycle);
    };
}

    

int main () 
{
  deallog.depth_console (0);
  deallog << std::setprecision(2);
  logfile << std::setprecision(2);
  
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);  

  LaplaceProblem<2> laplace_problem_2d;
  laplace_problem_2d.run ();
  
  return 0;
}
