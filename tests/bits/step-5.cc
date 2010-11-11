//----------------------------  step-5.cc  ---------------------------
//    $Id$
//    Version: $Name$ 
//
//    Copyright (C) 2005, 2006, 2007, 2008, 2010 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  step-5.cc  ---------------------------


// a un-hp-ified version of hp/step-5


#include "../tests.h"
#include <base/logstream.h>
#include <fstream>
std::ofstream logfile("step-5/output");


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
#include <hp/dof_handler.h>
#include <grid/tria_accessor.h>
#include <grid/tria_iterator.h>
#include <dofs/dof_accessor.h>
#include <dofs/dof_tools.h>
#include <fe/fe_q.h>
#include <hp/fe_values.h>
#include <numerics/vectors.h>
#include <numerics/matrices.h>
#include <numerics/data_out.h>

#include <grid/grid_in.h>

#include <grid/tria_boundary_lib.h>

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
    void solve ();
    void output_results (const unsigned int cycle) const;

    Triangulation<dim>   triangulation;
    FE_Q<dim>            fe;
    DoFHandler<dim>      dof_handler;

    SparsityPattern      sparsity_pattern;
    SparseMatrix<double> system_matrix;

    Vector<double>       solution;
    Vector<double>       system_rhs;
};



template <int dim>
class Coefficient : public Function<dim> 
{
  public:
    Coefficient ()  : Function<dim>() {}
    
    virtual double value (const Point<dim>   &p,
			  const unsigned int  component = 0) const;
    
    virtual void value_list (const std::vector<Point<dim> > &points,
			     std::vector<double>            &values,
			     const unsigned int              component = 0) const;
};



template <int dim>
double Coefficient<dim>::value (const Point<dim> &p,
				const unsigned int /*component*/) const 
{
  if (p.square() < 0.5*0.5)
    return 20;
  else
    return 1;
}




template <int dim>
void Coefficient<dim>::value_list (const std::vector<Point<dim> > &points,
				   std::vector<double>            &values,
				   const unsigned int              component) const 
{
  Assert (values.size() == points.size(), 
	  ExcDimensionMismatch (values.size(), points.size()));
  
  Assert (component == 0, 
	  ExcIndexRange (component, 0, 1));

  const unsigned int n_points = points.size();

  for (unsigned int i=0; i<n_points; ++i)
    {
      if (points[i].square() < 0.5*0.5)
	values[i] = 20;
      else
	values[i] = 1;
    }
}




template <int dim>
LaplaceProblem<dim>::LaplaceProblem () :
                fe (1),
		dof_handler (triangulation)
{}




template <int dim>
void LaplaceProblem<dim>::setup_system ()
{
  dof_handler.distribute_dofs (fe);

  deallog << "   Number of degrees of freedom: "
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
  QGauss<dim>  quadrature_formula (2);

  FEValues<dim> x_fe_values (fe, quadrature_formula, 
			   update_values    |  update_gradients |
			   update_q_points  |  update_JxW_values);

  const unsigned int   dofs_per_cell = fe.dofs_per_cell;
  const unsigned int   n_q_points    = quadrature_formula.size();

  FullMatrix<double>   cell_matrix (dofs_per_cell, dofs_per_cell);
  Vector<double>       cell_rhs (dofs_per_cell);

  std::vector<unsigned int> local_dof_indices (dofs_per_cell);

  const Coefficient<dim> coefficient;
  std::vector<double>    coefficient_values (n_q_points);

  typename DoFHandler<dim>::active_cell_iterator
    cell = dof_handler.begin_active(),
    endc = dof_handler.end();
  for (; cell!=endc; ++cell)
    {
      cell_matrix = 0;
      cell_rhs = 0;

      x_fe_values.reinit (cell);
      const FEValues<2> &fe_values = x_fe_values.get_present_fe_values();

      coefficient.value_list (fe_values.get_quadrature_points(),
			      coefficient_values);
      
      for (unsigned int q_point=0; q_point<n_q_points; ++q_point)
	for (unsigned int i=0; i<dofs_per_cell; ++i)
	  {
	    for (unsigned int j=0; j<dofs_per_cell; ++j)
	      cell_matrix(i,j) += (coefficient_values[q_point] *
				   fe_values.shape_grad(i,q_point) *
				   fe_values.shape_grad(j,q_point) *
				   fe_values.JxW(q_point));

	    cell_rhs(i) += (fe_values.shape_value(i,q_point) *
			    1.0 *
			    fe_values.JxW(q_point));
	  }


      cell->get_dof_indices (local_dof_indices);
      for (unsigned int i=0; i<dofs_per_cell; ++i)
	{
	  for (unsigned int j=0; j<dofs_per_cell; ++j)
	    system_matrix.add (local_dof_indices[i],
			       local_dof_indices[j],
			       cell_matrix(i,j));
	  
	  system_rhs(local_dof_indices[i]) += cell_rhs(i);
	}
    }

  std::map<unsigned int,double> boundary_values;
  VectorTools::interpolate_boundary_values (dof_handler,
					    0,
					    ZeroFunction<dim>(),
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

  PreconditionSSOR<> preconditioner;
  preconditioner.initialize(system_matrix, 1.2);

  cg.solve (system_matrix, solution, system_rhs,
	    preconditioner);

  deallog << "   " << solver_control.last_step()
	    << " CG iterations needed to obtain convergence."
	    << std::endl;
}



template <int dim>
void LaplaceProblem<dim>::output_results (const unsigned int cycle) const
{
                                   // reduce output a bit
  if (cycle >= 2)
    return;
  
  DataOut<dim,DoFHandler<dim> > data_out;

  data_out.attach_dof_handler (dof_handler);
  data_out.add_data_vector (solution, "solution");

  data_out.build_patches ();

  DataOutBase::EpsFlags eps_flags;
  eps_flags.z_scaling = 4;
  eps_flags.azimut_angle = 40;
  eps_flags.turn_angle   = 10;
  data_out.set_flags (eps_flags);

  data_out.write_eps (deallog.get_file_stream());
}




template <int dim>
void LaplaceProblem<dim>::run () 
{
  for (unsigned int cycle=0; cycle<6; ++cycle)
    {
      deallog << "Cycle " << cycle << ':' << std::endl;

      if (cycle != 0)
	triangulation.refine_global (1);
      else
	{
	  GridIn<dim> grid_in;
	  grid_in.attach_triangulation (triangulation);
	  std::ifstream input_file("step-5/circle-grid.inp");
	  Assert (dim==2, ExcInternalError());
	  
	  grid_in.read_ucd (input_file);

	  static const HyperBallBoundary<dim> boundary;
	  triangulation.set_boundary (0, boundary);
	}

      deallog << "   Number of active cells: "
		<< triangulation.n_active_cells()
		<< std::endl
		<< "   Total number of cells: "
		<< triangulation.n_cells()
		<< std::endl;

      setup_system ();
      assemble_system ();
      solve ();
      output_results (cycle);
    }
}



int main () 
{
  deallog << std::setprecision(2);
  logfile << std::setprecision(2);
    
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);  

  deallog.depth_console (0);

  LaplaceProblem<2> laplace_problem_2d;
  laplace_problem_2d.run ();

/*  
  Coefficient<2>    coefficient;
  std::vector<Point<2> > points (2);
  std::vector<double>    coefficient_values (1);
  coefficient.value_list (points, coefficient_values);
*/
  
  return 0;
}
