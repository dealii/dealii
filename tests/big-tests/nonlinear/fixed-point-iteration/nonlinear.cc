/* $Id$ */
/* Copyright W. Bangerth, University of Heidelberg, 1998 */




#include <lac/vector.h>
#include <grid/tria.h>
#include <dofs/dof_handler.h>
#include <grid/tria_accessor.h>
#include <dofs/dof_accessor.h>
#include <grid/tria_iterator.h>
#include <grid/tria_boundary.h>
#include <grid/grid_generator.h>
#include <grid/grid_refinement.h>
#include <numerics/data_out.h>
#include <base/function.h>
#include <fe/fe_lib.lagrange.h>
#include <base/quadrature_lib.h>
#include "../../problem_base.h"
#include <numerics/assembler.h>
#include <numerics/error_estimator.h>
#include <base/logstream.h>

#include <map>
#include <fstream>
#include <cmath>
#include <cstdlib>





template <int dim>
class RightHandSide :  public Function<dim> 
{
  public:
    virtual double value (const Point<dim> &p,
			  const unsigned int) const 
      {
	double x = 80;
	for (unsigned int d=0; d<dim; ++d)
	  if (p(d) < 0.5)
	    x *= -p(d);
	  else
	    x *= (1-p(d));
	
	return x;
      };
};



template <int dim>
class PoissonEquation :  public Equation<dim> {
  public:
    PoissonEquation (const Function<dim>  &rhs,
		     const Vector<double> &last_solution) :
		    Equation<dim>(1),
		    right_hand_side (rhs),
		    last_solution(last_solution) {};

    virtual void assemble (FullMatrix<double>  &cell_matrix,
			   Vector<double>      &rhs,
			   const FEValues<dim> &fe_values,
			   const DoFHandler<dim>::cell_iterator &cell) const;
    virtual void assemble (FullMatrix<double>  &cell_matrix,
			   const FEValues<dim> &fe_values,
			   const DoFHandler<dim>::cell_iterator &cell) const;
    virtual void assemble (Vector<double>      &rhs,
			   const FEValues<dim> &fe_values,
			   const DoFHandler<dim>::cell_iterator &cell) const;
  protected:
    const Function<dim>  &right_hand_side;
    const Vector<double> &last_solution;
};






template <int dim>
class NonlinearProblem : public ProblemBase<dim> {
  public:
    NonlinearProblem ();
    void run ();

  protected:
    Triangulation<dim> *tria;
    DoFHandler<dim>    *dof;

    Vector<double>      last_solution;
};




template <int dim>
void PoissonEquation<dim>::assemble (FullMatrix<double>  &cell_matrix,
				     Vector<double>      &rhs,
				     const FEValues<dim> &fe_values,
				     const DoFHandler<dim>::cell_iterator &) const {
  const vector<vector<Tensor<1,dim> > >&gradients = fe_values.get_shape_grads ();
  const FullMatrix<double> &values    = fe_values.get_shape_values ();
  vector<double>            rhs_values (fe_values.n_quadrature_points);
  const vector<double> &weights   = fe_values.get_JxW_values ();
  
  vector<Tensor<1,dim> >   last_solution_grads(fe_values.n_quadrature_points);
  fe_values.get_function_grads (last_solution, last_solution_grads);

  
  right_hand_side.value_list (fe_values.get_quadrature_points(), rhs_values);
   
  for (unsigned int point=0; point<fe_values.n_quadrature_points; ++point)
    for (unsigned int i=0; i<fe_values.dofs_per_cell; ++i) 
      {
	for (unsigned int j=0; j<fe_values.dofs_per_cell; ++j)
	  cell_matrix(i,j) += (gradients[i][point] *
			       gradients[j][point]) *
			      weights[point] /
			      sqrt(1+last_solution_grads[i]*last_solution_grads[i]);
	rhs(i) += values(i,point) *
		  rhs_values[point] *
		  weights[point];
      };
};




template <int dim>
void PoissonEquation<dim>::assemble (FullMatrix<double>  &,
				     const FEValues<dim> &,
				     const DoFHandler<dim>::cell_iterator &) const {
  Assert (false, typename Equation<dim>::ExcPureVirtualFunctionCalled());
};



template <int dim>
void PoissonEquation<dim>::assemble (Vector<double>      &,
				     const FEValues<dim> &,
				     const DoFHandler<dim>::cell_iterator &) const {
  Assert (false, typename Equation<dim>::ExcPureVirtualFunctionCalled());
};




template <int dim>
NonlinearProblem<dim>::NonlinearProblem () :
		tria(0), dof(0) {};



template <int dim>
void NonlinearProblem<dim>::run () {

				   // first reset everything to a virgin state
  clear ();
  
  tria = new Triangulation<dim>();
  dof = new DoFHandler<dim> (*tria);
  set_tria_and_dof (tria, dof);


  RightHandSide<dim>    rhs;
  ZeroFunction<dim>     boundary_values;
  StraightBoundary<dim> boundary;
  
  FEQ1<dim>                       fe;
  PoissonEquation<dim>            equation (rhs, last_solution);
  QGauss2<dim>                    quadrature;
  
  ProblemBase<dim>::FunctionMap dirichlet_bc;
  dirichlet_bc[0] = &boundary_values;


  GridGenerator::hyper_cube (*tria);
  tria->refine_global (4);

  for (unsigned int refinement_step=0; refinement_step<5; ++refinement_step)
    {
      cout << "Refinement step " << refinement_step << endl
	   << "  Grid has " << tria->n_active_cells() << " active cells." << endl;
  
      cout << "    Distributing dofs... "; 
      dof->distribute_dofs (fe);
      cout << dof->n_dofs() << " degrees of freedom." << endl;
      
				       // set the starting values for the iteration
				       // to a constant value of 1
      last_solution.reinit (dof->n_dofs());
      for (unsigned int i=0; i<dof->n_dofs(); ++i)
	last_solution(i) = 1;
  

				       // here comes the fixed point iteration
      for (unsigned int nonlinear_step=0; nonlinear_step<10; ++nonlinear_step)
	{
	  cout << "    Nonlinear step " << nonlinear_step << endl;
	  cout << "        Assembling matrices..." << endl;
	  assemble (equation, quadrature,
		    UpdateFlags(update_values | update_gradients |
				update_JxW_values | update_q_points),
		    dirichlet_bc);
	  
	  cout << "        Solving..." << endl;
	  solve ();
	  
	  if (nonlinear_step % 2 == 0)
	    {
	      string filename = "nonlinear.";
	      filename += ('0' + refinement_step);
	      filename += '.';
	      filename += ('0' + (nonlinear_step/2));
	      filename += ".gnuplot";
	      cout << "        Writing to file <" << filename << ">..." << endl;
	      
	      DataOut<dim> out;
	      ofstream gnuplot(filename.c_str());
	      fill_data (out);
	      out.build_patches ();
	      out.write_gnuplot (gnuplot);
	      gnuplot.close ();
	    };
	  
	  last_solution = solution;
	};

      Vector<float> error_indicator;
      KellyErrorEstimator<dim> ee;
      QSimpson<dim-1> eq;
      ee.estimate (*dof, eq,
		   KellyErrorEstimator<dim>::FunctionMap(),
		   solution,
		   error_indicator);
      GridRefinement::refine_and_coarsen_fixed_number (*tria, error_indicator, 0.3, 0);
      tria->execute_coarsening_and_refinement ();
    };
  
  
  delete dof;
  delete tria;
  
  cout << endl;
};




int main ()
{
  deallog.depth_console (0);
  
  NonlinearProblem<2> problem;
  problem.run ();
};
