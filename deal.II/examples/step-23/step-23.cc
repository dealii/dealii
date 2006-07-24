/* $Id: step-4.cc,v 1.34 2006/02/06 21:33:10 wolf Exp $ */
/* Author: Wolfgang Bangerth, University of Heidelberg, 1999 */

/*    $Id: step-4.cc,v 1.34 2006/02/06 21:33:10 wolf Exp $       */
/*    Version: $Name:  $                                          */
/*                                                                */
/*    Copyright (C) 1999, 2000, 2001, 2002, 2003, 2004, 2005, 2006 by the deal.II authors */
/*                                                                */
/*    This file is subject to QPL and may not be  distributed     */
/*    without copyright and license information. Please refer     */
/*    to the file deal.II/doc/license.html for the  text  and     */
/*    further information on this license.                        */


#include <grid/tria.h>
#include <dofs/dof_handler.h>
#include <grid/grid_generator.h>
#include <grid/tria_accessor.h>
#include <grid/tria_iterator.h>
#include <dofs/dof_accessor.h>
#include <fe/fe_q.h>
#include <dofs/dof_tools.h>
#include <fe/fe_values.h>
#include <base/quadrature_lib.h>
#include <base/function.h>
#include <numerics/vectors.h>
#include <numerics/matrices.h>
#include <lac/vector.h>
#include <lac/full_matrix.h>
#include <lac/sparse_matrix.h>
#include <lac/solver_cg.h>
#include <lac/precondition.h>
#include <dofs/dof_constraints.h>

#include <numerics/matrices.h>
#include <numerics/vectors.h>

#include <numerics/data_out.h>
#include <fstream>
#include <iostream>
#include <sstream>

#include <base/logstream.h>



template <int dim>
class WaveEquationProblem 
{
  public:
    WaveEquationProblem ();
    void run ();
    
  private:
    void make_grid_and_dofs ();
    void assemble_system ();
    void solve_u ();
    void solve_v ();
    void output_results (const unsigned int timestep_number) const;

    Triangulation<dim>   triangulation;
    FE_Q<dim>            fe;
    DoFHandler<dim>      dof_handler;

    SparsityPattern      sparsity_pattern;
    SparseMatrix<double> system_matrix;
    SparseMatrix<double> mass_matrix;
    SparseMatrix<double> laplace_matrix;

    double time_step;
    double theta;

    Vector<double>       solution_u, solution_v;
    Vector<double>       old_solution_u, old_solution_v;
    Vector<double>       system_rhs;
};



template <int dim>
class RightHandSide : public Function<dim> 
{
  public:
    RightHandSide () : Function<dim>() {};
    
    virtual double value (const Point<dim>   &p,
			  const unsigned int  component = 0) const;
};


template <int dim>
class InitialValuesU : public Function<dim> 
{
  public:
    InitialValuesU () : Function<dim>() {};
    
    virtual double value (const Point<dim>   &p,
			  const unsigned int  component = 0) const;
};



template <int dim>
class BoundaryValues : public Function<dim> 
{
  public:
    BoundaryValues () : Function<dim>() {};
    
    virtual double value (const Point<dim>   &p,
			  const unsigned int  component = 0) const;
};




template <int dim>
double RightHandSide<dim>::value (const Point<dim> &/*p*/,
				  const unsigned int /*component*/) const 
{
//   if (get_time() <= 0.25)
//     if ((p[0] <=0) && (p[1] <= 0))
//       return 1;

  return 0;
}


template <int dim>
double InitialValuesU<dim>::value (const Point<dim> &p,
				  const unsigned int /*component*/) const 
{
  //  return std::sqrt(p.square()) * std::exp (-p.square()) / 3;
  if ((p[0] <=0) && (p[1] <= 0))
    return 1;

  return 0;}


template <int dim>
double BoundaryValues<dim>::value (const Point<dim> &/*p*/,
				   const unsigned int /*component*/) const 
{
  return 0;
}







template <int dim>
WaveEquationProblem<dim>::WaveEquationProblem () :
                fe (1),
		dof_handler (triangulation),
		time_step (1./64),
		theta (0.5)
{}



template <int dim>
void WaveEquationProblem<dim>::make_grid_and_dofs ()
{
  GridGenerator::hyper_cube (triangulation, -1, 1);
  triangulation.refine_global (7);
  
  std::cout << "   Number of active cells: "
	    << triangulation.n_active_cells()
	    << std::endl
	    << "   Total number of cells: "
	    << triangulation.n_cells()
	    << std::endl;

  dof_handler.distribute_dofs (fe);

  std::cout << "   Number of degrees of freedom: "
	    << dof_handler.n_dofs()
	    << std::endl;

  sparsity_pattern.reinit (dof_handler.n_dofs(),
			   dof_handler.n_dofs(),
			   dof_handler.max_couplings_between_dofs());
  DoFTools::make_sparsity_pattern (dof_handler, sparsity_pattern);
  sparsity_pattern.compress();

  system_matrix.reinit (sparsity_pattern);
  mass_matrix.reinit (sparsity_pattern);
  laplace_matrix.reinit (sparsity_pattern);

  solution_u.reinit (dof_handler.n_dofs());
  solution_v.reinit (dof_handler.n_dofs());
  old_solution_u.reinit (dof_handler.n_dofs());
  old_solution_v.reinit (dof_handler.n_dofs());
  system_rhs.reinit (dof_handler.n_dofs());
}



template <int dim>
void WaveEquationProblem<dim>::assemble_system () 
{  
  MatrixCreator::create_mass_matrix (dof_handler, QGauss<dim>(3),
				     mass_matrix);
  MatrixCreator::create_laplace_matrix (dof_handler, QGauss<dim>(3),
					laplace_matrix);
  
  system_matrix.copy_from (mass_matrix);
  system_matrix.add (theta * theta * time_step * time_step, laplace_matrix);
}



template <int dim>
void WaveEquationProblem<dim>::solve_u () 
{
  SolverControl           solver_control (1000, 1e-8*system_rhs.l2_norm());
  SolverCG<>              cg (solver_control);
  cg.solve (system_matrix, solution_u, system_rhs,
	    PreconditionIdentity());

  std::cout << "   u-equation: " << solver_control.last_step()
	    << " CG iterations."
	    << std::endl;
}


template <int dim>
void WaveEquationProblem<dim>::solve_v () 
{
  SolverControl           solver_control (1000, 1e-8*system_rhs.l2_norm());
  SolverCG<>              cg (solver_control);
  cg.solve (mass_matrix, solution_v, system_rhs,
	    PreconditionIdentity());

  std::cout << "   v-equation: " << solver_control.last_step()
	    << " CG iterations."
	    << std::endl;
}



template <int dim>
void WaveEquationProblem<dim>::output_results (const unsigned int timestep_number) const
{
  DataOut<dim> data_out;

  data_out.attach_dof_handler (dof_handler);
  data_out.add_data_vector (solution_u, "U");
  data_out.add_data_vector (solution_v, "V");

  data_out.build_patches ();

  std::ostringstream filename;
  filename << "solution-"
	   << timestep_number
	   << ".gnuplot";
  std::ofstream output (filename.str().c_str());
  data_out.write_gnuplot (output);
}




template <int dim>
void WaveEquationProblem<dim>::run () 
{
  std::cout << "Solving problem in " << dim << " space dimensions." << std::endl;
  
  make_grid_and_dofs();
  assemble_system ();

  ConstraintMatrix constraints;
  constraints.close();
  VectorTools::project (dof_handler, constraints, QGauss<dim>(3),
			InitialValuesU<dim>(),
			old_solution_u);
  VectorTools::project (dof_handler, constraints, QGauss<dim>(3),
			ZeroFunction<dim>(),
			old_solution_v);

  unsigned int timestep_number = 1;
  for (double time = time_step; time<=5; time+=time_step, ++timestep_number)
    {
      std::cout << "Time step " << timestep_number
		<< " at t=" << time
		<< std::endl;
      
      Vector<double> tmp (solution_u.size());

      mass_matrix.vmult (system_rhs, old_solution_u);

      mass_matrix.vmult (tmp, old_solution_v);
      system_rhs.add (time_step, tmp);

      laplace_matrix.vmult (tmp, old_solution_u);
      system_rhs.add (-theta * (1-theta) * time_step * time_step, tmp);

      RightHandSide<dim> rhs_function;
      rhs_function.set_time (time);
      VectorTools::create_right_hand_side (dof_handler, QGauss<dim>(2),
					   rhs_function, tmp);
      system_rhs.add (theta * theta * time_step * time_step, tmp);

      rhs_function.set_time (time-time_step);
      VectorTools::create_right_hand_side (dof_handler, QGauss<dim>(2),
					   rhs_function, tmp);
      system_rhs.add (theta * (1-theta) * time_step * time_step, tmp);


      std::map<unsigned int,double> boundary_values;
      VectorTools::interpolate_boundary_values (dof_handler,
						0,
						BoundaryValues<dim>(),
						boundary_values);
      MatrixTools::apply_boundary_values (boundary_values,
					  system_matrix,
					  solution_u,
					  system_rhs);
      solve_u ();


      laplace_matrix.vmult (system_rhs, solution_u);
      system_rhs *= -theta * time_step;

      mass_matrix.vmult (tmp, old_solution_v);
      system_rhs += tmp;

      laplace_matrix.vmult (tmp, old_solution_u);
      system_rhs.add (-time_step * (1-theta), tmp);

      rhs_function.set_time (time);
      VectorTools::create_right_hand_side (dof_handler, QGauss<dim>(2),
					   rhs_function, tmp);
      system_rhs.add (theta * time_step, tmp);

      rhs_function.set_time (time-time_step);
      VectorTools::create_right_hand_side (dof_handler, QGauss<dim>(2),
					   rhs_function, tmp);
      system_rhs.add ((1-theta) * time_step, tmp);

      solve_v ();

      output_results (timestep_number);

      old_solution_u = solution_u;
      old_solution_v = solution_v;
    }
}



int main () 
{
  deallog.depth_console (0);
  {
    WaveEquationProblem<2> wave_equation_problem_2d;
    wave_equation_problem_2d.run ();
  }
  
  return 0;
}
