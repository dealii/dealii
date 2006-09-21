//----------------------------  hp_hanging_node_constraints_02.cc  ---------------------------
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
//----------------------------  hp_hanging_node_constraints_02.cc  ---------------------------


// modified step-3 to be a L_2 projection. This function
// L_2-projects a bilinear, biquadratic or bicubic function
// onto a square shaped domain which is filled with
// hp elements. Afterwards the L_2 error is computed, which
// should be in the range of the FP approximation order.

#include <base/logstream.h>
#include <grid/tria.h>
#include <dofs/dof_handler.h>
#include <dofs/dof_constraints.h>
#include <grid/grid_generator.h>

#include <grid/tria_accessor.h>
#include <grid/tria_iterator.h>
#include <dofs/dof_accessor.h>

#include <fe/fe_q.h>

#include <dofs/dof_tools.h>

#include <fe/hp_fe_values.h>
#include <base/quadrature_lib.h>

#include <base/function.h>
#include <numerics/vectors.h>
#include <numerics/matrices.h>

#include <lac/vector.h>
#include <lac/full_matrix.h>
#include <lac/sparse_matrix.h>
#include <lac/solver_cg.h>
#include <lac/precondition.h>

#include <numerics/data_out.h>
#include <fstream>
#include <iostream>

std::ofstream logfile("hp_hanging_nodes_02/output");


template <int dim>
class Linear : public Function<dim> 
{
  public:
    Linear () : Function<dim>() {};
    
    virtual double value (const Point<dim>   &p,
			  const unsigned int  component = 0) const;    
};



template <int dim>
double Linear<dim>::value (const Point<dim> &p,
			   const unsigned int) const 
{
  double res = p(0);
  for (unsigned int d = 1; d < dim; ++d)
    res *= p(d);
  
  return res;  
}



template <int dim>
class Quadratic : public Function<dim> 
{
  public:
    Quadratic () : Function<dim>() {};
    
    virtual double value (const Point<dim>   &p,
			  const unsigned int  component = 0) const;    
};



template <int dim>
double Quadratic<dim>::value (const Point<dim> &p,
			   const unsigned int) const 
{
  double res = p(0) * p(0);
  for (unsigned int d = 1; d < dim; ++d)
    res *= p(d) * p(d);
  
  return res;  
}


template <int dim>
class Cubic : public Function<dim> 
{
  public:
    Cubic () : Function<dim>() {};
    
    virtual double value (const Point<dim>   &p,
			  const unsigned int  component = 0) const;    
};



template <int dim>
double Cubic<dim>::value (const Point<dim> &p,
			   const unsigned int) const 
{
  double res = p(0) * p(0) * p(0);
  for (unsigned int d = 1; d < dim; ++d)
    res *= p(d) * p(d) * p(d);
  
  return res;  
}


template <int dim>
class LaplaceProblem 
{
  public:
    LaplaceProblem ();

    void run (const Function<dim> &f_test,
	      bool random,
	      unsigned int *indx);
    
  private:
    void make_grid_and_dofs (const bool random_p);
    void assemble_system (const Function<dim> &f_test);
    void solve ();
    void eval_error (const Function<dim> &f_test);
    void output_results () const;

    Triangulation<dim>     triangulation;
    hp::FECollection<dim>              fe;
    hp::DoFHandler<dim>        dof_handler;

    SparsityPattern      sparsity_pattern;
    SparseMatrix<double> system_matrix;

				     // Although we do not have h-refinement,
				     // hanging nodes will inevitably appear
				     // due to different polynomial degrees.
    ConstraintMatrix     hanging_node_constraints;
    
    Vector<double>       solution;
    Vector<double>       system_rhs;
};


template <int dim>
LaplaceProblem<dim>::LaplaceProblem () :
		dof_handler (triangulation)
{}


template <int dim>
void LaplaceProblem<dim>::make_grid_and_dofs (const bool random_p)
{
  GridGenerator::hyper_cube (triangulation, -1, 1);
  triangulation.refine_global (5-dim);
  deallog << "Number of active cells: "
	    << triangulation.n_active_cells()
	    << std::endl;
  deallog << "Total number of cells: "
	    << triangulation.n_cells()
	    << std::endl;

				   // Now to the p-Method. Assign
				   // random active_fe_indices to the
				   // different cells.
  typename hp::DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active (),
						     endc = dof_handler.end ();
  if (random_p)
    {      
      for (; cell != endc; ++cell)
	{      
	  cell->set_active_fe_index ((int)(4.0 * (double) random () / (double) RAND_MAX));
	}
    }
  else
    {      
      unsigned int cell_no = 0;  
      for (; cell != endc; ++cell)
	{
	  if (cell_no >= triangulation.n_active_cells () / 2)
	    cell->set_active_fe_index (1);
	  else
	    cell->set_active_fe_index (0);
	}  
    }

  
  dof_handler.distribute_dofs (fe);
  deallog << "Number of degrees of freedom: "
	    << dof_handler.n_dofs()
	    << std::endl;

  solution.reinit (dof_handler.n_dofs());
  system_rhs.reinit (dof_handler.n_dofs());

				   // Create sparsity pattern.  
  sparsity_pattern.reinit (dof_handler.n_dofs(),
			   dof_handler.n_dofs(),
			   dof_handler.max_couplings_between_dofs());
  DoFTools::make_sparsity_pattern (dof_handler, sparsity_pattern);

				   // Create constraints which stem from
				   // the different polynomial degrees on
				   // the different elements.
  hanging_node_constraints.clear ();
  DoFTools::make_hanging_node_constraints (dof_handler,
					   hanging_node_constraints);

  hanging_node_constraints.print (deallog.get_file_stream ());

  hanging_node_constraints.close ();
  hanging_node_constraints.condense (sparsity_pattern);
  
  sparsity_pattern.compress();
  system_matrix.reinit (sparsity_pattern);
}



template <int dim>
void LaplaceProblem<dim>::assemble_system (const Function<dim> &f_test) 
{
  hp::QCollection<dim>  quadrature_formula(QGauss<dim>(6));
  hp::FEValues<dim> x_fe_values (fe, quadrature_formula, 
			 update_values | update_gradients |
			       update_JxW_values | update_q_points);
  
  const unsigned int   max_dofs_per_cell = fe.max_dofs_per_cell ();
  const unsigned int   n_q_points    = quadrature_formula[0].n_quadrature_points;

  FullMatrix<double>   cell_matrix (max_dofs_per_cell, max_dofs_per_cell);
  Vector<double>       cell_rhs (max_dofs_per_cell);

  std::vector<unsigned int> local_dof_indices (max_dofs_per_cell);

  typename hp::DoFHandler<dim>::active_cell_iterator
    cell = dof_handler.begin_active(),
    endc = dof_handler.end();
  for (; cell!=endc; ++cell)
    {
      x_fe_values.reinit (cell);

      const FEValues<dim> &fe_values = x_fe_values.get_present_fe_values();
      
      cell_matrix = 0;
      cell_rhs = 0;

      const unsigned int dofs_per_cell = cell->get_fe ().dofs_per_cell;
      
      for (unsigned int i=0; i<dofs_per_cell; ++i)
	for (unsigned int j=0; j<dofs_per_cell; ++j)
	  for (unsigned int q_point=0; q_point<n_q_points; ++q_point)
	    cell_matrix(i,j) += (fe_values.shape_value (i, q_point) *
				 fe_values.shape_value (j, q_point) *
				 fe_values.JxW (q_point));

      for (unsigned int i=0; i<dofs_per_cell; ++i)
	for (unsigned int q_point=0; q_point<n_q_points; ++q_point)
	  {	    
	    Point<dim> p_q = fe_values.quadrature_point (q_point);
	    double f = f_test.value (p_q);	    
	    cell_rhs(i) += (fe_values.shape_value (i, q_point) *
			    f *
			    fe_values.JxW (q_point));
	  }      
	    
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

				   // Include hanging nodes.
  hanging_node_constraints.condense (system_matrix);
  hanging_node_constraints.condense (system_rhs);
}



template <int dim>
void LaplaceProblem<dim>::solve () 
{
  SolverControl           solver_control (1000, 1e-12);
  SolverCG<>              cg (solver_control);

  cg.solve (system_matrix, solution, system_rhs,
	    PreconditionIdentity());

  hanging_node_constraints.distribute (solution);
}



template <int dim>
void LaplaceProblem<dim>::eval_error (const Function<dim> &f_test) 
{
  hp::QCollection<dim>  quadrature(QGauss<dim>(6));
  Vector<double> cellwise_errors (triangulation.n_active_cells());
  VectorTools::integrate_difference (dof_handler, solution, f_test,
                                     cellwise_errors, quadrature,
                                     VectorTools::L2_norm);
  const double p_l2_error = cellwise_errors.l2_norm();

  deallog << "L2_Error : " << p_l2_error << std::endl;
}



template <int dim>
void LaplaceProblem<dim>::output_results () const
{
  DataOut<dim,hp::DoFHandler<dim> > data_out;
  data_out.attach_dof_handler (dof_handler);
  data_out.add_data_vector (solution, "solution");
  data_out.build_patches ();

  data_out.write_gnuplot (deallog.get_file_stream());
}



template <int dim>
void LaplaceProblem<dim>::run (const Function<dim> &f_test,
			       bool random,
			       unsigned int *indx) 
{
  FE_Q<dim> fe_1 (indx[0]),
    fe_2 (indx[1]),
    fe_3 (indx[2]),
    fe_4 (indx[3]);

  fe.push_back (fe_1);
  fe.push_back (fe_2);
  fe.push_back (fe_3);
  fe.push_back (fe_4);
  
  make_grid_and_dofs (random);
  assemble_system (f_test);
  solve ();
  eval_error (f_test);
  output_results ();
}


template <int dim>
void run_test (const Function<dim> &f_test,
	       unsigned int *indx)
{
  LaplaceProblem<dim> laplace_problem_1;
  laplace_problem_1.run (f_test, true, indx);

  LaplaceProblem<dim> laplace_problem_2;
  laplace_problem_2.run (f_test, false, indx);
}



int main () 
{
  logfile.precision(2);
  
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);  

  unsigned int index[] = 
    {
	  1,2,3,4,5,6,7
    };
  
  
  Linear<2> test2d_1;
  Quadratic<2> test2d_2;
  Cubic<2> test2d_3;

  Linear<3> test3d_1;
  Quadratic<3> test3d_2;
  Cubic<3> test3d_3;

  deallog << "Testing Order 1" << std::endl;
  run_test<2> (test2d_1, &(index[0]));
  run_test<3> (test3d_1, &(index[0]));

  deallog << "Testing Order 2" << std::endl;
  run_test<2> (test2d_2, &(index[1]));
  run_test<3> (test3d_2, &(index[1]));

  deallog << "Testing Order 3" << std::endl;
  run_test<2> (test2d_3, &(index[2]));
  run_test<3> (test3d_3, &(index[2]));
      
  return 0;
}
