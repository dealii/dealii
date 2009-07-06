/* $Id$ */
/* Author: Toby D. Young, Polish Academy of Sciences,             */
/*         Wolfgang Bangerth, Texas A&M University                */
/*    $Id$        */
/*                                                                */
/*    Copyright (C) 2009 by the deal.II authors */
/*                                                                */
/*    This file is subject to QPL and may not be  distributed     */
/*    without copyright and license information. Please refer     */
/*    to the file deal.II/doc/license.html for the  text  and     */
/*    further information on this license.                        */

                                // @sect3{Include files} We bundle the
                                // "usual" deal.II include files as we
                                // did in step-4:
#include <base/logstream.h>
#include <base/quadrature_lib.h>
#include <base/function.h>
#include <base/function_parser.h>
#include <base/parameter_handler.h>
#include <base/utilities.h>
#include <lac/full_matrix.h>
#include <lac/compressed_simple_sparsity_pattern.h>
#include <grid/tria.h>
#include <grid/grid_generator.h>
#include <grid/tria_accessor.h>
#include <grid/tria_iterator.h>
#include <dofs/dof_handler.h>
#include <dofs/dof_accessor.h>
#include <dofs/dof_tools.h>
#include <fe/fe_q.h>
#include <fe/fe_values.h>
#include <numerics/vectors.h>
#include <numerics/matrices.h>
#include <numerics/data_out.h>

                                // PETSc appears here because SLEPc
                                // depends on this library:
#include <lac/petsc_sparse_matrix.h>
#include <lac/petsc_vector.h>

                                // and then we need to actually import
                                // the interfaces for solvers that
                                // SLEPc provides:
#include <lac/slepc_solver.h>

                                // and some standard C++:
#include <fstream>
#include <iostream>

				// and the finally, as in previous
				// programs, we import all the deal.II
				// class and function names into the
				// global namespace:
using namespace dealii;

                                // @sect1{The
                                // <code>EigenvalueProblem</code>
                                // class template}

template <int dim>
class EigenvalueProblem 
{
public:
  EigenvalueProblem (const std::string &prm_file);
  void run ();
  
private:
  void make_grid_and_dofs ();
  void assemble_system ();
  void solve ();
  void output_results () const;
  
  Triangulation<dim> triangulation;
  FE_Q<dim>          fe;
  DoFHandler<dim>    dof_handler;
  
                                // These are data types pertaining to
                                // the generalized problem.
  PETScWrappers::SparseMatrix        stiffness_matrix, mass_matrix;
  std::vector<PETScWrappers::Vector> eigenfunctions;
  std::vector<double>                eigenvalues;   
  
  ParameterHandler parameters;
};

                                // @sect2{Implementation of the
                                // <code>EigenvalueProblem</code>
                                // class}

                                // @sect3{EigenvalueProblem::EigenvalueProblem}

template <int dim>
EigenvalueProblem<dim>::EigenvalueProblem (const std::string &prm_file)
  :
  fe (1),
  dof_handler (triangulation)
{

                                // Declare some of the needed data
                                // from a file; you can always change
                                // this!
  parameters.declare_entry ("Number of eigenvalues/eigenfunctions", "10",
			    Patterns::Integer (0, 100),
			    "The number of eigenvalues/eigenfunctions "
			    "to be computed.");
  parameters.declare_entry ("Potential", "0",
			    Patterns::Anything(),
			    "A functional description of the potential.");
  
                                // Entries are declared, so now we
                                // read them in...
  parameters.read_input (prm_file);
}


                                // @sect3{EigenvalueProblem::make_grid_and_dofs}
template <int dim>
void EigenvalueProblem<dim>::make_grid_and_dofs ()
{
  GridGenerator::hyper_cube (triangulation, -0.5, 0.5);
  triangulation.refine_global (5);
  dof_handler.distribute_dofs (fe);

  CompressedSimpleSparsityPattern csp (dof_handler.n_dofs(),
       				       dof_handler.n_dofs());
  DoFTools::make_sparsity_pattern (dof_handler, csp);
  csp.compress ();

  // What is going on here?

  // This does not work!
  //   stiffness_matrix.reinit (csp);
  //   mass_matrix.reinit (csp);
  
  // But this does... TODO: Fix it!
  stiffness_matrix.reinit (dof_handler.n_dofs(), dof_handler.n_dofs(),
  			   dof_handler.max_couplings_between_dofs());
  mass_matrix.reinit (dof_handler.n_dofs(), dof_handler.n_dofs(),
  		      dof_handler.max_couplings_between_dofs());

                                // with this done we stream-out the
                                // sparsity pattern
  std::ofstream out ("constrained_sparsity_pattern.gpl");
  csp.print_gnuplot (out);
  
                                // The next step is to take care of
                                // the eigenspectrum. In this case,
                                // the outputs are eigenfunctions and
                                // eigenvalues. Set the collective
                                // eigenfunction block to be as big as
                                // we wanted!
  eigenfunctions
    .resize (parameters.get_integer ("Number of eigenvalues/eigenfunctions"));
  for (unsigned int i=0; i<eigenfunctions.size (); ++i)
    eigenfunctions[i].reinit (dof_handler.n_dofs ());

                                // and do the same for the eigenvalue
                                // output, which had better be the
                                // same size as the eigenfunction
                                // block
  eigenvalues
    .resize (eigenfunctions.size ());

}


                                // @sect3{EigenvalueProblem::assemble_system}
template <int dim>
void EigenvalueProblem<dim>::assemble_system () 
{  
  QGauss<dim>   quadrature_formula(2);
  
  FEValues<dim> fe_values (fe, quadrature_formula, 
			   update_values | update_gradients |
                           update_quadrature_points | update_JxW_values);
  
  const unsigned int dofs_per_cell = fe.dofs_per_cell;
  const unsigned int n_q_points    = quadrature_formula.size();
  
  FullMatrix<double> cell_stiffness_matrix (dofs_per_cell, dofs_per_cell);
  FullMatrix<double> cell_mass_matrix (dofs_per_cell, dofs_per_cell);
  
  std::vector<unsigned int> local_dof_indices (dofs_per_cell);
  
  FunctionParser<dim> potential;
  potential.initialize (FunctionParser<dim>::default_variable_names (),
			parameters.get ("Potential"),
			typename FunctionParser<dim>::ConstMap());

                                // Here we call the actual
                                // quadrature-point values of the
                                // potential defined in the prm file
                                // which we read in earlier
  std::vector<double> potential_values (n_q_points);
  
                                // Initialize objects denoting zero
                                // boundary constraints for the
                                // present grid.
  ConstraintMatrix constraints;
  constraints.clear ();
  DoFTools::make_zero_boundary_constraints (dof_handler, constraints);
  constraints.close ();
  
  typename DoFHandler<dim>::active_cell_iterator
    cell = dof_handler.begin_active (),
    endc = dof_handler.end ();
  for (; cell!=endc; ++cell)
    {
      fe_values.reinit (cell);
      cell_stiffness_matrix = 0;
      cell_mass_matrix      = 0;
      
      potential.value_list (fe_values.get_quadrature_points(),
			    potential_values);
      
      for (unsigned int q_point=0; q_point<n_q_points; 
	   ++q_point)
	for (unsigned int i=0; i<dofs_per_cell; 
	     ++i)
	  for (unsigned int j=0; j<dofs_per_cell; 
	       ++j)
	    {
	      cell_stiffness_matrix (i, j)
		+= (fe_values.shape_grad (i, q_point) *
		    0.5                               *
		    fe_values.shape_grad (j, q_point) 
		    + 
		    fe_values.shape_value (i, q_point) *
		    //		     potential_values[q_point] *
		    0 * // infinite potential well
		    fe_values.shape_value (j, q_point)
		    ) * fe_values.JxW (q_point);
	      
	      cell_mass_matrix (i, j)
		+= (fe_values.shape_value (i, q_point) *
		    fe_values.shape_value (j, q_point) 
		    ) * fe_values.JxW (q_point);
	    }

	                        // Now we have the local system we
	                        // transfer it into the global objects
	                        // and take care of zero boundary
	                        // constraints,
      cell->get_dof_indices (local_dof_indices);
      
      constraints
	.distribute_local_to_global (cell_stiffness_matrix,
				     local_dof_indices,
				     stiffness_matrix);
      constraints
	.distribute_local_to_global (cell_mass_matrix,
				     local_dof_indices,
				     mass_matrix);
    }

                                // and finally, set the matrices and
                                // vectors in an assembled state.
  stiffness_matrix.compress();
  mass_matrix.compress();

				// make sure that the diagonal entries
				// of constrained degrees of freedom
				// are non-zero to ensure that the
				// matrix is actually invertible
  //   for (unsigned int i=0; i<dof_handler.n_dofs(); ++i)
  //     if (constraints.is_constrained(i))
  //       {
  // 	stiffness_matrix.set (i, i, 1);
  // 	mass_matrix.set (i, i, 1);
  //       }
  
                                // finally set the matrices in an
                                // assembled state so SLEPc likes:
  //   stiffness_matrix.compress ();
  //   mass_matrix.compress ();  
}


                                // @sect3{EigenvalueProblem::solve}
                                // Now that the system is set up, here
                                // is a good time to actually solve
                                // the problem: As with other examples
                                // this is done using a "solve"
                                // routine
template <int dim>
void EigenvalueProblem<dim>::solve () 
{
                                // We start by assigning the accuracy
                                // to which we would like to solve the
                                // system,
  SolverControl solver_control (1000, 1e-6);

                                // and assign our solver of
                                // choice. Here we want to use the
                                // Krylov-Schur solver, which is
                                // pretty darn fast and robust:
  SLEPcWrappers::SolverKrylovSchur eigensolver (solver_control);

                                // Lets assign the solver which part
                                // of the spectrum we want to solve
  eigensolver.set_which_eigenpairs (EPS_SMALLEST_MAGNITUDE);

                                // Finally, we actually solve the
                                // generalized eigenproblem:
  eigensolver.solve (stiffness_matrix, mass_matrix, 
   		     eigenvalues, eigenfunctions, 
   		     eigenfunctions.size());
}


                                // @sect3{EigenvalueProblem::output_results}
template <int dim>
void EigenvalueProblem<dim>::output_results () const
{
  DataOut<dim> data_out;

  data_out.attach_dof_handler (dof_handler);

  for (unsigned int i=0; i<eigenfunctions.size(); ++i)
    data_out.add_data_vector (eigenfunctions[i],
			      std::string("solution") +
			      Utilities::int_to_string(i));

                                // How does this work?
  Vector<double> projected_potential (dof_handler.n_dofs());
  FunctionParser<dim> potential;
  potential.initialize (FunctionParser<dim>::default_variable_names (),
			parameters.get ("Potential"),
			typename FunctionParser<dim>::ConstMap());
  VectorTools::interpolate (dof_handler, potential, projected_potential);
  data_out.add_data_vector (projected_potential, "interpolated_potential");
  
  data_out.build_patches ();

  std::ofstream output ("eigenvectors.vtk");
  data_out.write_vtk (output);

  for (unsigned int i=0; i<eigenvalues.size(); ++i)
    std::cout << std::endl 
	      << "      eigenvalue " << i 
	      << " : " << eigenvalues[i];

}


                                 // @sect3{EigenvalueProblem::run}

                                 // This is the function which has the
				 // top-level control over
				 // everything. It is very similar as
				 // for the previous examples.
template <int dim>
void EigenvalueProblem<dim>::run () 
{
  std::cout << "Solving problem in " << dim << " space dimensions." << std::endl;
  
  make_grid_and_dofs ();

                                 // While we are here, lets count the
                                 // number of active cells and degrees
                                 // of freedom like we always do.
  std::cout << "   Number of active cells:       "
	    << triangulation.n_active_cells()
	    << std::endl
	    << "   Number of degrees of freedom: "
	    << dof_handler.n_dofs();
  
  assemble_system ();
  solve ();
  output_results ();
}


                                 // @sect3{The <code>main</code> function}
int main (int argc, char **argv) 
{
  try
    {

                                 // Here is another difference from
                                 // other steps: We initialize the
                                 // SLEPc work space which inherently
                                 // initializes the PETSc work space
      SlepcInitialize (&argc,&argv,0,0);

      {
	deallog.depth_console (0);
	
	EigenvalueProblem<2> problem ("step-36.prm");
	problem.run ();
      }

                                 // and then unitialize the SLEPc
                                 // work space when the job is done:
      SlepcFinalize ();
    }

                                 // or panic if something goes wrong.
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
  
                                 // or show that we happy and didn't
                                 // crap out on the calculation...
  std::cout << std::endl 
	    << "Completed" 
	    << std::endl;

  return 0;
}
