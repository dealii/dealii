/* $Id$ */
/* Author: Toby Young, Polish Academy of Sciences,                */
/*         Wolfgang Bangerth, Texas A&M University                */
/*    $Id$        */
/*                                                                */
/*    Copyright (C) 2009 by the deal.II authors */
/*                                                                */
/*    This file is subject to QPL and may not be  distributed     */
/*    without copyright and license information. Please refer     */
/*    to the file deal.II/doc/license.html for the  text  and     */
/*    further information on this license.                        */

                                 // @sect3{Include files}

#include <base/logstream.h>
#include <base/quadrature_lib.h>
#include <base/function.h>
#include <base/function_parser.h>
#include <base/parameter_handler.h>
#include <base/utilities.h>
#include <lac/full_matrix.h>
#include <lac/petsc_sparse_matrix.h>
#include <lac/petsc_vector.h>
#include <lac/solver_cg.h>
#include <lac/precondition.h>
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

#include <fstream>
#include <iostream>


				 // The final step, as in previous
				 // programs, is to import all the
				 // deal.II class and function names
				 // into the global namespace:
using namespace dealii;

                                 // @sect3{The <code>EigenvalueProblem</code> class template}

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

    Triangulation<dim>   triangulation;
    FE_Q<dim>            fe;
    DoFHandler<dim>      dof_handler;

    PETScWrappers::SparseMatrix        stiffness_matrix, mass_matrix;
    std::vector<PETScWrappers::Vector> eigenfunctions;

    ParameterHandler parameters;
};



                                 // @sect3{Implementation of the <code>EigenvalueProblem</code> class}

                                 // @sect4{EigenvalueProblem::EigenvalueProblem}

template <int dim>
EigenvalueProblem<dim>::EigenvalueProblem (const std::string &prm_file)
		:
                fe (1),
		dof_handler (triangulation)
{
  parameters.declare_entry ("Number of eigenvalues/eigenfunctions", "10",
			    Patterns::Integer (0, 100),
			    "The number of eigenvalues/eigenfunctions "
			    "to be computed.");
  parameters.declare_entry ("Potential", "0",
			    Patterns::Anything(),
			    "A functional description of the potential.");

  parameters.read_input (prm_file);
}


                                 // @sect4{EigenvalueProblem::make_grid_and_dofs}
template <int dim>
void EigenvalueProblem<dim>::make_grid_and_dofs ()
{
  GridGenerator::hyper_cube (triangulation, -1, 1);
  triangulation.refine_global (4);
  dof_handler.distribute_dofs (fe);

  CompressedSimpleSparsityPattern csp (dof_handler.n_dofs(),
				       dof_handler.n_dofs());
  DoFTools::make_sparsity_pattern (dof_handler, csp);
  stiffness_matrix.reinit (csp);
  mass_matrix.reinit (csp);

  eigenfunctions
    .resize (parameters.get_integer ("Number of eigenvalues/eigenfunctions"));
  for (unsigned int i=0; i<eigenfunctions.size(); ++i)
    eigenfunctions[i].reinit (dof_handler.n_dofs());
}


                                 // @sect4{EigenvalueProblem::assemble_system}
template <int dim>
void EigenvalueProblem<dim>::assemble_system () 
{  
  QGauss<dim>  quadrature_formula(2);

  FEValues<dim> fe_values (fe, quadrature_formula, 
			   update_values   | update_gradients |
                           update_quadrature_points | update_JxW_values);

  const unsigned int dofs_per_cell = fe.dofs_per_cell;
  const unsigned int n_q_points    = quadrature_formula.size();

  FullMatrix<double> cell_stiffness_matrix (dofs_per_cell, dofs_per_cell);
  FullMatrix<double> cell_mass_matrix (dofs_per_cell, dofs_per_cell);

  std::vector<unsigned int> local_dof_indices (dofs_per_cell);

  FunctionParser<dim> potential;
  potential.initialize ((dim == 2 ?
			 "x,y" :
			 "x,y,z"),
			parameters.get ("Potential"),
			typename FunctionParser<dim>::ConstMap());
  std::vector<double> potential_values (n_q_points);

  ConstraintMatrix constraints;
  VectorTools::interpolate_boundary_values (dof_handler,
					    0,
					    ZeroFunction<dim>(),
					    constraints);
  constraints.close ();
  
  typename DoFHandler<dim>::active_cell_iterator
    cell = dof_handler.begin_active(),
    endc = dof_handler.end();
  for (; cell!=endc; ++cell)
    {
      fe_values.reinit (cell);
      cell_stiffness_matrix = 0;
      cell_mass_matrix = 0;

      potential.value_list (fe_values.get_quadrature_points(),
			    potential_values);
      
      for (unsigned int q_point=0; q_point<n_q_points; ++q_point)
	for (unsigned int i=0; i<dofs_per_cell; ++i)
	  for (unsigned int j=0; j<dofs_per_cell; ++j)
	    {
	      cell_stiffness_matrix(i,j)
		+= ((fe_values.shape_grad (i, q_point) *
		     fe_values.shape_grad (j, q_point)
		     +
		     potential_values[q_point] *
		     fe_values.shape_value (i, q_point) *
		     fe_values.shape_value (j, q_point)) *
		    fe_values.JxW (q_point));
	      cell_mass_matrix(i,j)
		+= (fe_values.shape_value (i, q_point) *
		    fe_values.shape_value (j, q_point) *
		    fe_values.JxW (q_point));
	    }
      
      cell->get_dof_indices (local_dof_indices);

      constraints.distribute_local_to_global (cell_stiffness_matrix,
					      local_dof_indices,
					      stiffness_matrix);
      constraints.distribute_local_to_global (cell_mass_matrix,
					      local_dof_indices,
					      mass_matrix);
    }
}


                                 // @sect4{EigenvalueProblem::solve}
template <int dim>
void EigenvalueProblem<dim>::solve () 
{
// do whatever is necessary here
}


                                 // @sect4{EigenvalueProblem::output_results}

template <int dim>
void EigenvalueProblem<dim>::output_results () const
{
  DataOut<dim> data_out;

  data_out.attach_dof_handler (dof_handler);

  for (unsigned int i=0; i<eigenfunctions.size(); ++i)
    data_out.add_data_vector (eigenfunctions[i],
			      std::string("solution") +
			      Utilities::int_to_string(i));

  data_out.build_patches ();

  std::ofstream output ("eigenvectors.vtk");
  data_out.write_vtk (output);
}



                                 // @sect4{EigenvalueProblem::run}

                                 // This is the function which has the
				 // top-level control over everything. It is
				 // the same as for the previous example.
template <int dim>
void EigenvalueProblem<dim>::run () 
{
  std::cout << "Solving problem in " << dim << " space dimensions." << std::endl;
  
  make_grid_and_dofs();
  assemble_system ();
  solve ();
  output_results ();
}


                                 // @sect3{The <code>main</code> function}

int main (int argc, char **argv) 
{
  try
    {
      PetscInitialize(&argc,&argv,0,0);
      deallog.depth_console (0);
      
      EigenvalueProblem<2> problem ("step-36.prm");
      problem.run ();
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
