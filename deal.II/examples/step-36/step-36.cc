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

				 // @sect3{Include files}

				 // As mentioned in the introduction, this
				 // program is essentially only a slightly
				 // revised version of step-4. As a
				 // consequence, most of the following include
				 // files are as used there, or at least as
				 // used already in previous tutorial
				 // programs:
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

				 // And then we need to actually import
				 // the interfaces for solvers that
				 // SLEPc provides:
#include <lac/slepc_solver.h>

				 // We also need some standard C++:
#include <fstream>
#include <iostream>

				 // Finally, as in previous programs, we
				 // import all the deal.II class and function
				 // names into the global namespace:
using namespace dealii;

				 // @sect3{The <code>EigenvalueProblem</code> class template}

				 // Following is the class declaration for the
				 // main class template. It looks pretty much
				 // exactly like what has already been shown
				 // in step-4:
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

				     // With these exceptions: For our
				     // eigenvalue problem, we need both a
				     // stiffness matrix for the left hand
				     // side as well as a mass matrix for the
				     // right hand side. We also need not just
				     // one solution function, but a whole set
				     // of those for the eigenfunctions we
				     // want to compute, along with the
				     // corresponding eigenvectors:
    PETScWrappers::SparseMatrix        stiffness_matrix, mass_matrix;
    std::vector<PETScWrappers::Vector> eigenfunctions;
    std::vector<double>                eigenvalues;   

				     // And then we need an object that will
				     // store several run-time parameters that
				     // we will specify in an input file:
    ParameterHandler parameters;

				     // Finally, we will have an object that
				     // contains "constraints" on our degrees
				     // of freedom. This could include hanging
				     // node constraints if we had adaptively
				     // refined meshes (which we don't have in
				     // the current program). Here, we will
				     // store the constraints for boundary
				     // nodes $U_i=0$.
    ConstraintMatrix constraints;
};

				 // @sect3{Implementation of the <code>EigenvalueProblem</code> class}

				 // @sect4{EigenvalueProblem::EigenvalueProblem}

				 // First up, the constructor. The main, new
				 // part is handling the run-time input
				 // parameters. We need to declare their
				 // existence first, and then read their
				 // values from the input file whose name is
				 // specified as an argument to this function:
template <int dim>
EigenvalueProblem<dim>::EigenvalueProblem (const std::string &prm_file)
		:
		fe (1),
		dof_handler (triangulation)
{
  parameters.declare_entry ("Global mesh refinement steps", "5",
			    Patterns::Integer (0, 20),
			    "The number number of times the 1-cell coarse mesh should "
			    "be refined globally for our computations.");
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

				 // The next function creates a mesh on the
				 // domain $[-1,1]^d$, refines it as many
				 // times as the input file calls for, and
				 // then attaches a DoFHandler to it and
				 // initializes the matrices and vectors to
				 // their correct sizes. We also build the
				 // constraints that correspond to the
				 // boundary values $u|_{\partial\Omega}=0$.
				 //
				 // For the matrices, we use the PETSc
				 // wrappers. These have the ability to
				 // allocate memory as necessary as non-zero
				 // entries are added. This seems inefficient:
				 // we could as well first compute the
				 // sparsity pattern, initialize the matrices
				 // with it, and as we then insert entries we
				 // can be sure that we do not need to
				 // re-allocate memory and free the one used
				 // previously. One way to do that would be to
				 // use code like this:
				 // @code
				 //   CompressedSimpleSparsityPattern
				 //      csp (dof_handler.n_dofs(),
				 //           dof_handler.n_dofs());
				 //   DoFTools::make_sparsity_pattern (dof_handler, csp);
				 //   csp.compress ();
				 //   stiffness_matrix.reinit (csp);
				 //   mass_matrix.reinit (csp);
				 // @code
				 // instead of the two <code>reinit()</code>
				 // calls for the stiffness and mass matrices.
				 //
				 // This doesn't quite work,
				 // unfortunately. The code above may lead to
				 // a few entries in the non-zero pattern to
				 // which we only ever write zero entries;
				 // most notably, this holds true for
				 // off-diagonal entries for those rows and
				 // columns that belong to boundary
				 // nodes. This shouldn't be a problem, but
				 // for whatever reason, PETSc's ILU
				 // preconditioner, which we use to solve
				 // linear systems in the eigenvalue solver,
				 // doesn't like these extra entries and
				 // aborts with an error message.
				 //
				 // Absent any obvious way to avoid this, we
				 // simply settle for the second best option,
				 // which is have PETSc allocate memory as
				 // necessary. That said, since this is not a
				 // time critical part, this whole affair is
				 // of no further importance.
template <int dim>
void EigenvalueProblem<dim>::make_grid_and_dofs ()
{
  GridGenerator::hyper_cube (triangulation, -1, 1);
  triangulation.refine_global (parameters.get_integer ("Global mesh refinement steps"));
  dof_handler.distribute_dofs (fe);

  DoFTools::make_zero_boundary_constraints (dof_handler, constraints);
  constraints.close ();
  
  stiffness_matrix.reinit (dof_handler.n_dofs(),
			   dof_handler.n_dofs(),
  			   dof_handler.max_couplings_between_dofs());
  mass_matrix.reinit (dof_handler.n_dofs(),
		      dof_handler.n_dofs(),
  		      dof_handler.max_couplings_between_dofs());

				   // The next step is to take care of the
				   // eigenspectrum. In this case, the outputs
				   // are eigenfunctions and eigenvalues, so
				   // we set the size of the list of
				   // eigenfunctions and eigenvalues to be as
				   // large as asked for in the input file:
  eigenfunctions
    .resize (parameters.get_integer ("Number of eigenvalues/eigenfunctions"));
  for (unsigned int i=0; i<eigenfunctions.size (); ++i)
    eigenfunctions[i].reinit (dof_handler.n_dofs ());

  eigenvalues.resize (eigenfunctions.size ());
}


				 // @sect4{EigenvalueProblem::assemble_system}

				 // Here, we assemble the global stiffness and
				 // mass matrices from local contributions
				 // $A^K_{ij} = \int_K \nabla\varphi_i(\mathbf
				 // x) \cdot \nabla\varphi_j(\mathbf x) +
				 // V(\mathbf x)\varphi_i(\mathbf
				 // x)\varphi_j(\mathbf x)$. The function
				 // should be immediately familiar if you've
				 // seen previous tutorial programs. The only
				 // thing new would be setting up an object
				 // that described the potential $V(\mathbf
				 // x)$ using the expression that we got from
				 // the input file. We then need to evaluate
				 // this object at the quadrature points on
				 // each cell. If you've seen how to evaluate
				 // function objects (see, for example the
				 // coefficient in step-5), the code here will
				 // also look rather familiar.
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
  
  std::vector<double> potential_values (n_q_points);
  
  
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
      
      for (unsigned int q_point=0; q_point<n_q_points; ++q_point)
	for (unsigned int i=0; i<dofs_per_cell; ++i)
	  for (unsigned int j=0; j<dofs_per_cell; ++j)
	    {
	      cell_stiffness_matrix (i, j)
		+= (fe_values.shape_grad (i, q_point) *
		    fe_values.shape_grad (j, q_point) 
		    + 
		    potential_values[q_point] *
		    fe_values.shape_value (i, q_point) *
		    fe_values.shape_value (j, q_point)
		) * fe_values.JxW (q_point);
	      
	      cell_mass_matrix (i, j)
		+= (fe_values.shape_value (i, q_point) *
		    fe_values.shape_value (j, q_point) 
		) * fe_values.JxW (q_point);
	    }

				       // Now that we have the local matrix
				       // contributions, we transfer them into
				       // the global objects and take care of
				       // zero boundary constraints:
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

				   // At the end of the function, we tell
				   // PETSc that the matrices have now been
				   // fully assembled and that the sparse
				   // matrix representation can now be
				   // compressed as no more entries will be
				   // added:
  stiffness_matrix.compress();
  mass_matrix.compress();
}


				 // @sect4{EigenvalueProblem::solve}
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

				   // Now rescale eigenfunctions so that they
				   // have $\|\phi_i(\mathbf
				   // x)\|_{L^\infty(\Omega)}=1$ instead of
				   // $\|\Phi\|_{l_2}=1$:
  for (unsigned int i=0; i<eigenfunctions.size(); ++i)
    eigenfunctions[i] /= eigenfunctions[i].linfty_norm ();
}


				 // @sect4{EigenvalueProblem::output_results}
template <int dim>
void EigenvalueProblem<dim>::output_results () const
{
  DataOut<dim> data_out;

  data_out.attach_dof_handler (dof_handler);

  for (unsigned int i=0; i<eigenfunctions.size(); ++i)
    data_out.add_data_vector (eigenfunctions[i],
			      std::string("eigenfunction_") +
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
}


                                 // @sect4{EigenvalueProblem::run}

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

  for (unsigned int i=0; i<eigenvalues.size(); ++i)
    std::cout << std::endl 
	      << "      eigenvalue " << i 
	      << " : " << eigenvalues[i];
}


                                 // @sect3{The <code>main</code> function}
int main (int argc, char **argv) 
{
  try
    {

				       // Here is another difference from
				       // other steps: We initialize the SLEPc
				       // work space which inherently
				       // initializes the PETSc work space,
				       // run the whole program, ...
      SlepcInitialize (&argc,&argv,0,0);

      {
	deallog.depth_console (0);
	
	EigenvalueProblem<2> problem ("step-36.prm");
	problem.run ();
      }

				       // ...and then unitialize the SLEPc
				       // work space when the job is done:
      SlepcFinalize ();
    }

				   // All the while, we are watching out if
				   // any exceptions should have been
				   // generated. If that is so, we panic...
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
  
				   // ...or show that we are happy:
  std::cout << std::endl 
	    << "Done." 
	    << std::endl;

  return 0;
}
