/* $Id$ */
/* Author: Guido Kanschat, University of Heidelberg, 2003  */
/*         Baerbel Janssen, University of Heidelberg, 2010 */
/*         Wolfgang Bangerth, Texas A&M University, 2010   */

/*    $Id$       */
/*                                                                */
/*    Copyright (C) 2003, 2004, 2006, 2007, 2008, 2009, 2010 by the deal.II authors                   */
/*                                                                */
/*    This file is subject to QPL and may not be  distributed     */
/*    without copyright and license information. Please refer     */
/*    to the file deal.II/doc/license.html for the  text  and     */
/*    further information on this license.                        */

				 // Again, the first few include files
				 // are already known, so we won't
				 // comment on them:
#include <base/quadrature_lib.h>
#include <base/function.h>
#include <base/logstream.h>
#include <base/utilities.h>

#include <lac/constraint_matrix.h>
#include <lac/vector.h>
#include <lac/full_matrix.h>
#include <lac/sparse_matrix.h>
#include <lac/solver_cg.h>
#include <lac/precondition.h>

#include <grid/tria.h>
#include <grid/tria_accessor.h>
#include <grid/tria_iterator.h>
#include <grid/grid_generator.h>
#include <grid/grid_refinement.h>
#include <grid/tria_boundary_lib.h>

#include <dofs/dof_accessor.h>
#include <dofs/dof_tools.h>

#include <fe/fe_q.h>
#include <fe/fe_values.h>

#include <numerics/vectors.h>
#include <numerics/data_out.h>
#include <numerics/error_estimator.h>

//These are the same include files
//as in step-16 necessary for the
//multi-level methods
#include <multigrid/multigrid.h>
#include <multigrid/mg_dof_handler.h>
#include <multigrid/mg_dof_accessor.h>
#include <multigrid/mg_transfer.h>
#include <multigrid/mg_tools.h>
#include <multigrid/mg_coarse.h>
#include <multigrid/mg_smoother.h>
#include <multigrid/mg_matrix.h>

				 // This is C++:
#include <fstream>
#include <sstream>

				 // The last step is as in all
				 // previous programs:
using namespace dealii;


//This class is basically the same
//class as in step-16. The only
//difference is that here we solve Laplace's
//problem on an adaptively refined grid.
template <int dim>
class LaplaceProblem
{
  public:
    LaplaceProblem (const unsigned int deg);
    void run ();

  private:
    void setup_system ();
    void assemble_system ();
    void assemble_multigrid ();
    void solve ();
    void refine_grid ();
    void output_results (const unsigned int cycle) const;

    Triangulation<dim>   triangulation;
    FE_Q<dim>            fe;
    MGDoFHandler<dim>    mg_dof_handler;

    SparsityPattern      sparsity_pattern;
    SparseMatrix<double> system_matrix;

				     //This object holds the information f
				     //or the hanging nodes.
    ConstraintMatrix     constraints;

    MGLevelObject<SparsityPattern> mg_sparsity;
    MGLevelObject<SparseMatrix<double> > mg_matrices;

				     /* The matrices at the interface
				      * between two refinement levels,
				      * coupling coarse to fine.*/
    MGLevelObject<SparseMatrix<double> > mg_interface_matrices_up;

    Vector<double>       solution;
    Vector<double>       system_rhs;

    const unsigned int degree;
};


template <int dim>
LaplaceProblem<dim>::LaplaceProblem (const unsigned int deg)
		:
		triangulation (Triangulation<dim>::limit_level_difference_at_vertices),
		fe (deg),
		mg_dof_handler (triangulation),
		degree(deg)
{}


template <int dim>
void LaplaceProblem<dim>::setup_system ()
{
  mg_dof_handler.distribute_dofs (fe);

				   // Here we output not only the
				   // degrees of freedom on the finest
				   // level, but also in the
				   // multilevel structure
  deallog << "Number of degrees of freedom: "
	  << mg_dof_handler.n_dofs();

  for (unsigned int l=0;l<triangulation.n_levels();++l)
    deallog << "   " << 'L' << l << ": "
	    << mg_dof_handler.n_dofs(l);
  deallog  << std::endl;

  sparsity_pattern.reinit (mg_dof_handler.n_dofs(),
			   mg_dof_handler.n_dofs(),
			   mg_dof_handler.max_couplings_between_dofs());
  DoFTools::make_sparsity_pattern (
    static_cast<const DoFHandler<dim>&>(mg_dof_handler),
    sparsity_pattern);

  solution.reinit (mg_dof_handler.n_dofs());
  system_rhs.reinit (mg_dof_handler.n_dofs());

  constraints.clear ();
  DoFTools::make_hanging_node_constraints (mg_dof_handler, constraints);
  VectorTools::interpolate_boundary_values (mg_dof_handler,
					    0,
					    ZeroFunction<dim>(),
					    constraints);
  constraints.close ();
  constraints.condense (sparsity_pattern);
  sparsity_pattern.compress();
  system_matrix.reinit (sparsity_pattern);

				   // The multi-level objects are
				   // resized to hold matrices for
				   // every level. The coarse level is
				   // zero (this is mandatory right
				   // now but may change in a future
				   // revision). Remark, that the
				   // finest level is nlevels-1.
				   // We first have to resize the
				   // container holding the
				   // SparseMatrix classes, since they
				   // have to release their
				   // SparsityPattern before it can be
				   // destroyed.
  const unsigned int nlevels = triangulation.n_levels();

  mg_interface_matrices_up.resize(0, nlevels-1);
  mg_interface_matrices_up.clear ();
  mg_matrices.resize(0, nlevels-1);
  mg_matrices.clear ();
  mg_sparsity.resize(0, nlevels-1);

				   // Now, we have to build a matrix
				   // on each level. Technically, we
				   // could use the matrix initialized
				   // above on the finest
				   // level. Beware that this is not
				   // true anymore with local
				   // refinement!
  for (unsigned int level=0;level<nlevels;++level)
    {
      mg_sparsity[level].reinit (mg_dof_handler.n_dofs(level),
				 mg_dof_handler.n_dofs(level),
				 mg_dof_handler.max_couplings_between_dofs());
      MGTools::make_sparsity_pattern (mg_dof_handler, mg_sparsity[level], level);
      CompressedSparsityPattern ci_sparsity;
      if(level>0)
	{
	  ci_sparsity.reinit(mg_dof_handler.n_dofs(level),
			     mg_dof_handler.n_dofs(level));
	  MGTools::make_sparsity_pattern(mg_dof_handler, ci_sparsity, level);
	}
    }

//And the same for the mg matrices
//for the interface. Note that there
//is no such interface on the coarsest level
  for(unsigned int level=0; level<nlevels; ++level)
    {
      mg_sparsity[level].compress();
      mg_matrices[level].reinit(mg_sparsity[level]);
      mg_interface_matrices_up[level].reinit(mg_sparsity[level]);
    }
}

				 // This is the standard assemble function
				 // for the Poisson equation you have seen a
				 // lot of times before.
template <int dim>
void LaplaceProblem<dim>::assemble_system ()
{
  QGauss<dim>  quadrature_formula(1+degree);

  FEValues<dim> fe_values (fe, quadrature_formula,
			   update_values   | update_gradients |
			   update_quadrature_points | update_JxW_values);

  const unsigned int   dofs_per_cell = fe.dofs_per_cell;
  const unsigned int   n_q_points    = quadrature_formula.size();

  FullMatrix<double>   cell_matrix (dofs_per_cell, dofs_per_cell);
  Vector<double>       cell_rhs (dofs_per_cell);

  std::vector<unsigned int> local_dof_indices (dofs_per_cell);

  typename DoFHandler<dim>::active_cell_iterator
    cell = mg_dof_handler.begin_active(),
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
	      cell_matrix(i,j) += (fe_values.shape_grad(i,q_point)
				   * fe_values.shape_grad(j,q_point))
				  * fe_values.JxW(q_point);

	    cell_rhs(i) += (fe_values.shape_value(i,q_point)
			    * 1.0 * fe_values.JxW(q_point));
	  }

      cell->get_dof_indices (local_dof_indices);
      constraints.distribute_local_to_global (cell_matrix, cell_rhs,
					      local_dof_indices,
					      system_matrix, system_rhs);
    }
}


				 // Here is another assemble
				 // function. The integration core is
				 // the same as above. Only the loop
				 // goes over all existing cells now
				 // and the results must be entered
				 // into the correct matrix. Here comes
                                 // the difference to global refinement
                                 // into play. We have to fill the interface
                                 // matrices correctly.

				 // Since we only do multi-level
				 // preconditioning, no right-hand
				 // side is assembled here.
template <int dim>
void LaplaceProblem<dim>::assemble_multigrid ()
{
  QGauss<dim>  quadrature_formula(1+degree);

  FEValues<dim> fe_values (fe, quadrature_formula,
			   update_values   | update_gradients |
			   update_quadrature_points | update_JxW_values);

  const unsigned int   dofs_per_cell   = fe.dofs_per_cell;
  const unsigned int   n_q_points      = quadrature_formula.size();

  FullMatrix<double>   cell_matrix (dofs_per_cell, dofs_per_cell);

  std::vector<unsigned int> local_dof_indices (dofs_per_cell);

  std::vector<std::vector<bool> > interface_dofs;
  std::vector<std::vector<bool> > boundary_interface_dofs;
  for (unsigned int level = 0; level<triangulation.n_levels(); ++level)
    {
      std::vector<bool> tmp (mg_dof_handler.n_dofs(level));
      interface_dofs.push_back (tmp);
      boundary_interface_dofs.push_back (tmp);
    }
  MGTools::extract_inner_interface_dofs (mg_dof_handler,
					 interface_dofs,
					 boundary_interface_dofs);

  typename FunctionMap<dim>::type      dirichlet_boundary;
  ZeroFunction<dim>                    homogeneous_dirichlet_bc (1);
  dirichlet_boundary[0] = &homogeneous_dirichlet_bc;

  std::vector<std::set<unsigned int> > boundary_indices(triangulation.n_levels());
  MGTools::make_boundary_list (mg_dof_handler, dirichlet_boundary,
			       boundary_indices);

  std::vector<ConstraintMatrix> boundary_constraints (triangulation.n_levels());
  std::vector<ConstraintMatrix> boundary_interface_constraints (triangulation.n_levels());
  for (unsigned int level=0; level<triangulation.n_levels(); ++level)
    {
      boundary_constraints[level].add_lines (interface_dofs[level]);
      boundary_constraints[level].add_lines (boundary_indices[level]);
      boundary_constraints[level].close ();

      boundary_interface_constraints[level]
	.add_lines (boundary_interface_dofs[level]);
      boundary_interface_constraints[level].close ();
    }

  typename MGDoFHandler<dim>::cell_iterator cell = mg_dof_handler.begin(),
					    endc = mg_dof_handler.end();

  for (; cell!=endc; ++cell)
    {
      cell_matrix = 0;

				       // Compute the values specified
				       // by update flags above.
      fe_values.reinit (cell);

				       // This is exactly the
				       // integration loop of the cell
				       // matrix above.
      for (unsigned int q_point=0; q_point<n_q_points; ++q_point)
	for (unsigned int i=0; i<dofs_per_cell; ++i)
	  {
	    for (unsigned int j=0; j<dofs_per_cell; ++j)
	      cell_matrix(i,j) += (fe_values.shape_grad(i,q_point)
				   * fe_values.shape_grad(j,q_point))
				  * fe_values.JxW(q_point);
	  }

				       // Oops! This is a tiny
				       // difference easily
				       // forgotten. The indices we
				       // want here are the ones for
				       // that special level, not for
				       // the global
				       // matrix. Therefore, a little
				       // 'mg' entered into the
				       // function call.
      cell->get_mg_dof_indices (local_dof_indices);

      const unsigned int level = cell->level();
      boundary_constraints[level]
	.distribute_local_to_global (cell_matrix,
				     local_dof_indices,
				     mg_matrices[level]);

      for (unsigned int i=0; i<dofs_per_cell; ++i)
	for (unsigned int j=0; j<dofs_per_cell; ++j)
	  if( !(interface_dofs[level][local_dof_indices[i]]==true &&
		interface_dofs[level][local_dof_indices[j]]==false))
	    cell_matrix(i,j) = 0;

      boundary_interface_constraints[level]
	.distribute_local_to_global (cell_matrix,
				     local_dof_indices,
				     mg_interface_matrices_up[level]);
    }
}

template <int dim>
void LaplaceProblem<dim>::solve ()
{
				   // Create a memory handler for
				   // regular vectors. Note, that
				   // GrowingVectorMemory is more time
				   // efficient than the
				   // PrimitiveVectorMemory class.
  GrowingVectorMemory<>   vector_memory;

				   // Now, create an object handling
				   // the transfer of functions
				   // between different grid
				   // levels.
  MGTransferPrebuilt<Vector<double> > mg_transfer(constraints);
  mg_transfer.build_matrices(mg_dof_handler);

				   // Next, we need a coarse grid
				   // solver. Since our coarse grid is
				   // VERY coarse, we decide for a
				   // direct solver, even if its
				   // implementation is not very
				   // clever.
  FullMatrix<double> coarse_matrix;
  coarse_matrix.copy_from (mg_matrices[0]);
  MGCoarseGridHouseholder<double, Vector<double> > mg_coarse;
  mg_coarse.initialize(coarse_matrix);

				   // The final ingredient for the
				   // multilevel preconditioner is the
				   // smoother. It is very customary
				   // to use a relaxation method
				   // here. Names are getting quite
				   // long here, so we help with
				   // typedefs.
  typedef PreconditionSOR<SparseMatrix<double> > RELAXATION;
//  typedef PreconditionJacobi<SparseMatrix<double> > RELAXATION;
//  typedef SparseILU<double> RELAXATION;
  MGSmootherRelaxation<SparseMatrix<double>, RELAXATION, Vector<double> >
    mg_smoother(vector_memory);

				   // Initialize the smoother with our
				   // level matrices and the required,
				   // additional data for the
				   // relaxaton method with default
				   // values.
  RELAXATION::AdditionalData smoother_data;//(0, 9,false);
  mg_smoother.initialize(mg_matrices, smoother_data);

				   // Do two smoothing steps per level
  mg_smoother.set_steps(1);
				   // Since the SOR method is not
				   // symmetric, but we use conjugate
				   // gradient iteration below, here
				   // is a trick to make the
				   // multilevel preconditioner a
				   // symmetric operator even for
				   // nonsymmetric smoothers.
  mg_smoother.set_symmetric(true);
  mg_smoother.set_variable(true);


				   // We must wrap our matrices in an
				   // object having the required
				   // multiplication functions.
  MGMatrix<SparseMatrix<double>, Vector<double> >
    mg_matrix(&mg_matrices);
				   //do the same for the interface matrices
  MGMatrix<SparseMatrix<double>, Vector<double> >
    mg_interface_up(&mg_interface_matrices_up);
  MGMatrix<SparseMatrix<double>, Vector<double> >
    mg_interface_down(&mg_interface_matrices_up);
				   // Now, we are ready to set up the
				   // V-cycle operator and the
				   // multilevel preconditioner.
  Multigrid<Vector<double> > mg(mg_dof_handler,
				mg_matrix,
				mg_coarse,
				mg_transfer,
				mg_smoother,
				mg_smoother);
  mg.set_edge_matrices(mg_interface_down, mg_interface_up);

  PreconditionMG<dim, Vector<double>,
    MGTransferPrebuilt<Vector<double> > >
  preconditioner(mg_dof_handler, mg, mg_transfer);

				   // Finally, create the solver
				   // object and solve the system
  ReductionControl solver_control (100, 1.e-20, 1.e-10, true, true);
  SolverCG<>    cg (solver_control);

  solution = 0;

  cg.solve (system_matrix, solution, system_rhs,
	    preconditioner);
  constraints.distribute (solution);

  std::cout << "   " << solver_control.last_step()
	    << " CG iterations needed to obtain convergence."
	    << std::endl;
}

template <int dim>
void LaplaceProblem<dim>::refine_grid ()
{
  Vector<float> estimated_error_per_cell (triangulation.n_active_cells());

  KellyErrorEstimator<dim>::estimate (static_cast<DoFHandler<dim>&>(mg_dof_handler),
				      QGauss<dim-1>(3),
				      typename FunctionMap<dim>::type(),
				      solution,
				      estimated_error_per_cell);
  GridRefinement::refine_and_coarsen_fixed_number (triangulation,
						   estimated_error_per_cell,
						   0.3, 0.03);
  triangulation.execute_coarsening_and_refinement ();
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
	   << ".vtk";

  std::ofstream output (filename.str().c_str());
  data_out.write_vtk (output);
}

template <int dim>
void LaplaceProblem<dim>::run ()
{
  for (unsigned int cycle=0; cycle<8; ++cycle)
    {
      std::cout << "Cycle " << cycle << ':' << std::endl;

      if (cycle == 0)
	{
	  GridGenerator::hyper_ball (triangulation);

	  static const HyperBallBoundary<dim> boundary;
	  triangulation.set_boundary (0, boundary);

	  triangulation.refine_global (1);
	}
      else
	refine_grid ();


      std::cout << "   Number of active cells:       "
		<< triangulation.n_active_cells()
		<< std::endl;

      setup_system ();

      std::cout << "   Number of degrees of freedom: "
		<< mg_dof_handler.n_dofs()
		<< std::endl;

      assemble_system ();
      assemble_multigrid ();

      solve ();
      output_results (cycle);
    }
}



int main ()
{
  try
    {
      deallog.depth_console (0);

      LaplaceProblem<2> laplace_problem(1);
      laplace_problem.run ();
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
