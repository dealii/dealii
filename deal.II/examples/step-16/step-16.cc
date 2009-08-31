/* $Id$ */
/* Author: Guido Kanschat, University of Heidelberg, 2003 */

/*    $Id$       */
/*                                                                */
/*    Copyright (C) 2003, 2004, 2005, 2006, 2007, 2008, 2009 by the deal.II authors */
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
#include <lac/vector.h>
#include <lac/full_matrix.h>
#include <lac/sparse_matrix.h>
#include <lac/solver_cg.h>
#include <lac/precondition.h>
#include <grid/tria.h>
#include <grid/tria_accessor.h>
#include <grid/tria_iterator.h>
#include <grid/tria_boundary_lib.h>
#include <grid/grid_generator.h>
#include <dofs/dof_handler.h>
#include <dofs/dof_accessor.h>
#include <dofs/dof_tools.h>
#include <fe/fe_q.h>
#include <fe/fe_values.h>
#include <numerics/vectors.h>
#include <numerics/matrices.h>
#include <numerics/data_out.h>
				 // These are the new include files
				 // required for multi-level methods.
				 // First, the file defining the
				 // multigrid method itself.
#include <multigrid/multigrid.h>
				 // The DoFHandler is replaced by an
				 // MGDoFHandler which is defined
				 // here.
#include <multigrid/mg_dof_handler.h>
#include <multigrid/mg_dof_accessor.h>

				 // Then, we need some pre-made
				 // transfer routines between grids.
#include <multigrid/mg_transfer.h>
                                 // And a file in which equivalents to the
                                 // DoFTools class are declared:
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

				 // This class is based on the same
				 // class in step-5. Remark that we
				 // replaced the DoFHandler by
				 // MGDoFHandler. since this inherits
				 // from DoFHandler, the new object
				 // incorporates the old functionality
				 // plus the new functions for degrees
				 // of freedom on different
				 // levels. Furthermore, we added
				 // MGLevelObjects for sparsity
				 // patterns and matrices.
template <int dim>
class LaplaceProblem 
{
  public:
    LaplaceProblem ();
    void run ();
    
  private:
    void setup_system ();
    void assemble_system ();
				     // We add this function for
				     // assembling the multilevel
				     // matrices.
    void assemble_multigrid ();
    void solve ();
    void output_results (const unsigned int cycle) const;

    Triangulation<dim>   triangulation;
    FE_Q<dim>            fe;
    MGDoFHandler<dim>    mg_dof_handler;

    SparsityPattern      sparsity_pattern;
    SparseMatrix<double> system_matrix;

				     // Here are the new objects for
				     // handling level matrices: sparsity
				     // patterns and matrices. We use number
				     // type float to save memory. It's only
				     // a preconditioner!
    MGLevelObject<SparsityPattern>      mg_sparsity;
    MGLevelObject<SparseMatrix<float> > mg_matrices;
    
    Vector<double>       solution;
    Vector<double>       system_rhs;
};



				 // This function is as before.
template <int dim>
LaplaceProblem<dim>::LaplaceProblem () :
                fe (1),
		mg_dof_handler (triangulation)
{}



				 // This is the function of step-5
				 // augmented by the setup of the
				 // multi-level matrices in the end.
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
  DoFTools::make_sparsity_pattern (static_cast<const DoFHandler<dim>&>(mg_dof_handler),
                                   sparsity_pattern);
  sparsity_pattern.compress();

  system_matrix.reinit (sparsity_pattern);

  solution.reinit (mg_dof_handler.n_dofs());
  system_rhs.reinit (mg_dof_handler.n_dofs());

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
  mg_matrices.resize(0, nlevels-1);
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
      mg_sparsity[level].compress();
      mg_matrices[level].reinit(mg_sparsity[level]);
    }
}

				 // This is the standard assemble function
				 // for the Poisson equation you have seen a
				 // lot of times before.
template <int dim>
void LaplaceProblem<dim>::assemble_system () 
{  
  QGauss<dim>  quadrature_formula(2);

  FEValues<dim> fe_values (fe, quadrature_formula, 
			   update_values   | update_gradients |
                           update_quadrature_points | update_JxW_values);

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

				       // As before, we want the FEValues
				       // object to compute the quantities
				       // which we told him to compute in
				       // the constructor using the update
				       // flags. Then, we loop over all
				       // quadrature points and the local
				       // matrix rows and columns for
				       // computing the element
				       // contribution. This is the same as
				       // in step-4. For the right hand
				       // side, we use a constant value of
				       // 1.
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

				   // The Dirichlet boundary
				   // conditions on the finest level
				   // are handled as usual.
  std::map<unsigned int,double> boundary_values;  

  VectorTools::interpolate_boundary_values (mg_dof_handler,
					    0,
					    ZeroFunction<dim>(),
					    boundary_values);

  MatrixTools::apply_boundary_values (boundary_values,
				      system_matrix,
				      solution,
				      system_rhs);
}


				 // Here is another assemble
				 // function. The integration core is
				 // the same as above. Only the loop
				 // goes over all existing cells now
				 // and the results must be entered
				 // into the correct matrix.

				 // Since we only do multi-level
				 // preconditioning, no right-hand
				 // side is assembled here.
template <int dim>
void LaplaceProblem<dim>::assemble_multigrid () 
{  
  QGauss<dim>  quadrature_formula(2);

  FEValues<dim> fe_values (fe, quadrature_formula, 
			   update_values   | update_gradients |
                           update_quadrature_points | update_JxW_values);

  const unsigned int   dofs_per_cell = fe.dofs_per_cell;
  const unsigned int   n_q_points    = quadrature_formula.size();

  FullMatrix<double>   cell_matrix (dofs_per_cell, dofs_per_cell);

  std::vector<unsigned int> local_dof_indices (dofs_per_cell);

  typename MGDoFHandler<dim>::cell_iterator cell = mg_dof_handler.begin(),
					    endc = mg_dof_handler.end();
  for (; cell!=endc; ++cell)
    {
				       // Remember the level of the
				       // current cell.
      const unsigned int level = cell->level();
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
      for (unsigned int i=0; i<dofs_per_cell; ++i)
	{
	  for (unsigned int j=0; j<dofs_per_cell; ++j)
					     // And now add everything
					     // to the matrix on the
					     // right level.
	    mg_matrices[level].add (local_dof_indices[i],
                                    local_dof_indices[j],
                                    cell_matrix(i,j));
	}
    }

				   // Now we have to eliminate the
				   // boundary nodes on all
				   // levels. This is done exactly in
				   // the same fashion as on the
				   // finest level. Therefore, we need
				   // a function map first. On the
				   // other hand, since we use
				   // multigrid only as a
				   // preconditioner, we always use
				   // homogeneous boundary conditions
				   // and the ZeroFunction will be
				   // sufficient in all cases; using a
				   // different function here would
				   // not hurt on the other hand,
				   // since the values are ignored.

  typename FunctionMap<dim>::type      dirichlet_boundary;
  ZeroFunction<dim>                    homogeneous_dirichlet_bc (1);
  dirichlet_boundary[0] = &homogeneous_dirichlet_bc;
  
				   // Next we generate the set of
				   // dof indices to be eliminated on
				   // each level.
  std::vector<std::set<unsigned int> > boundary_indices(triangulation.n_levels());
  MGTools::make_boundary_list (mg_dof_handler,
			       dirichlet_boundary,
			       boundary_indices);

				   // And finally we eliminate these
				   // degrees of freedom from every
				   // matrix in the multilevel hierarchy.
 const unsigned int nlevels = triangulation.n_levels();
 for (unsigned int level=0;level<nlevels;++level)
   {
     MGTools::apply_boundary_values(boundary_indices[level],
				    mg_matrices[level],
				    true);
   }
}



				 // The solution process again looks
				 // mostly like in the previous
				 // examples. However, we will now use
				 // a preconditioned conjugate
				 // gradient algorithm. It is not very
				 // difficult to make this change:
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
  MGTransferPrebuilt<Vector<double> > mg_transfer;
  mg_transfer.build_matrices(mg_dof_handler);

				   // Next, we need a coarse grid
				   // solver. Since our coarse grid is
				   // VERY coarse, we decide for a
				   // direct solver, even if its
				   // implementation is not very
				   // clever.
  FullMatrix<float> coarse_matrix;
  coarse_matrix.copy_from (mg_matrices[0]);
  MGCoarseGridHouseholder<float, Vector<double> > mg_coarse;
  mg_coarse.initialize(coarse_matrix);

				   // The final ingredient for the
				   // multilevel preconditioner is the
				   // smoother. It is very customary
				   // to use a relaxation method
				   // here. Names are getting quite
				   // long here, so we help with
				   // typedefs.
  typedef PreconditionSOR<SparseMatrix<float> > RELAXATION;
  MGSmootherRelaxation<SparseMatrix<float>, RELAXATION, Vector<double> >
    mg_smoother(vector_memory);

				   // Initialize the smoother with our
				   // level matrices and the (required)
				   // additional data for the
				   // relaxation method with default
				   // values.
  RELAXATION::AdditionalData smoother_data;
  mg_smoother.initialize(mg_matrices, smoother_data);
  
				   // Do two smoothing steps per level
  mg_smoother.set_steps(2);
				   // Since the SOR method is not
				   // symmetric, but we use conjugate
				   // gradient iteration below, here
				   // is a trick to make the
				   // multilevel preconditioner a
				   // symmetric operator even for
				   // nonsymmetric smoothers.
  mg_smoother.set_symmetric(true);

				   // We must wrap our matrices in an
				   // object having the required
				   // multiplication functions.
  MGMatrix<SparseMatrix<float>, Vector<double> >
    mg_matrix(&mg_matrices);
				   // Now, we are ready to set up the
				   // V-cycle operator and the
				   // multilevel preconditioner.
  Multigrid<Vector<double> > mg(mg_dof_handler,
				mg_matrix,
				mg_coarse,
				mg_transfer,
				mg_smoother,
				mg_smoother);
  PreconditionMG<dim, Vector<double>,
    MGTransferPrebuilt<Vector<double> > >
    preconditioner(mg_dof_handler, mg, mg_transfer);
  
				   // Finally, create the solver
				   // object and solve the system
  SolverControl           solver_control (1000, 1e-12);
  SolverCG<>              cg (solver_control);

  
  cg.solve (system_matrix, solution, system_rhs,
	    preconditioner);
}



				 // Here is the data output, which is
				 // a simplified version of step-5. We
				 // do a standard gnuplot output for
				 // each grid produced in the
				 // refinement process.
template <int dim>
void LaplaceProblem<dim>::output_results (const unsigned int cycle) const
{
				   // Construct and initialize a DataOut object
  DataOut<dim> data_out;

  data_out.attach_dof_handler (mg_dof_handler);
  data_out.add_data_vector (solution, "solution");
  data_out.build_patches ();

				   // The following block generates
				   // the file name and opens the
				   // file:
  std::ostringstream filename;
  filename << "solution-"
	   << cycle
	   << ".gnuplot";

  std::ofstream output (filename.str().c_str());
  data_out.write_gnuplot (output);
}



template <int dim>
void LaplaceProblem<dim>::run () 
{
  for (unsigned int cycle=0; cycle<6; ++cycle)
    {
      deallog << "Cycle " << cycle << std::endl;

      if (cycle == 0)
	{
					   // Generate a simple hyperball grid.
	  GridGenerator::hyper_ball(triangulation);
	  static const HyperBallBoundary<dim> boundary;
	  triangulation.set_boundary (0, boundary);
	}
      triangulation.refine_global (1);
      setup_system ();
      assemble_system ();
      assemble_multigrid ();
      solve ();
      output_results (cycle);
    };
}

    

				 // The main function looks mostly
				 // like the one in the previous
				 // example, so we won't comment on it
				 // further.
int main () 
{
  LaplaceProblem<2> laplace_problem_2d;
  laplace_problem_2d.run ();
  
  return 0;
}
