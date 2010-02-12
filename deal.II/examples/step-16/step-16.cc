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

				 // As discussed in the introduction, most of
				 // this program is copied almost verbatim
				 // from step-6, which itself is only a slight
				 // modification of step-5. Consequently, a
				 // significant part of this program is not
				 // new if you've read all the material up to
				 // step-6, and we won't comment on that part
				 // of the functionality that is
				 // unchanged. Rather, we will focus on those
				 // aspects of the program that have to do
				 // with the multigrid functionality which
				 // forms the new aspect of this tutorial
				 // program.

                                 // @sect3{Include files}

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

				 // These, now, are the include necessary for
				 // the multi-level methods. The first two
				 // declare classes that allow us to enumerate
				 // degrees of freedom not only on the finest
				 // mesh level, but also on intermediate
				 // levels (that's what the MGDoFHandler class
				 // does) as well as allow to access this
				 // information (iterators and accessors over
				 // these cells).
				 //
				 // The rest of the include files deals with
				 // the mechanics of multigrid as a linear
				 // operator (solver or preconditioner).
#include <multigrid/mg_dof_handler.h>
#include <multigrid/mg_dof_accessor.h>
#include <multigrid/multigrid.h>
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


                                 // @sect3{The <code>LaplaceProblem</code> class template}

				 // This main class is basically the same
				 // class as in step-6. As far as member
				 // functions is concerned, the only addition
				 // is the <code>assemble_multigrid</code>
				 // function that assembles the matrices that
				 // correspond to the discrete operators on
				 // intermediate levels:
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

    ConstraintMatrix     constraints;

    Vector<double>       solution;
    Vector<double>       system_rhs;

    const unsigned int degree;

				     // The following three objects are the
				     // only additional member variables,
				     // compared to step-6. They represent the
				     // operators that act on individual
				     // levels of the multilevel hierarchy,
				     // rather than on the finest mesh as do
				     // the objects above.
				     //
				     // To facilitate having objects on each
				     // level of a multilevel hierarchy,
				     // deal.II has the MGLevelObject class
				     // template that provides storage for
				     // objects on each level. What we need
				     // here are matrices on each level, which
				     // implies that we also need sparsity
				     // patterns on each level. As outlined in
				     // the @ref mg_paper, the operators
				     // (matrices) that we need are actually
				     // twofold: one on the interior of each
				     // level, and one at the interface
				     // between each level and that part of
				     // the domain where the mesh is
				     // coarser. In fact, we will need the
				     // latter in two versions: for the
				     // direction from coarse to fine mesh and
				     // from fine to coarse. Fortunately,
				     // however, we here have a self-adjoint
				     // problem for which one of these is the
				     // transpose of the other, and so we only
				     // have to build one; we choose the one
				     // from coarse to fine.
    MGLevelObject<SparsityPattern>       mg_sparsity_patterns;
    MGLevelObject<SparseMatrix<double> > mg_matrices;
    MGLevelObject<SparseMatrix<double> > mg_interface_matrices;
};



                                 // @sect3{Nonconstant coefficients}

				 // The implementation of nonconstant
				 // coefficients is copied verbatim
				 // from step-5 and step-6:

template <int dim>
class Coefficient : public Function<dim>
{
  public:
    Coefficient () : Function<dim>() {}

    virtual double value (const Point<dim>   &p,
			  const unsigned int  component = 0) const;

    virtual void value_list (const std::vector<Point<dim> > &points,
			     std::vector<double>            &values,
			     const unsigned int              component = 0) const;
};



template <int dim>
double Coefficient<dim>::value (const Point<dim> &p,
				const unsigned int) const
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
  const unsigned int n_points = points.size();

  Assert (values.size() == n_points,
	  ExcDimensionMismatch (values.size(), n_points));

  Assert (component == 0,
	  ExcIndexRange (component, 0, 1));

  for (unsigned int i=0; i<n_points; ++i)
    values[i] = Coefficient<dim>::value (points[i]);
}


                                 // @sect3{The <code>LaplaceProblem</code> class implementation}

                                 // @sect4{LaplaceProblem::LaplaceProblem}

				 // The constructor is left mostly
				 // unchanged. We take the polynomial degree
				 // of the finite elements to be used as a
				 // constructor argument and store it in a
				 // member variable.
				 //
				 // By convention, all adaptively refined
				 // triangulations in deal.II never change by
				 // more than one level across a face between
				 // cells. For our multigrid algorithms,
				 // however, we need a slightly stricter
				 // guarantee, namely that the mesh also does
				 // not change by more than refinement level
				 // across vertices that might connect two
				 // cells. In other words, we must prevent the
				 // following situation:
				 //
				 // @image html limit_level_difference_at_vertices.png ""
				 //
				 // This is achieved by passing the
				 // Triangulation::limit_level_difference_at_vertices
				 // flag to the constructor of the
				 // triangulation class.
template <int dim>
LaplaceProblem<dim>::LaplaceProblem (const unsigned int degree)
		:
		triangulation (Triangulation<dim>::
			       limit_level_difference_at_vertices),
		fe (degree),
		mg_dof_handler (triangulation),
		degree(degree)
{}



                                 // @sect4{LaplaceProblem::setup_system}

				 // The following function extends what the
				 // corresponding one in step-6 did. The top
				 // part, apart from the additional output,
				 // does the same:
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

				   // But it starts to be a wee bit different
				   // here, although this still doesn't have
				   // anything to do with multigrid
				   // methods. step-6 took care of boundary
				   // values and hanging nodes in a separate
				   // step after assembling the global matrix
				   // from local contributions. This works,
				   // but the same can be done in a slightly
				   // simpler way if we already take care of
				   // these constraints at the time of copying
				   // local contributions into the global
				   // matrix. To this end, we here do not just
				   // compute the constraints do to hanging
				   // nodes, but also due to zero boundary
				   // conditions. Both kinds of constraints
				   // can be put into the same object
				   // (<code>constraints</code>), and we will
				   // use this set of constraints later on to
				   // help us copy local contributions
				   // correctly into the global linear system
				   // right away, without the need for a later
				   // clean-up stage:
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

				   // Now for the things that concern the
				   // multigrid data structures. First, we
				   // resize the multi-level objects to hold
				   // matrices and sparsity patterns for every
				   // level. The coarse level is zero (this is
				   // mandatory right now but may change in a
				   // future revision). Note that these
				   // functions take a complete, inclusive
				   // range here (not a starting index and
				   // size), so the finest level is
				   // <code>n_levels-1</code>.  We first have
				   // to resize the container holding the
				   // SparseMatrix classes, since they have to
				   // release their SparsityPattern before the
				   // can be destroyed upon resizing.
  const unsigned int n_levels = triangulation.n_levels();

  mg_interface_matrices.resize(0, n_levels-1);
  mg_interface_matrices.clear ();
  mg_matrices.resize(0, n_levels-1);
  mg_matrices.clear ();
  mg_sparsity_patterns.resize(0, n_levels-1);

				   // Now, we have to provide a matrix on each
				   // level. To this end, we first use the
				   // MGTools::make_sparsity_pattern function
				   // to first generate a preliminary
				   // compressed sparsity pattern on each
				   // level (see the @ref Sparsity module for
				   // more information on this topic) and then
				   // copy it over to the one we really
				   // want. The next step is to initialize
				   // both kinds of level matrices with these
				   // sparsity patterns.
				   //
				   // It may be worth pointing out that the
				   // interface matrices only have entries for
				   // degrees of freedom that sit at or next
				   // to the interface between coarser and
				   // finer levels of the mesh. They are
				   // therefore even sparser than the matrices
				   // on the individual levels of our
				   // multigrid hierarchy. If we were more
				   // concerned about memory usage (and
				   // possibly the speed with which we can
				   // multiply with these matrices), we should
				   // use separate and different sparsity
				   // patterns for these two kinds of
				   // matrices.
  for (unsigned int level=0; level<n_levels; ++level)
    {
      CompressedSparsityPattern csp;
      csp.reinit(mg_dof_handler.n_dofs(level),
		 mg_dof_handler.n_dofs(level));
      MGTools::make_sparsity_pattern(mg_dof_handler, csp, level);

      mg_sparsity_patterns[level].copy_from (csp);
      
      mg_matrices[level].reinit(mg_sparsity_patterns[level]);
      mg_interface_matrices[level].reinit(mg_sparsity_patterns[level]);
    }
}


                                 // @sect4{LaplaceProblem::assemble_system}

				 // The following function assembles the
				 // linear system on the finesh level of the
				 // mesh. It is almost exactly the same as in
				 // step-6, with the exception that we don't
				 // eliminate hanging nodes and boundary
				 // values after assembling, but while copying
				 // local contributions into the global
				 // matrix. This is not only simpler but also
				 // more efficient for large problems.
template <int dim>
void LaplaceProblem<dim>::assemble_system ()
{
  const QGauss<dim>  quadrature_formula(degree+1);

  FEValues<dim> fe_values (fe, quadrature_formula,
			   update_values    |  update_gradients |
			   update_quadrature_points  |  update_JxW_values);

  const unsigned int   dofs_per_cell = fe.dofs_per_cell;
  const unsigned int   n_q_points    = quadrature_formula.size();

  FullMatrix<double>   cell_matrix (dofs_per_cell, dofs_per_cell);
  Vector<double>       cell_rhs (dofs_per_cell);

  std::vector<unsigned int> local_dof_indices (dofs_per_cell);

  const Coefficient<dim> coefficient;
  std::vector<double>    coefficient_values (n_q_points);

  typename MGDoFHandler<dim>::active_cell_iterator
    cell = mg_dof_handler.begin_active(),
    endc = mg_dof_handler.end();
  for (; cell!=endc; ++cell)
    {
      cell_matrix = 0;
      cell_rhs = 0;

      fe_values.reinit (cell);

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
      constraints.distribute_local_to_global (cell_matrix, cell_rhs,
					      local_dof_indices,
					      system_matrix, system_rhs);
    }
}


                                 // @sect4{LaplaceProblem::assemble_multigrid}

				 // The next function is the one that builds
				 // the linear operators (matrices) that
				 // define the multigrid method on each level
				 // of the mesh. The integration core is the
				 // same as above, but the loop below will go
				 // over all existing cells instead of just
				 // the active ones, and the results must be
				 // entered into the correct matrix. Note also
				 // that since we only do multi-level
				 // preconditioning, no right-hand side needs
				 // to be assembled here.
				 //
				 // Before we go there, however, we have to
				 // take care of a significant amount of book
				 // keeping:
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

				   // Next a few things that are specific to
				   // building the multigrid data structures
				   // (since we only need them in the current
				   // function, rather than also elsewhere, we
				   // build them here instead of the
				   // <code>setup_system</code> function). 
  std::vector<std::vector<bool> > interface_dofs;
  std::vector<std::vector<bool> > boundary_interface_dofs;
  for (unsigned int level = 0; level<triangulation.n_levels(); ++level)
    {
      interface_dofs.push_back (std::vector<bool>
				(mg_dof_handler.n_dofs(level)));
      boundary_interface_dofs.push_back (std::vector<bool>
					 (mg_dof_handler.n_dofs(level)));
    }
  MGTools::extract_inner_interface_dofs (mg_dof_handler,
					 interface_dofs,
					 boundary_interface_dofs);

  typename FunctionMap<dim>::type      dirichlet_boundary;
  ZeroFunction<dim>                    homogeneous_dirichlet_bc (1);
  dirichlet_boundary[0] = &homogeneous_dirichlet_bc;

  std::vector<std::set<unsigned int> > boundary_indices (triangulation.n_levels());
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

  const Coefficient<dim> coefficient;
  std::vector<double>    coefficient_values (n_q_points);

  typename MGDoFHandler<dim>::cell_iterator cell = mg_dof_handler.begin(),
					    endc = mg_dof_handler.end();

  for (; cell!=endc; ++cell)
    {
      cell_matrix = 0;

				       // Compute the values specified
				       // by update flags above.
      fe_values.reinit (cell);

      coefficient.value_list (fe_values.get_quadrature_points(),
			      coefficient_values);

				       // This is exactly the
				       // integration loop of the cell
				       // matrix above.
      for (unsigned int q_point=0; q_point<n_q_points; ++q_point)
	for (unsigned int i=0; i<dofs_per_cell; ++i)
	  for (unsigned int j=0; j<dofs_per_cell; ++j)
	    cell_matrix(i,j) += (coefficient_values[q_point] *
				 fe_values.shape_grad(i,q_point) *
				 fe_values.shape_grad(j,q_point) *
				 fe_values.JxW(q_point));

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
				     mg_interface_matrices[level]);
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
  MGMatrix<SparseMatrix<double>, Vector<double> >
    mg_matrix(&mg_matrices);
				   //do the same for the interface matrices
  MGMatrix<SparseMatrix<double>, Vector<double> >
    mg_interface_up(&mg_interface_matrices);
  MGMatrix<SparseMatrix<double>, Vector<double> >
    mg_interface_down(&mg_interface_matrices);
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



                                 // @sect4{Postprocessing}

				 // The following two functions postprocess a
				 // solution once it is computed. In
				 // particular, the first one refines the mesh
				 // at the beginning of each cycle while the
				 // second one outputs results at the end of
				 // each such cycle. The functions are almost
				 // unchanged from those in step-6, with the
				 // exception of two minor differences: The
				 // KellyErrorEstimator::estimate function
				 // wants an argument of type DoFHandler, not
				 // MGDoFHandler, and so we have to cast from
				 // derived to base class; and we generate
				 // output in VTK format, to use the more
				 // modern visualization programs available
				 // today compared to those that were
				 // available when step-6 was written.
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


                                 // @sect4{LaplaceProblem::run}

				 // Like several of the functions above, this
				 // is almost exactly a copy of of the
				 // corresponding function in step-6. The only
				 // difference is the call to
				 // <code>assemble_multigrid</code> that takes
				 // care of forming the matrices on every
				 // level that we need in the multigrid
				 // method.
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
		<< " (by level: ";
      for (unsigned int level=0; level<triangulation.n_levels(); ++level)
	std::cout << mg_dof_handler.n_dofs(level)
		  << (level == triangulation.n_levels()-1
		      ? ")" : ", ");
      std::cout << std::endl;

      assemble_system ();
      assemble_multigrid ();

      solve ();
      output_results (cycle);
    }
}


				 // @sect3{The main() function}
				 //
				 // This is again the same function as
				 // in step-6:
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
