/* $Id$ */
/* Author: Guido Kanschat, Texas A&M University, 2009 */

/*    $Id$       */
/*                                                                */
/*    Copyright (C) 2001, 2002, 2003, 2004, 2005, 2006, 2007, 2009 by the deal.II authors */
/*                                                                */
/*    This file is subject to QPL and may not be  distributed     */
/*    without copyright and license information. Please refer     */
/*    to the file deal.II/doc/license.html for the  text  and     */
/*    further information on this license.                        */

				 // The first few files have already
				 // been covered in example 12
				 // and will thus not be further
				 // commented on.
#include <base/quadrature_lib.h>
#include <base/function.h>
#include <lac/vector.h>
#include <lac/sparse_matrix.h>
#include <grid/tria.h>
#include <grid/grid_generator.h>
#include <grid/grid_out.h>
#include <grid/grid_refinement.h>
#include <grid/tria_accessor.h>
#include <grid/tria_iterator.h>
#include <fe/fe_values.h>
#include <dofs/dof_handler.h>
#include <dofs/dof_accessor.h>
#include <dofs/dof_tools.h>
#include <numerics/data_out.h>
#include <fe/mapping_q1.h>
#include <fe/fe_dgq.h>
#include <lac/solver_richardson.h>
#include <lac/precondition_block.h>
#include <numerics/derivative_approximation.h>
#include <base/timer.h>

				 // Here come the new include files
				 // for using the MeshWorker framework
#include <numerics/mesh_worker.h>
#include <numerics/mesh_worker_loop.h>

#include <iostream>
#include <fstream>

				 // The last step is as in all
				 // previous programs:
using namespace dealii;

				 // @sect3{Equation data}
				 //
				 // First we define the classes
				 // representing the equation-specific
				 // functions. Both classes, <code>RHS</code>
				 // and <code>BoundaryValues</code>, are
				 // derived from the <code>Function</code>
				 // class. Only the <code>value_list</code>
				 // function are implemented because
				 // only lists of function values are
				 // computed rather than single
				 // values.
template <int dim>
class RHS:  public Function<dim>
{
  public:
    RHS () {};
    virtual void value_list (const std::vector<Point<dim> > &points,
			     std::vector<double> &values,
			     const unsigned int component=0) const;
};


template <int dim>
class BoundaryValues:  public Function<dim>
{
  public:
    BoundaryValues () {};
    virtual void value_list (const std::vector<Point<dim> > &points,
			     std::vector<double> &values,
			     const unsigned int component=0) const;
};


				 // The class <code>Beta</code> represents the
				 // vector valued flow field of the
				 // linear transport equation and is
				 // not derived from the <code>Function</code>
				 // class as we prefer to get function
				 // values of type <code>Point</code> rather
				 // than of type
				 // <code>Vector@<double@></code>. This, because
				 // there exist scalar products
				 // between <code>Point</code> and <code>Point</code> as
				 // well as between <code>Point</code> and
				 // <code>Tensor</code>, simplifying terms like
				 // $\beta\cdot n$ and
				 // $\beta\cdot\nabla v$.
                                 //
                                 // An unnecessary empty constructor
                                 // is added to the class to work
                                 // around a bug in Compaq's cxx
                                 // compiler which otherwise reports
                                 // an error about an omitted
                                 // initializer for an object of
                                 // this class further down.
template <int dim>
class Beta
{
  public:
    Beta () {}
    void value_list (const std::vector<Point<dim> > &points,
		     std::vector<Point<dim> > &values) const;
};


				 // The implementation of the
				 // <code>value_list</code> functions of these
				 // classes are rather simple.  For
				 // simplicity the right hand side is
				 // set to be zero but will be
				 // assembled anyway.
template <int dim>
void RHS<dim>::value_list(const std::vector<Point<dim> > &points,
			  std::vector<double> &values,
			  const unsigned int) const
{
				   // Usually we check whether input
				   // parameter have the right sizes.
  Assert(values.size()==points.size(),
	 ExcDimensionMismatch(values.size(),points.size()));

  for (unsigned int i=0; i<values.size(); ++i)
    values[i]=0;
}


				 // The flow field is chosen to be
				 // circular, counterclockwise, and with
				 // the origin as midpoint.
template <int dim>
void Beta<dim>::value_list(const std::vector<Point<dim> > &points,
			   std::vector<Point<dim> > &values) const
{
  Assert(values.size()==points.size(),
	 ExcDimensionMismatch(values.size(),points.size()));

  for (unsigned int i=0; i<points.size(); ++i)
    {
      values[i](0) = -points[i](1);
      values[i](1) = points[i](0);
      values[i] /= std::sqrt(values[i].square());
    }
}


				 // Hence the inflow boundary of the
				 // unit square [0,1]^2 are the right
				 // and the lower boundaries. We
				 // prescribe discontinuous boundary
				 // values 1 and 0 on the x-axis and
				 // value 0 on the right boundary. The
				 // values of this function on the
				 // outflow boundaries will not be
				 // used within the DG scheme.
template <int dim>
void BoundaryValues<dim>::value_list(const std::vector<Point<dim> > &points,
				       std::vector<double> &values,
				       const unsigned int) const
{
  Assert(values.size()==points.size(),
	 ExcDimensionMismatch(values.size(),points.size()));

  for (unsigned int i=0; i<values.size(); ++i)
    {
      if (points[i](0)<0.5)
	values[i]=1.;
      else
	values[i]=0.;
    }
}


				 // @sect3{Integrating cell and face matrices}
				 //
				 // We define a class that fits into
				 // the MeshWorker framework. Since it
				 // will be used by
				 // MeshWorker::AssemblingIntegrator,
				 // it needs the functions for cell,
				 // boundary and interior face
				 // integration specified exactly as
				 // below.

				 // The base class Subscriptor is
				 // needed so that
				 // MeshWorker::AssemblingIntegrator
				 // can store a SmartPointer to an
				 // object of this class.
template <int dim>
class DGIntegrator : public Subscriptor
{
  public:
				     // First, we define the types of
				     // the two info objects handed to
				     // the local integration
				     // functions in order to make our
				     // life easier below.
    typedef typename MeshWorker::IntegrationWorker<dim>::CellInfo CellInfo;
    typedef typename MeshWorker::IntegrationWorker<dim>::FaceInfo FaceInfo;

				     // The following three functions
				     // are the ones that get called
				     // inside the generic loop over all
				     // cells and faces. They are the
				     // ones doing the actual
				     // integration.
    void cell(CellInfo& info) const;
    void bdry(FaceInfo& info) const;
    void face(FaceInfo& info1, FaceInfo& info2) const;

				     // Additionally, like in step-12,
				     // we have objects of the
				     // functions used in this class.
  private:
    const Beta<dim> beta_function;
    const RHS<dim> rhs_function;
    const BoundaryValues<dim> boundary_function;
};

				 // @sect4{The local integrators}

				 // These functions are analogous to
				 // step-12 and differ only in the
				 // data structures. Instead of
				 // providing the local matrices
				 // explicitly in the argument list,
				 // they are part of the info object.

				 // Note that here we still have the
				 // local integration loop inside the
				 // following functions. The program
				 // would be even shorter, if we used
				 // pre-made operators from the
				 // Operators namespace (which will be
				 // added soon).

template <int dim>
void DGIntegrator<dim>::cell(CellInfo& info) const
{
				   // First, let us retrieve some of
				   // the objects used here from
				   // @p info. Note that these objects
				   // can handle much more complex
				   // structures, thus the access here
				   // looks more complicated than
				   // might seem necessary.
  const FEValuesBase<dim>& fe_v = info.fe();
  FullMatrix<double>& local_matrix = info.M1[0].matrix;
  Vector<double>& local_vector = info.R[0].block(0);

				   // With these objects, we continue
				   // local integration like in step-12.
  const std::vector<double> &JxW = fe_v.get_JxW_values ();
  
  std::vector<Point<dim> > beta (fe_v.n_quadrature_points);
  std::vector<double> rhs (fe_v.n_quadrature_points);
  
  beta_function.value_list (fe_v.get_quadrature_points(), beta);
  rhs_function.value_list (fe_v.get_quadrature_points(), rhs);
  
  for (unsigned int point=0; point<fe_v.n_quadrature_points; ++point)
    for (unsigned int i=0; i<fe_v.dofs_per_cell; ++i) 
      {
	for (unsigned int j=0; j<fe_v.dofs_per_cell; ++j)
	  local_matrix(i,j) -= beta[point]*fe_v.shape_grad(i,point)*
			      fe_v.shape_value(j,point) *
			      JxW[point];
	
	local_vector(i) += rhs[point] * fe_v.shape_value(i,point) * JxW[point];
      }
}

				 // Now the same for the boundary
				 // terms. Note that now we use
				 // FEFaceValuesBase in order to get
				 // normal vectors.
template <int dim>
void DGIntegrator<dim>::bdry(FaceInfo& info) const
{
  const FEFaceValuesBase<dim>& fe_v = info.fe();
  FullMatrix<double>& local_matrix = info.M1[0].matrix;
  Vector<double>& local_vector = info.R[0].block(0);

  const std::vector<double> &JxW = fe_v.get_JxW_values ();
  const std::vector<Point<dim> > &normals = fe_v.get_normal_vectors ();
  
  std::vector<Point<dim> > beta (fe_v.n_quadrature_points);
  std::vector<double> g(fe_v.n_quadrature_points);
  
  beta_function.value_list (fe_v.get_quadrature_points(), beta);
  boundary_function.value_list (fe_v.get_quadrature_points(), g);

  for (unsigned int point=0; point<fe_v.n_quadrature_points; ++point)
    {
      const double beta_n=beta[point] * normals[point];      
      if (beta_n>0)
	for (unsigned int i=0; i<fe_v.dofs_per_cell; ++i)
	  for (unsigned int j=0; j<fe_v.dofs_per_cell; ++j)
	    local_matrix(i,j) += beta_n *
				 fe_v.shape_value(j,point) *
				 fe_v.shape_value(i,point) *
				 JxW[point];
      else
	for (unsigned int i=0; i<fe_v.dofs_per_cell; ++i)
	  local_vector(i) -= beta_n *
			     g[point] *
			     fe_v.shape_value(i,point) *
			     JxW[point];
    }
}

				 // Finally, the interior face
				 // terms. The difference here is that
				 // we receive two info objects, one
				 // for each cell adjacent to the face
				 // and we assemble four matrices, one
				 // for each cell and two for coupling
				 // back and forth.
template <int dim>
void DGIntegrator<dim>::face(FaceInfo& info1, FaceInfo& info2) const
{
				   // For quadrature points, weights,
				   // etc., we use the
				   // FEFaceValuesBase object of the
				   // first argument.
  const FEFaceValuesBase<dim>& fe_v = info1.fe();

				   // For additional shape functions,
				   // we have to ask the neighbors
				   // FEFaceValuesBase.
  const FEFaceValuesBase<dim>& fe_v_neighbor = info2.fe();
  
				   // Then we get references to the
				   // four local matrices. The letters
				   // u and v refer to trial and test
				   // functions, respectively. The
				   // numbers indicate the cells
				   // provided by info1 and info2. By
				   // convention, the two matrices in
				   // each info object refer to the
				   // test functions on the respective
				   // cell. The first matrix contains the
				   // interior couplings of that cell,
				   // while the second contains the
				   // couplings between cells.
  FullMatrix<double>& u1_v1_matrix = info1.M1[0].matrix;
  FullMatrix<double>& u2_v1_matrix = info1.M2[0].matrix;
  FullMatrix<double>& u1_v2_matrix = info2.M2[0].matrix;
  FullMatrix<double>& u2_v2_matrix = info2.M1[0].matrix;
  
				   // Here, following the previous
				   // functions, we would have the
				   // local right hand side
				   // vectors. Fortunately, the
				   // interface terms only involve the
				   // solution and the right hand side
				   // does not obtain a contribution.
  
  const std::vector<double> &JxW = fe_v.get_JxW_values ();
  const std::vector<Point<dim> > &normals = fe_v.get_normal_vectors ();

  std::vector<Point<dim> > beta (fe_v.n_quadrature_points);
  
  beta_function.value_list (fe_v.get_quadrature_points(), beta);

  for (unsigned int point=0; point<fe_v.n_quadrature_points; ++point)
    {
      const double beta_n=beta[point] * normals[point];
      if (beta_n>0)
	{
					   // This term we've already
					   // seen.
	  for (unsigned int i=0; i<fe_v.dofs_per_cell; ++i)
	    for (unsigned int j=0; j<fe_v.dofs_per_cell; ++j)
	      u1_v1_matrix(i,j) += beta_n *
				 fe_v.shape_value(j,point) *
				 fe_v.shape_value(i,point) *
				 JxW[point];

					   // We additionally assemble
					   // the term $(\beta\cdot n
					   // u,\hat v)_{\partial
					   // \kappa_+}$,
	  for (unsigned int k=0; k<fe_v_neighbor.dofs_per_cell; ++k)
	    for (unsigned int j=0; j<fe_v.dofs_per_cell; ++j)
	      u1_v2_matrix(k,j) -= beta_n *
				  fe_v.shape_value(j,point) *
				  fe_v_neighbor.shape_value(k,point) *
				  JxW[point];
	}
      else
	{
					   // This one we've already
					   // seen, too.
	  for (unsigned int i=0; i<fe_v.dofs_per_cell; ++i)
	    for (unsigned int l=0; l<fe_v_neighbor.dofs_per_cell; ++l)
	      u2_v1_matrix(i,l) += beta_n *
				  fe_v_neighbor.shape_value(l,point) *
				  fe_v.shape_value(i,point) *
				  JxW[point];

					   // And this is another new
					   // one: $(\beta\cdot n \hat
					   // u,\hat v)_{\partial
					   // \kappa_-}$.
	  for (unsigned int k=0; k<fe_v_neighbor.dofs_per_cell; ++k)
	    for (unsigned int l=0; l<fe_v_neighbor.dofs_per_cell; ++l)
	      u2_v2_matrix(k,l) -= beta_n *
				   fe_v_neighbor.shape_value(l,point) *
				   fe_v_neighbor.shape_value(k,point) *
				   JxW[point];
	}
    }
}


				 // @sect3{Class: DGMethod}
				 //
				 // After these preparations, we
				 // proceed with the main part of this
				 // program. The main class, here
				 // called <code>DGMethod</code> is basically
				 // the main class of step-6. One of
				 // the differences is that there's no
				 // ConstraintMatrix object. This is,
				 // because there are no hanging node
				 // constraints in DG discretizations.
template <int dim>
class DGMethod
{
  public:
    DGMethod ();
    ~DGMethod ();

    void run ();
    
  private:
    void setup_system ();
    void assemble_system ();
    void solve (Vector<double> &solution);
    void refine_grid ();
    void output_results (const unsigned int cycle) const;
    
    Triangulation<dim>   triangulation;
    const MappingQ1<dim> mapping;
    
				     // Furthermore we want to use DG
				     // elements of degree 1 (but this
				     // is only specified in the
				     // constructor). If you want to
				     // use a DG method of a different
				     // degree the whole program stays
				     // the same, only replace 1 in
				     // the constructor by the wanted
				     // degree.
    FE_DGQ<dim>          fe;
    DoFHandler<dim>      dof_handler;

    SparsityPattern      sparsity_pattern;
    SparseMatrix<double> system_matrix;

				     // And there are two solution
				     // vectors, that store the
				     // solutions to the problems
				     // corresponding to the two
				     // different assembling routines
				     // <code>assemble_system1</code> and
				     // <code>assemble_system2</code>;
    Vector<double>       solution;
    Vector<double>       right_hand_side;    
};


template <int dim>
DGMethod<dim>::DGMethod ()
		:
						 // Change here for DG
						 // methods of
						 // different degrees.
                fe (1),
		dof_handler (triangulation)
{}


template <int dim>
DGMethod<dim>::~DGMethod () 
{
  dof_handler.clear ();
}


template <int dim>
void DGMethod<dim>::setup_system ()
{
				   // First we need to distribute the
				   // DoFs.
  dof_handler.distribute_dofs (fe);

				   // The DoFs of a cell are coupled
				   // with all DoFs of all neighboring
				   // cells.  Therefore the maximum
				   // number of matrix entries per row
				   // is needed when all neighbors of
				   // a cell are once more refined
				   // than the cell under
				   // consideration.
  sparsity_pattern.reinit (dof_handler.n_dofs(),
			   dof_handler.n_dofs(),
			   (GeometryInfo<dim>::faces_per_cell
			    *GeometryInfo<dim>::max_children_per_face+1)*fe.dofs_per_cell);
  
				   // For DG discretizations we call
				   // the function analogue to
				   // DoFTools::make_sparsity_pattern.
  DoFTools::make_flux_sparsity_pattern (dof_handler, sparsity_pattern);
  
				   // All following function calls are
				   // already known.
  sparsity_pattern.compress();
  
  system_matrix.reinit (sparsity_pattern);

  solution.reinit (dof_handler.n_dofs());
  right_hand_side.reinit (dof_handler.n_dofs());
}

				 // @sect4{Function: assemble_system}
				 // Here we see the major difference
				 // to assembling by hand. Instead of
				 // writing loops over cells and
				 // faces, we leave all this to the
				 // MeshWorker framework. In order to
				 // do so, we just have to define
				 // local integration objects and use
				 // one of the classes in Assembler to
				 // build the global system. 
template <int dim>
void DGMethod<dim>::assemble_system () 
{
				   // Here we generate an object of
				   // our own integration class, which
				   // knows how to compute cell and
				   // face contributions for the
				   // matrix and the residual.
  const DGIntegrator<dim> dg;

				   // This is the magic object, which
				   // knows everything about the data
				   // structures and local integration
				   // (the latter through our object
				   // @p dg). This is the object doing
				   // the work in the function
				   // MeshWorker::loop(), which is
				   // implicitly called by
				   // MeshWorker::integration_loop()
				   // below.
				   // After @p dg did the local
				   // integration, the
				   // MeshWorker::Assembler::SystemSimple
				   // object distributes these into
				   // the global sparse matrix and the
				   // right hand side vector.
				   //
				   // It turns out,
				   // MeshWorker::AssemblingIntegrator
				   // itself is not that clever at
				   // all, but all its capabilities
				   // are provided by the two later
				   // template arguments. By
				   // exchanging
				   // MeshWorker::Assembler::SystemSimple,
				   // we could for instance assemble a
				   // BlockMatrix or just a Vector
				   // instead.
  MeshWorker::AssemblingIntegrator<dim, MeshWorker::Assembler::SystemSimple<SparseMatrix<double>, Vector<double> >, DGIntegrator<dim> >
    integrator(dg);
  
				   // First, we initialize the
				   // quadrature formulae and the
				   // update flags in the worker base
				   // class. For quadrature, we play
				   // safe and use a QGauss formula
				   // with number of points one higher
				   // than the polynomial degree
				   // used. Since the quadratures for
				   // cells, boundary and interior
				   // faces can be selected
				   // independently, we have to hand
				   // over this value three times.
  const unsigned int n_gauss_points = dof_handler.get_fe().degree+1;
  integrator.initialize_gauss_quadrature(n_gauss_points, n_gauss_points, n_gauss_points);

				   // These are the types of values we
				   // need for integrating our
				   // system. They are added to the
				   // flags used on cells, boundary
				   // and interior faces, as well as
				   // interior neighbor faces, which is
				   // forced by the four @p true values.
  UpdateFlags update_flags = update_quadrature_points | update_values | update_gradients;
  integrator.add_update_flags(update_flags, true, true, true, true);

				   // Finally, we have to tell the
				   // assembler base class, where to
				   // put the local data. These will
				   // be our system matrix and the
				   // right hand side.
  integrator.initialize(system_matrix, right_hand_side);

				   // Finally, we get to the
				   // integration loop. @p info_box is
				   // an object, that generates the
				   // extended iterators for cells and
				   // faces of type
				   // MeshWorker::IntegrationInfo. Since
				   // we need five different of them,
				   // this is a handy shortcut. It
				   // receives all the stuff we
				   // created so far.
  MeshWorker::IntegrationInfoBox<dim> info_box(dof_handler);
  info_box.initialize(integrator, fe, mapping);

				   // Finally, the integration loop
				   // over all active cells
				   // (determined by the first
				   // argument, which is an active iterator).
  MeshWorker::integration_loop(dof_handler.begin_active(), dof_handler.end(), info_box, integrator);
}


				 // @sect3{All the rest}
				 //
				 // For this simple problem we use the
				 // simplest possible solver, called
				 // Richardson iteration, that
				 // represents a simple defect
				 // correction. This, in combination
				 // with a block SSOR preconditioner,
				 // that uses the special block matrix
				 // structure of system matrices
				 // arising from DG
				 // discretizations. The size of these
				 // blocks are the number of DoFs per
				 // cell. Here, we use a SSOR
				 // preconditioning as we have not
				 // renumbered the DoFs according to
				 // the flow field. If the DoFs are
				 // renumbered downstream the flow,
				 // then a block Gauss-Seidel
				 // preconditioner (see the
				 // PreconditionBlockSOR class with
				 // relaxation=1) makes a much better
				 // job.
template <int dim>
void DGMethod<dim>::solve (Vector<double> &solution) 
{
  SolverControl           solver_control (1000, 1e-12, false, false);
  SolverRichardson<>      solver (solver_control);

				   // Here we create the
				   // preconditioner,
  PreconditionBlockSSOR<SparseMatrix<double> > preconditioner;

				   // we assigned the matrix to it and
				   // set the right block size.
  preconditioner.initialize(system_matrix, fe.dofs_per_cell);

				   // After these preparations we are
				   // ready to start the linear solver.
  solver.solve (system_matrix, solution, right_hand_side,
		preconditioner);
}


				 // We refine the grid according to a
				 // very simple refinement criterion,
				 // namely an approximation to the
				 // gradient of the solution. As here
				 // we consider the DG(1) method
				 // (i.e. we use piecewise bilinear
				 // shape functions) we could simply
				 // compute the gradients on each
				 // cell. But we do not want to base
				 // our refinement indicator on the
				 // gradients on each cell only, but
				 // want to base them also on jumps of
				 // the discontinuous solution
				 // function over faces between
				 // neighboring cells. The simplest
				 // way of doing that is to compute
				 // approximative gradients by
				 // difference quotients including the
				 // cell under consideration and its
				 // neighbors. This is done by the
				 // <code>DerivativeApproximation</code> class
				 // that computes the approximate
				 // gradients in a way similar to the
				 // <code>GradientEstimation</code> described
				 // in step-9 of this tutorial. In
				 // fact, the
				 // <code>DerivativeApproximation</code> class
				 // was developed following the
				 // <code>GradientEstimation</code> class of
				 // step-9. Relating to the
				 // discussion in step-9, here we
				 // consider $h^{1+d/2}|\nabla_h
				 // u_h|$. Furthermore we note that we
				 // do not consider approximate second
				 // derivatives because solutions to
				 // the linear advection equation are
				 // in general not in $H^2$ but in $H^1$
				 // (to be more precise, in $H^1_\beta$)
				 // only.
template <int dim>
void DGMethod<dim>::refine_grid ()
{
				   // The <code>DerivativeApproximation</code>
				   // class computes the gradients to
				   // float precision. This is
				   // sufficient as they are
				   // approximate and serve as
				   // refinement indicators only.
  Vector<float> gradient_indicator (triangulation.n_active_cells());

				   // Now the approximate gradients
				   // are computed
  DerivativeApproximation::approximate_gradient (mapping,
						 dof_handler,
						 solution,
						 gradient_indicator);

				   // and they are cell-wise scaled by
				   // the factor $h^{1+d/2}$
  typename DoFHandler<dim>::active_cell_iterator
    cell = dof_handler.begin_active(),
    endc = dof_handler.end();
  for (unsigned int cell_no=0; cell!=endc; ++cell, ++cell_no)
    gradient_indicator(cell_no)*=std::pow(cell->diameter(), 1+1.0*dim/2);

				   // Finally they serve as refinement
				   // indicator.
  GridRefinement::refine_and_coarsen_fixed_number (triangulation,
						   gradient_indicator,
						   0.3, 0.1);

  triangulation.execute_coarsening_and_refinement ();
}


				 // The output of this program
				 // consists of eps-files of the
				 // adaptively refined grids and the
				 // numerical solutions given in
				 // gnuplot format. This was covered
				 // in previous examples and will not
				 // be further commented on.
template <int dim>
void DGMethod<dim>::output_results (const unsigned int cycle) const
{
				   // Write the grid in eps format.
  std::string filename = "grid-";
  filename += ('0' + cycle);
  Assert (cycle < 10, ExcInternalError());
  
  filename += ".eps";
  std::cout << "Writing grid to <" << filename << ">..." << std::endl;
  std::ofstream eps_output (filename.c_str());

  GridOut grid_out;
  grid_out.write_eps (triangulation, eps_output);
  
				   // Output of the solution in
				   // gnuplot format.
  filename = "sol-";
  filename += ('0' + cycle);
  Assert (cycle < 10, ExcInternalError());
  
  filename += ".gnuplot";
  std::cout << "Writing solution to <" << filename << ">..."
	    << std::endl << std::endl;
  std::ofstream gnuplot_output (filename.c_str());
  
  DataOut<dim> data_out;
  data_out.attach_dof_handler (dof_handler);
  data_out.add_data_vector (solution, "u");

  data_out.build_patches ();
  
  data_out.write_gnuplot(gnuplot_output);
}


				 // The following <code>run</code> function is
				 // similar to previous examples. The
				 // only difference is that the
				 // problem is assembled and solved
				 // twice on each refinement step;
				 // first by <code>assemble_system1</code> that
				 // implements the first version and
				 // then by <code>assemble_system2</code> that
				 // implements the second version of
				 // writing the DG
				 // discretization. Furthermore the
				 // time needed by each of the two
				 // assembling routines is measured.
template <int dim>
void DGMethod<dim>::run () 
{
  for (unsigned int cycle=0; cycle<6; ++cycle)
    {
      std::cout << "Cycle " << cycle << ':' << std::endl;

      if (cycle == 0)
	{
	  GridGenerator::hyper_cube (triangulation);

	  triangulation.refine_global (3);
	}
      else
	refine_grid ();
      

      std::cout << "   Number of active cells:       "
		<< triangulation.n_active_cells()
		<< std::endl;

      setup_system ();

      std::cout << "   Number of degrees of freedom: "
		<< dof_handler.n_dofs()
		<< std::endl;

      assemble_system ();
      solve (solution);

				       // Finally we perform the
				       // output.
      output_results (cycle);
    }
}

				 // The following <code>main</code> function is
				 // similar to previous examples and
				 // need not to be commented on.
int main () 
{
  try
    {
      DGMethod<2> dgmethod;
      dgmethod.run ();
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
    };
  
  return 0;
}


