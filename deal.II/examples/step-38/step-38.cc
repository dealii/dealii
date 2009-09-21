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
    virtual void value_list (const std::vector<Point<dim> > &points,
			     std::vector<double> &values,
			     const unsigned int component=0) const;
};


template <int dim>
class BoundaryValues:  public Function<dim>
{
  public:
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


				 // @sect3{Class: DGTransportEquation}
				 //
				 // Next we define the
				 // equation-dependent and
				 // DG-method-dependent class
				 // <code>DGTransportEquation</code>. Its
				 // member functions were already
				 // mentioned in the Introduction and
				 // will be explained
				 // below. Furthermore it includes
				 // objects of the previously defined
				 // <code>Beta</code>, <code>RHS</code> and
				 // <code>BoundaryValues</code> function
				 // classes.
template <int dim>
class DGTransportEquation
{
  public:
    DGTransportEquation();

    void assemble_cell_term(const FEValues<dim>& fe_v,
			    FullMatrix<double> &ui_vi_matrix,
			    Vector<double> &cell_vector) const;
    
    void assemble_boundary_term(const FEFaceValues<dim>& fe_v,
				FullMatrix<double> &ui_vi_matrix,
				Vector<double> &cell_vector) const;
    
    void assemble_face_term1(const FEFaceValuesBase<dim>& fe_v,
			     const FEFaceValuesBase<dim>& fe_v_neighbor,
			     FullMatrix<double> &ui_vi_matrix,
			     FullMatrix<double> &ue_vi_matrix) const;

    void assemble_face_term2(const FEFaceValuesBase<dim>& fe_v,
			     const FEFaceValuesBase<dim>& fe_v_neighbor,
			     FullMatrix<double> &ui_vi_matrix,
			     FullMatrix<double> &ue_vi_matrix,
			     FullMatrix<double> &ui_ve_matrix,
			     FullMatrix<double> &ue_ve_matrix) const;
  private:
    const Beta<dim> beta_function;
    const RHS<dim> rhs_function;
    const BoundaryValues<dim> boundary_function;
};


template <int dim>
DGTransportEquation<dim>::DGTransportEquation ()
		:
		beta_function (),
		rhs_function (),
		boundary_function ()
{}


				 // @sect4{Function: assemble_cell_term}
				 //
				 // The <code>assemble_cell_term</code>
				 // function assembles the cell terms
				 // of the discretization.
				 // <code>ui_vi_matrix</code> is a cell matrix,
				 // i.e. for a DG method of degree 1,
				 // it is of size 4 times 4, and
				 // <code>cell_vector</code> is of size 4.
				 // When this function is invoked,
				 // <code>fe_v</code> is already reinit'ed with the
				 // current cell before and includes
				 // all shape values needed.
template <int dim>
void DGTransportEquation<dim>::assemble_cell_term(
  const FEValues<dim> &fe_v,
  FullMatrix<double> &ui_vi_matrix,
  Vector<double> &cell_vector) const
{
				   // First we ask <code>fe_v</code> for the
				   // quadrature weights,
  const std::vector<double> &JxW = fe_v.get_JxW_values ();

				   // Then the flow field beta and the
				   // <code>rhs_function</code> are evaluated at
				   // the quadrature points,
  std::vector<Point<dim> > beta (fe_v.n_quadrature_points);
  std::vector<double> rhs (fe_v.n_quadrature_points);
  
  beta_function.value_list (fe_v.get_quadrature_points(), beta);
  rhs_function.value_list (fe_v.get_quadrature_points(), rhs);
  
				   // and the cell matrix and cell
				   // vector are assembled due to the
				   // terms $-(u,\beta\cdot\nabla
				   // v)_K$ and $(f,v)_K$.
  for (unsigned int point=0; point<fe_v.n_quadrature_points; ++point)
    for (unsigned int i=0; i<fe_v.dofs_per_cell; ++i) 
      {
	for (unsigned int j=0; j<fe_v.dofs_per_cell; ++j)
	  ui_vi_matrix(i,j) -= beta[point]*fe_v.shape_grad(i,point)*
			      fe_v.shape_value(j,point) *
			      JxW[point];
	
	cell_vector(i) += rhs[point] * fe_v.shape_value(i,point) * JxW[point];
      }
}


				 // @sect4{Function: assemble_boundary_term}
				 //
				 // The <code>assemble_boundary_term</code>
				 // function assembles the face terms
				 // at boundary faces.  When this
				 // function is invoked, <code>fe_v</code> is
				 // already reinit'ed with the current
				 // cell and current face. Hence it
				 // provides the shape values on that
				 // boundary face.
template <int dim>
void DGTransportEquation<dim>::assemble_boundary_term(
  const FEFaceValues<dim>& fe_v,    
  FullMatrix<double> &ui_vi_matrix,
  Vector<double> &cell_vector) const
{
				   // Again, as in the previous
				   // function, we ask the
				   // <code>FEValues</code> object for the
				   // quadrature weights
  const std::vector<double> &JxW = fe_v.get_JxW_values ();
				   // but here also for the normals.
  const std::vector<Point<dim> > &normals = fe_v.get_normal_vectors ();

				   // We evaluate the flow field
				   // and the boundary values at the
				   // quadrature points.
  std::vector<Point<dim> > beta (fe_v.n_quadrature_points);
  std::vector<double> g(fe_v.n_quadrature_points);
  
  beta_function.value_list (fe_v.get_quadrature_points(), beta);
  boundary_function.value_list (fe_v.get_quadrature_points(), g);

				   // Then we assemble cell vector and
				   // cell matrix according to the DG
				   // method given in the
				   // introduction.
  for (unsigned int point=0; point<fe_v.n_quadrature_points; ++point)
    {
      const double beta_n=beta[point] * normals[point];      
					 // We assemble the term
					 // $(\beta\cdot n
					 // u,v)_{\partial\kappa_+}$,
      if (beta_n>0)
	for (unsigned int i=0; i<fe_v.dofs_per_cell; ++i)
	  for (unsigned int j=0; j<fe_v.dofs_per_cell; ++j)
	    ui_vi_matrix(i,j) += beta_n *
			       fe_v.shape_value(j,point) *
			       fe_v.shape_value(i,point) *
			       JxW[point];
      else
					 // and the term $(\beta\cdot
					 // n g,v)_{\partial
					 // \kappa_-\cap\partial\Omega}$,
	for (unsigned int i=0; i<fe_v.dofs_per_cell; ++i)
	  cell_vector(i) -= beta_n *
			    g[point] *
			    fe_v.shape_value(i,point) *
			    JxW[point];
    }
}


				 // @sect4{Function: assemble_face_term1}
				 //
				 // The <code>assemble_face_term1</code>
				 // function assembles the face terms
				 // corresponding to the first version
				 // of the DG method, cf. above. For
				 // that case, the face terms are
				 // given as a sum of integrals over
				 // all cell boundaries.
				 //
				 // When this function is invoked,
				 // <code>fe_v</code> and <code>fe_v_neighbor</code> are
				 // already reinit'ed with the current
				 // cell and the neighoring cell,
				 // respectively, as well as with the
				 // current face. Hence they provide
				 // the inner and outer shape values
				 // on the face.
				 //
				 // In addition to the cell matrix
				 // <code>ui_vi_matrix</code> this function
				 // gets a new argument
				 // <code>ue_vi_matrix</code>, that stores
				 // contributions to the system matrix
				 // that are based on exterior values
				 // of $u$ and interior values of
				 // $v$. Here we note that <code>ue</code> is
				 // the short notation for <code>u
				 // exterior</code> and represents $u_h^-$,
				 // see the introduction.
template <int dim>
void DGTransportEquation<dim>::assemble_face_term1(
  const FEFaceValuesBase<dim>& fe_v,
  const FEFaceValuesBase<dim>& fe_v_neighbor,      
  FullMatrix<double> &ui_vi_matrix,
  FullMatrix<double> &ue_vi_matrix) const
{
				   // Again, as in the previous
				   // function, we ask the FEValues
				   // objects for the quadrature
				   // weights and the normals
  const std::vector<double> &JxW = fe_v.get_JxW_values ();
  const std::vector<Point<dim> > &normals = fe_v.get_normal_vectors ();

				   // and we evaluate the flow field
				   // at the quadrature points.
  std::vector<Point<dim> > beta (fe_v.n_quadrature_points);
  beta_function.value_list (fe_v.get_quadrature_points(), beta);

				   // Then we assemble the cell
				   // matrices according to the DG
				   // method given in the
				   // introduction.
  for (unsigned int point=0; point<fe_v.n_quadrature_points; ++point)
    {
      const double beta_n=beta[point] * normals[point];
					 // We assemble the term
					 // $(\beta\cdot n
					 // u,v)_{\partial\kappa_+}$,
      if (beta_n>0)
	for (unsigned int i=0; i<fe_v.dofs_per_cell; ++i)
	  for (unsigned int j=0; j<fe_v.dofs_per_cell; ++j)
	    ui_vi_matrix(i,j) += beta_n *
			       fe_v.shape_value(j,point) *
			       fe_v.shape_value(i,point) *
			       JxW[point];
      else
					 // and the
					 // term $(\beta\cdot n
					 // \hat u,v)_{\partial
					 // \kappa_-}$.
	for (unsigned int i=0; i<fe_v.dofs_per_cell; ++i)
	  for (unsigned int k=0; k<fe_v_neighbor.dofs_per_cell; ++k)
	    ue_vi_matrix(i,k) += beta_n *
				fe_v_neighbor.shape_value(k,point) *
				fe_v.shape_value(i,point) *
				JxW[point];
    }
}


				 // @sect4{Function: assemble_face_term2}
				 //
				 // Now we look at the
				 // <code>assemble_face_term2</code> function
				 // that assembles the face terms
				 // corresponding to the second
				 // version of the DG method,
				 // cf. above. For that case the face
				 // terms are given as a sum of
				 // integrals over all faces.  Here we
				 // need two additional cell matrices
				 // <code>ui_ve_matrix</code> and
				 // <code>ue_ve_matrix</code> that will store
				 // contributions due to terms
				 // involving ui and ve as well as ue
				 // and ve.
template <int dim>
void DGTransportEquation<dim>::assemble_face_term2(
  const FEFaceValuesBase<dim>& fe_v,
  const FEFaceValuesBase<dim>& fe_v_neighbor,      
  FullMatrix<double> &ui_vi_matrix,
  FullMatrix<double> &ue_vi_matrix,
  FullMatrix<double> &ui_ve_matrix,
  FullMatrix<double> &ue_ve_matrix) const
{
				   // the first few lines are the same
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
	      ui_vi_matrix(i,j) += beta_n *
				 fe_v.shape_value(j,point) *
				 fe_v.shape_value(i,point) *
				 JxW[point];

					   // We additionally assemble
					   // the term $(\beta\cdot n
					   // u,\hat v)_{\partial
					   // \kappa_+}$,
	  for (unsigned int k=0; k<fe_v_neighbor.dofs_per_cell; ++k)
	    for (unsigned int j=0; j<fe_v.dofs_per_cell; ++j)
	      ui_ve_matrix(k,j) -= beta_n *
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
	      ue_vi_matrix(i,l) += beta_n *
				  fe_v_neighbor.shape_value(l,point) *
				  fe_v.shape_value(i,point) *
				  JxW[point];

					   // And this is another new
					   // one: $(\beta\cdot n \hat
					   // u,\hat v)_{\partial
					   // \kappa_-}$.
	  for (unsigned int k=0; k<fe_v_neighbor.dofs_per_cell; ++k)
	    for (unsigned int l=0; l<fe_v_neighbor.dofs_per_cell; ++l)
	      ue_ve_matrix(k,l) -= beta_n *
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
    void assemble_system1 ();
    void assemble_system2 ();
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

				     // We define the quadrature
				     // formulae for the cell and the
				     // face terms of the
				     // discretization.
    const QGauss<dim>   quadrature;
    const QGauss<dim-1> face_quadrature;
    
				     // And there are two solution
				     // vectors, that store the
				     // solutions to the problems
				     // corresponding to the two
				     // different assembling routines
				     // <code>assemble_system1</code> and
				     // <code>assemble_system2</code>;
    Vector<double>       solution1;
    Vector<double>       solution2;
    Vector<double>       right_hand_side;
    
				     // Finally this class includes an
				     // object of the
				     // DGTransportEquations class
				     // described above.
    const DGTransportEquation<dim> dg;
};


template <int dim>
DGMethod<dim>::DGMethod ()
		:
		mapping (),
						 // Change here for DG
						 // methods of
						 // different degrees.
                fe (1),
		dof_handler (triangulation),
		quadrature (4),
		face_quadrature (4),
		dg ()
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

  solution1.reinit (dof_handler.n_dofs());
  solution2.reinit (dof_handler.n_dofs());
  right_hand_side.reinit (dof_handler.n_dofs());
}


				 // @sect4{Function: assemble_system1}
				 //
				 // We proceed with the
				 // <code>assemble_system1</code> function that
				 // implements the DG discretization
				 // in its first version. This
				 // function repeatedly calls the
				 // <code>assemble_cell_term</code>,
				 // <code>assemble_boundary_term</code> and
				 // <code>assemble_face_term1</code> functions
				 // of the <code>DGTransportEquation</code>
				 // object.  The
				 // <code>assemble_boundary_term</code> covers
				 // the first case mentioned in the
				 // introduction.
				 //
				 // 1. face is at boundary
				 //
				 // This function takes a
				 // <code>FEFaceValues</code> object as
				 // argument.  In contrast to that
				 // <code>assemble_face_term1</code>
				 // takes two <code>FEFaceValuesBase</code>
				 // objects; one for the shape
				 // functions on the current cell and
				 // the other for shape functions on
				 // the neighboring cell under
				 // consideration. Both objects are
				 // either of class <code>FEFaceValues</code>
				 // or of class <code>FESubfaceValues</code>
				 // (both derived from
				 // <code>FEFaceValuesBase</code>) according to
				 // the remaining cases mentioned
				 // in the introduction:
				 //
				 // 2. neighboring cell is finer
				 // (current cell: <code>FESubfaceValues</code>,
				 // neighboring cell: <code>FEFaceValues</code>);
				 //
				 // 3. neighboring cell is of the same
				 // refinement level (both, current
				 // and neighboring cell:
				 // <code>FEFaceValues</code>);
				 //
				 // 4. neighboring cell is coarser
				 // (current cell: <code>FEFaceValues</code>,
				 // neighboring cell:
				 // <code>FESubfaceValues</code>).
				 //
				 // If we considered globally refined
				 // meshes then only case 3 would
				 // occur. But as we consider also
				 // locally refined meshes we need to
				 // distinguish all four cases making
				 // the following assembling function
				 // a bit longish.
template <int dim>
void DGMethod<dim>::assemble_system1 () 
{
  const unsigned int dofs_per_cell = dof_handler.get_fe().dofs_per_cell;
  std::vector<unsigned int> dofs (dofs_per_cell);
  std::vector<unsigned int> dofs_neighbor (dofs_per_cell);

				   // First we create the
				   // <code>update_flags</code> for the
				   // <code>FEValues</code> and the
				   // <code>FEFaceValues</code> objects.
  const UpdateFlags update_flags = update_values
                                   | update_gradients
                                   | update_quadrature_points
                                   | update_JxW_values;

				   // Note, that on faces we do not
				   // need gradients but we need
				   // normal vectors.
  const UpdateFlags face_update_flags = update_values
                                        | update_quadrature_points
                                        | update_JxW_values
                                        | update_normal_vectors;
  
				   // On the neighboring cell we only
				   // need the shape values. Given a
				   // specific face, the quadrature
				   // points and `JxW values' are the
				   // same as for the current cells,
				   // the normal vectors are known to
				   // be the negative of the normal
				   // vectors of the current cell.
  const UpdateFlags neighbor_face_update_flags = update_values;
   
				   // Then we create the <code>FEValues</code>
				   // object. Note, that since version
				   // 3.2.0 of deal.II the constructor
				   // of this class takes a
				   // <code>Mapping</code> object as first
				   // argument. Although the
				   // constructor without <code>Mapping</code>
				   // argument is still supported it
				   // is recommended to use the new
				   // constructor. This reduces the
				   // effect of `hidden magic' (the
				   // old constructor implicitely
				   // assumes a <code>MappingQ1</code> mapping)
				   // and makes it easier to change
				   // the mapping object later.
  FEValues<dim> fe_v (
    mapping, fe, quadrature, update_flags);
  
				   // Similarly we create the
				   // <code>FEFaceValues</code> and
				   // <code>FESubfaceValues</code> objects for
				   // both, the current and the
				   // neighboring cell. Within the
				   // following nested loop over all
				   // cells and all faces of the cell
				   // they will be reinited to the
				   // current cell and the face (and
				   // subface) number.
  FEFaceValues<dim> fe_v_face (
    mapping, fe, face_quadrature, face_update_flags);
  FESubfaceValues<dim> fe_v_subface (
    mapping, fe, face_quadrature, face_update_flags);
  FEFaceValues<dim> fe_v_face_neighbor (
    mapping, fe, face_quadrature, neighbor_face_update_flags);
  FESubfaceValues<dim> fe_v_subface_neighbor (
    mapping, fe, face_quadrature, neighbor_face_update_flags);

				   // Now we create the cell matrices
				   // and vectors. Here we need two
				   // cell matrices, both for face
				   // terms that include test
				   // functions <code>vi</code> (internal shape
				   // functions, i.e. shape functions
				   // of the current cell). To be more
				   // precise, the first matrix will
				   // include the `ui and vi terms'
				   // and the second will include the
				   // `ue and vi terms'. Here we
				   // recall the convention that `ui'
				   // is the shortcut for $u_h^+$ and
				   // `ue' represents $u_h^-$, see the
				   // introduction.
  FullMatrix<double> ui_vi_matrix (dofs_per_cell, dofs_per_cell);
  FullMatrix<double> ue_vi_matrix (dofs_per_cell, dofs_per_cell);

  Vector<double>  cell_vector (dofs_per_cell);

				   // Furthermore we need some cell
				   // iterators.
  typename DoFHandler<dim>::active_cell_iterator
    cell = dof_handler.begin_active(),
    endc = dof_handler.end();

				   // Now we start the loop over all
				   // active cells.
  for (;cell!=endc; ++cell) 
    {
				       // In the
				       // <code>assemble_face_term1</code>
				       // function contributions to
				       // the cell matrices and the
				       // cell vector are only
				       // ADDED. Therefore on each
				       // cell we need to reset the
				       // <code>ui_vi_matrix</code> and
				       // <code>cell_vector</code> to zero,
				       // before assembling the cell terms.
      ui_vi_matrix = 0;
      cell_vector = 0;
      
				       // Now we reinit the <code>FEValues</code>
				       // object for the current cell
      fe_v.reinit (cell);

				       // and call the function
				       // that assembles the cell
				       // terms. The first argument is
				       // the <code>FEValues</code> that was
				       // previously reinit'ed on the
				       // current cell.
      dg.assemble_cell_term(fe_v,
			    ui_vi_matrix,
			    cell_vector);

				       // As in previous examples the
				       // vector `dofs' includes the
				       // dof_indices.
      cell->get_dof_indices (dofs);

				       // This is the start of the
				       // nested loop over all faces.
      for (unsigned int face_no=0; face_no<GeometryInfo<dim>::faces_per_cell; ++face_no)
	{
					   // First we set the face
					   // iterator
	  typename DoFHandler<dim>::face_iterator face=cell->face(face_no);
	  
					   // and clear the
					   // <code>ue_vi_matrix</code> on each
					   // face.
	  ue_vi_matrix = 0;
		  
					   // Now we distinguish the
					   // four different cases in
					   // the ordering mentioned
					   // above. We start with
					   // faces belonging to the
					   // boundary of the domain.
	  if (face->at_boundary())
	    {
					       // We reinit the
					       // <code>FEFaceValues</code>
					       // object to the
					       // current face
	      fe_v_face.reinit (cell, face_no);

					       // and assemble the
					       // corresponding face
					       // terms.
	      dg.assemble_boundary_term(fe_v_face,
					ui_vi_matrix,
					cell_vector);
	    }
	  else
	    {
					       // Now we are not on
					       // the boundary of the
					       // domain, therefore
					       // there must exist a
					       // neighboring cell.
	      typename DoFHandler<dim>::cell_iterator neighbor=
		cell->neighbor(face_no);;

					       // We proceed with the
					       // second and most
					       // complicated case:
					       // the neighboring cell
					       // is more refined than
					       // the current cell. As
					       // in deal.II
					       // neighboring cells
					       // are restricted to
					       // have a level
					       // difference of not
					       // more than one, the
					       // neighboring cell is
					       // known to be at most
					       // ONCE more refined
					       // than the current
					       // cell. Furthermore
					       // also the face is
					       // more refined,
					       // i.e. it has
					       // children. Here we
					       // note that the
					       // following part of
					       // code will not work
					       // for <code>dim==1</code>.
	      if (face->has_children())
		{
						   // First we store
						   // which number the
						   // current cell has
						   // in the list of
						   // neighbors of the
						   // neighboring
						   // cell. Hence,
						   // neighbor-@>neighbor(neighbor2)
						   // equals the
						   // current cell
						   // <code>cell</code>.
		  const unsigned int neighbor2=
		    cell->neighbor_of_neighbor(face_no);
		  
		  
						   // We loop over
						   // subfaces
		  for (unsigned int subface_no=0;
		       subface_no<face->n_children(); ++subface_no)
		    {
						       // and set the
						       // cell
						       // iterator
						       // <code>neighbor_child</code>
						       // to the cell
						       // placed
						       // `behind' the
						       // current
						       // subface.
		      typename DoFHandler<dim>::active_cell_iterator
                        neighbor_child
                        = cell->neighbor_child_on_subface (face_no, subface_no);
		      
		      Assert (!neighbor_child->has_children(), ExcInternalError());

						       // We need to
						       // reset the
						       // <code>ue_vi_matrix</code>
						       // on each
						       // subface
						       // because on
						       // each subface
						       // the <code>un</code>
						       // belong to
						       // different
						       // neighboring
						       // cells.
		      ue_vi_matrix = 0;
		      
						       // As already
						       // mentioned
						       // above for
						       // the current
						       // case (case
						       // 2) we employ
						       // the
						       // <code>FESubfaceValues</code>
						       // of the
						       // current
						       // cell (here
						       // reinited for
						       // the current
						       // cell, face
						       // and subface)
						       // and we
						       // employ the
						       // FEFaceValues
						       // of the
						       // neighboring
						       // child cell.
		      fe_v_subface.reinit (cell, face_no, subface_no);
		      fe_v_face_neighbor.reinit (neighbor_child, neighbor2);

		      dg.assemble_face_term1(fe_v_subface,
					     fe_v_face_neighbor,
					     ui_vi_matrix,
					     ue_vi_matrix);
		      
						       // Then we get
						       // the dof
						       // indices of
						       // the
						       // neighbor_child
						       // cell
		      neighbor_child->get_dof_indices (dofs_neighbor);
		      						
						       // and
						       // distribute
						       // <code>ue_vi_matrix</code>
						       // to the
						       // system_matrix
		      for (unsigned int i=0; i<dofs_per_cell; ++i)
			for (unsigned int k=0; k<dofs_per_cell; ++k)
			  system_matrix.add(dofs[i], dofs_neighbor[k],
					    ue_vi_matrix(i,k));
		    }
						   // End of <code>if
						   // (face-@>has_children())</code>
		}
	      else
		{
						   // We proceed with
						   // case 3,
						   // i.e. neighboring
						   // cell is of the
						   // same refinement
						   // level as the
						   // current cell.
		  if (neighbor->level() == cell->level()) 
		    {
						       // Like before
						       // we store
						       // which number
						       // the current
						       // cell has in
						       // the list of
						       // neighbors of
						       // the
						       // neighboring
						       // cell.
		      const unsigned int neighbor2=cell->neighbor_of_neighbor(face_no);

						       // We reinit
						       // the
						       // <code>FEFaceValues</code>
						       // of the
						       // current and
						       // neighboring
						       // cell to the
						       // current face
						       // and assemble
						       // the
						       // corresponding
						       // face terms.
		      fe_v_face.reinit (cell, face_no);
		      fe_v_face_neighbor.reinit (neighbor, neighbor2);
		      
		      dg.assemble_face_term1(fe_v_face,
					     fe_v_face_neighbor,
					     ui_vi_matrix,
					     ue_vi_matrix);
						       // End of <code>if
						       // (neighbor-@>level()
						       // ==
						       // cell-@>level())</code>
		    }
		  else
		    {
						       // Finally we
						       // consider
						       // case 4. When
						       // the
						       // neighboring
						       // cell is not
						       // finer and
						       // not of the
						       // same
						       // refinement
						       // level as the
						       // current cell
						       // it must be
						       // coarser.
		      Assert(neighbor->level() < cell->level(), ExcInternalError());

						       // Find out the
						       // how many'th
						       // face_no and
						       // subface_no
						       // the current
						       // face is
						       // w.r.t. the
						       // neighboring
						       // cell.
		      const std::pair<unsigned int, unsigned int> faceno_subfaceno=
			cell->neighbor_of_coarser_neighbor(face_no);
		      const unsigned int neighbor_face_no=faceno_subfaceno.first,
				      neighbor_subface_no=faceno_subfaceno.second;

		      Assert (neighbor->neighbor_child_on_subface (neighbor_face_no,
                                                                   neighbor_subface_no)
                              == cell,
                              ExcInternalError());

						       // Reinit the
						       // appropriate
						       // <code>FEFaceValues</code>
						       // and assemble
						       // the face
						       // terms.
		      fe_v_face.reinit (cell, face_no);
		      fe_v_subface_neighbor.reinit (neighbor, neighbor_face_no,
						    neighbor_subface_no);
		      
		      dg.assemble_face_term1(fe_v_face,
					     fe_v_subface_neighbor,
					     ui_vi_matrix,
					     ue_vi_matrix);
		    }

						   // Now we get the
						   // dof indices of
						   // the
						   // <code>neighbor_child</code>
						   // cell,
		  neighbor->get_dof_indices (dofs_neighbor);
		 						
						   // and distribute the
						   // <code>ue_vi_matrix</code>.
		  for (unsigned int i=0; i<dofs_per_cell; ++i)
		    for (unsigned int k=0; k<dofs_per_cell; ++k)
		      system_matrix.add(dofs[i], dofs_neighbor[k],
					ue_vi_matrix(i,k));
		}
					       // End of <code>face not at boundary</code>:
	    }
					   // End of loop over all faces:
	}
      
				       // Finally we distribute the
				       // <code>ui_vi_matrix</code>
      for (unsigned int i=0; i<dofs_per_cell; ++i)
	for (unsigned int j=0; j<dofs_per_cell; ++j)
	  system_matrix.add(dofs[i], dofs[j], ui_vi_matrix(i,j));
      
				       // and the cell vector.
      for (unsigned int i=0; i<dofs_per_cell; ++i)
	right_hand_side(dofs[i]) += cell_vector(i);
    }
}


				 // @sect4{Function: assemble_system2}
				 //
				 // We proceed with the
				 // <code>assemble_system2</code> function that
				 // implements the DG discretization
				 // in its second version. This
				 // function is very similar to the
				 // <code>assemble_system1</code>
				 // function. Therefore, here we only
				 // discuss the differences between
				 // the two functions. This function
				 // repeatedly calls the
				 // <code>assemble_face_term2</code> function
				 // of the DGTransportEquation object,
				 // that assembles the face terms
				 // written as a sum of integrals over
				 // all faces. Therefore, we need to
				 // make sure that each face is
				 // treated only once. This is achieved
				 // by introducing the rule:
				 // 
				 // a) If the current and the
				 // neighboring cells are of the same
				 // refinement level we access and
				 // treat the face from the cell with
				 // lower index.
				 //
				 // b) If the two cells are of
				 // different refinement levels we
				 // access and treat the face from the
				 // coarser cell.
				 //
				 // Due to rule b) we do not need to
				 // consider case 4 (neighboring cell
				 // is coarser) any more.

template <int dim>
void DGMethod<dim>::assemble_system2 () 
{
  const unsigned int dofs_per_cell = dof_handler.get_fe().dofs_per_cell;
  std::vector<unsigned int> dofs (dofs_per_cell);
  std::vector<unsigned int> dofs_neighbor (dofs_per_cell);

  const UpdateFlags update_flags = update_values
                                   | update_gradients
                                   | update_quadrature_points
                                   | update_JxW_values;
  
  const UpdateFlags face_update_flags = update_values
                                        | update_quadrature_points
                                        | update_JxW_values
                                        | update_normal_vectors;
  
  const UpdateFlags neighbor_face_update_flags = update_values;

				   // Here we do not need
				   // <code>fe_v_face_neighbor</code> as case 4
				   // does not occur.
  FEValues<dim> fe_v (
    mapping, fe, quadrature, update_flags);
  FEFaceValues<dim> fe_v_face (
    mapping, fe, face_quadrature, face_update_flags);
  FESubfaceValues<dim> fe_v_subface (
    mapping, fe, face_quadrature, face_update_flags);
  FEFaceValues<dim> fe_v_face_neighbor (
    mapping, fe, face_quadrature, neighbor_face_update_flags);


  FullMatrix<double> ui_vi_matrix (dofs_per_cell, dofs_per_cell);
  FullMatrix<double> ue_vi_matrix (dofs_per_cell, dofs_per_cell);
  
				   // Additionally we need the
				   // following two cell matrices,
				   // both for face term that include
				   // test function <code>ve</code> (external
				   // shape functions, i.e. shape
				   // functions of the neighboring
				   // cell). To be more precise, the
				   // first matrix will include the `u
				   // and vn terms' and the second
				   // that will include the `un and vn
				   // terms'.
  FullMatrix<double> ui_ve_matrix (dofs_per_cell, dofs_per_cell);
  FullMatrix<double> ue_ve_matrix (dofs_per_cell, dofs_per_cell);
  
  Vector<double>  cell_vector (dofs_per_cell);

				   // The following lines are roughly
				   // the same as in the previous
				   // function.
  typename DoFHandler<dim>::active_cell_iterator
    cell = dof_handler.begin_active(),
    endc = dof_handler.end();
  for (;cell!=endc; ++cell) 
    {
      ui_vi_matrix = 0;
      cell_vector = 0;

      fe_v.reinit (cell);

      dg.assemble_cell_term(fe_v,
			    ui_vi_matrix,
			    cell_vector);
      
      cell->get_dof_indices (dofs);

      for (unsigned int face_no=0; face_no<GeometryInfo<dim>::faces_per_cell; ++face_no)
	{
	  typename DoFHandler<dim>::face_iterator face=
	    cell->face(face_no);

					   // Case 1:
	  if (face->at_boundary())
	    {
	      fe_v_face.reinit (cell, face_no);

	      dg.assemble_boundary_term(fe_v_face,
					ui_vi_matrix,
					cell_vector);
	    }
	  else
	    {
	      Assert (cell->neighbor(face_no).state() == IteratorState::valid,
		      ExcInternalError());
	      typename DoFHandler<dim>::cell_iterator neighbor=
		cell->neighbor(face_no);
					       // Case 2:
	      if (face->has_children())
		{
		  const unsigned int neighbor2=
		    cell->neighbor_of_neighbor(face_no);
		  
		  for (unsigned int subface_no=0;
		       subface_no<face->n_children(); ++subface_no)
		    {
		      typename DoFHandler<dim>::cell_iterator neighbor_child
                        = cell->neighbor_child_on_subface (face_no, subface_no);
		      Assert (!neighbor_child->has_children(), ExcInternalError());
		      
		      ue_vi_matrix = 0;
		      ui_ve_matrix = 0;
		      ue_ve_matrix = 0;
		      
		      fe_v_subface.reinit (cell, face_no, subface_no);
		      fe_v_face_neighbor.reinit (neighbor_child, neighbor2);

		      dg.assemble_face_term2(fe_v_subface,
					     fe_v_face_neighbor,
					     ui_vi_matrix,
					     ue_vi_matrix,
					     ui_ve_matrix,
					     ue_ve_matrix);
		  
		      neighbor_child->get_dof_indices (dofs_neighbor);
		      						
		      for (unsigned int i=0; i<dofs_per_cell; ++i)
			for (unsigned int j=0; j<dofs_per_cell; ++j)
			  {
			    system_matrix.add(dofs[i], dofs_neighbor[j],
					      ue_vi_matrix(i,j));
			    system_matrix.add(dofs_neighbor[i], dofs[j],
					      ui_ve_matrix(i,j));
			    system_matrix.add(dofs_neighbor[i], dofs_neighbor[j],
					      ue_ve_matrix(i,j));
			  }
		    }
		}
	      else
		{
						   // Case 3, with the
						   // additional rule
						   // a)
		  if (neighbor->level() == cell->level() &&
		      neighbor->index() > cell->index()) 
		    {
		      const unsigned int neighbor2=cell->neighbor_of_neighbor(face_no);
		      
		      ue_vi_matrix = 0;
		      ui_ve_matrix = 0;
		      ue_ve_matrix = 0;
		      
		      fe_v_face.reinit (cell, face_no);
		      fe_v_face_neighbor.reinit (neighbor, neighbor2);
		      
		      dg.assemble_face_term2(fe_v_face,
					     fe_v_face_neighbor,
					     ui_vi_matrix,
					     ue_vi_matrix,
					     ui_ve_matrix,
					     ue_ve_matrix);

		      neighbor->get_dof_indices (dofs_neighbor);

		      for (unsigned int i=0; i<dofs_per_cell; ++i)
			for (unsigned int j=0; j<dofs_per_cell; ++j)
			  {
			    system_matrix.add(dofs[i], dofs_neighbor[j],
					      ue_vi_matrix(i,j));
			    system_matrix.add(dofs_neighbor[i], dofs[j],
					      ui_ve_matrix(i,j));
			    system_matrix.add(dofs_neighbor[i], dofs_neighbor[j],
					      ue_ve_matrix(i,j));
			  }
		    }

						   // Due to rule b)
						   // we do not need
						   // to consider case
						   // 4.
		}
	    }
	}
      
      for (unsigned int i=0; i<dofs_per_cell; ++i)
	for (unsigned int j=0; j<dofs_per_cell; ++j)
	  system_matrix.add(dofs[i], dofs[j], ui_vi_matrix(i,j));
      
      for (unsigned int i=0; i<dofs_per_cell; ++i)
	right_hand_side(dofs[i]) += cell_vector(i);
    }
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
				 // neighboring cells. The simpliest
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
				 // u_h|$. Futhermore we note that we
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
						 solution2,
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
  data_out.add_data_vector (solution2, "u");

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

	  triangulation.refine_global (5);
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

				       // The constructor of the Timer
				       // class automatically starts
				       // the time measurement.
      Timer assemble_timer;
				       // First assembling routine.
      assemble_system1 ();
				       // The operator () accesses the
				       // current time without
				       // disturbing the time
				       // measurement.
      std::cout << "Time of assemble_system1: "
		<< assemble_timer()
		<< std::endl;
      solve (solution1);

				       // As preparation for the
				       // second assembling routine we
				       // reinit the system matrix, the
				       // right hand side vector and
				       // the Timer object.
      system_matrix = 0;
      right_hand_side = 0;
      assemble_timer.reset();

				       // We start the Timer,
      assemble_timer.start();
				       // call the second assembling routine
      assemble_system2 ();
				       // and access the current time.
      std::cout << "Time of assemble_system2: "
		<< assemble_timer()
		<< std::endl;
      solve (solution2);

				       // To make sure that both
				       // versions of the DG method
				       // yield the same
				       // discretization and hence the
				       // same solution we check the
				       // two solutions for equality.
      solution1-=solution2;
      const double difference=solution1.linfty_norm();
      if (difference>1e-12)
	std::cout << "solution1 and solution2 differ!!" << std::endl;
      else
	std::cout << "solution1 and solution2 coincide." << std::endl;
	
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


