/* $Id$ */
/* Author: Ralf Hartmann, University of Heidelberg, 2001 */

/*    $Id$       */
/*    Version: $Name$                                          */
/*                                                                */
/*    Copyright (C) 2001, 2002 by the deal.II authors */
/*                                                                */
/*    This file is subject to QPL and may not be  distributed     */
/*    without copyright and license information. Please refer     */
/*    to the file deal.II/doc/license.html for the  text  and     */
/*    further information on this license.                        */

				 // The first few files have already
				 // been covered in previous examples
				 // and will thus not be further
				 // commented on.
#include <base/quadrature_lib.h>
#include <base/function.h>
#include <lac/vector.h>
#include <lac/sparse_matrix.h>
#include <lac/vector_memory.h>
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

				 // This is the first new file. It
				 // declares the ``MappingQ1'' class
				 // that gives the standard bilinear
				 // mapping. For bilinear mappings use
				 // an object of this class rather
				 // than an object of the
				 // ``MappingQ(1)'' class, as the
				 // ``MappingQ1'' class is optimized
				 // due to the pre-knowledge of the
				 // actual polynomial degree 1.
#include <fe/mapping_q1.h>
				 // Here the discontinuous finite
				 // elements are defined. They are
				 // used as all other finite elements.
#include <fe/fe_dgq.h>
				 // We are going to use the simplest
				 // possible solver, called Richardson
				 // iteration, that represents a simple
				 // defect correction. This, in
				 // combination with a block SSOR
				 // preconditioner (defined in
				 // precondition_block.h), that uses
				 // the special block matrix structur
				 // of system matrices arising from DG
				 // discretizations.
#include <lac/solver_richardson.h>
#include <lac/precondition_block.h>
				 // We are going to use gradients as
				 // refinement indicator.
#include <numerics/derivative_approximation.h>
				 // Finally we do some time comparison
				 // using the ``Timer'' class.
#include <base/timer.h>

				 // And this again is C++:
#include <iostream>
#include <fstream>


				 // @sect3{Equation data}
				 //
				 // First we define the classes
				 // representing the equation-specific
				 // functions. Both classes, ``RHS''
				 // and ``BoundaryValues'', are
				 // derived from the ``Function''
				 // class. Only the ``value_list''
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


				 // The class ``Beta'' represents the
				 // vector valued flow field of the
				 // linear transport equation and is
				 // not derived from the ``Function''
				 // class as we prefer to get function
				 // values of type ``Point'' rather
				 // than of type
				 // ``Vector<double>''. This, because
				 // there exist scalar products
				 // between ``Point'' and ``Point'' as
				 // well as between ``Point'' and
				 // ``Tensor'', simplifying terms like
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
    Beta () {};
    void value_list (const std::vector<Point<dim> > &points,
		     std::vector<Point<dim> > &values) const;
};


				 // The implementation of the
				 // ``value_list'' functions of these
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
      const Point<dim> &p=points[i];
      Point<dim> &beta=values[i];

      beta(0) = -p(1);
      beta(1) = p(0);
      beta /= std::sqrt(beta.square());
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
				 // ``DGTransportEquation''. Its
				 // member functions were already
				 // mentioned in the Introduction and
				 // will be explained
				 // below. Furthermore it includes
				 // objects of the previously defined
				 // ``Beta'', ``RHS'' and
				 // ``BoundaryValues'' function
				 // classes.
template <int dim>
class DGTransportEquation
{
  public:
    DGTransportEquation() {};

    void assemble_cell_term(const FEValues<dim>& fe_v,
			    FullMatrix<double> &u_v_matrix,
			    Vector<double> &cell_vector) const;
    
    void assemble_boundary_term(const FEFaceValues<dim>& fe_v,
				FullMatrix<double> &u_v_matrix,
				Vector<double> &cell_vector) const;
    
    void assemble_face_term1(const FEFaceValuesBase<dim>& fe_v,
			     const FEFaceValuesBase<dim>& fe_v_neighbor,
			     FullMatrix<double> &u_v_matrix,
			     FullMatrix<double> &un_v_matrix) const;

    void assemble_face_term2(const FEFaceValuesBase<dim>& fe_v,
			     const FEFaceValuesBase<dim>& fe_v_neighbor,
			     FullMatrix<double> &u_v_matrix,
			     FullMatrix<double> &un_v_matrix,
			     FullMatrix<double> &u_vn_matrix,
			     FullMatrix<double> &un_vn_matrix) const;
  private:
    const Beta<dim> beta_function;
    const RHS<dim> rhs_function;
    const BoundaryValues<dim> boundary_function;
};


				 // @sect4{Function: assemble_cell_term}
				 //
				 // The ``assemble_cell_term''
				 // function assembles the cell terms
				 // of the discretization.
				 // ``u_v_matrix'' is a cell matrix,
				 // i.e. for a DG method of degree 1,
				 // it is of size 4 times 4, and
				 // ``cell_vector'' is of size 4.
				 // When this function is invoked,
				 // ``fe_v'' is already reinit'ed with the
				 // current cell before and includes
				 // all shape values needed.
template <int dim>
void DGTransportEquation<dim>::assemble_cell_term(
  const FEValues<dim> &fe_v,
  FullMatrix<double> &u_v_matrix,
  Vector<double> &cell_vector) const
{
				   // First we ask ``fe_v'' for the
				   // shape gradients, shape values and
				   // quadrature weights,
  const std::vector<std::vector<Tensor<1,dim> > > &grad_v = fe_v.get_shape_grads ();
  const FullMatrix<double> &v = fe_v.get_shape_values ();
  const std::vector<double> &JxW = fe_v.get_JxW_values ();

				   // Then the flow field beta and the
				   // ``rhs_function'' are evaluated at
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
	  u_v_matrix(i,j) -= beta[point]*grad_v[i][point]*
			      v[j][point] *
			      JxW[point];
	
	cell_vector(i) += rhs[point] * v[i][point] * JxW[point];
      }
}


				 // @sect4{Function: assemble_boundary_term}
				 //
				 // The ``assemble_boundary_term''
				 // function assembles the face terms
				 // at boundary faces.  When this
				 // function is invoked, ``fe_v'' is
				 // already reinit'ed with the current
				 // cell and current face. Hence it
				 // provides the shape values on that
				 // boundary face.
template <int dim>
void DGTransportEquation<dim>::assemble_boundary_term(
  const FEFaceValues<dim>& fe_v,    
  FullMatrix<double> &u_v_matrix,
  Vector<double> &cell_vector) const
{
				   // First we check whether the
				   // current face is really at the
				   // boundary.
  Assert(fe_v.get_face()->at_boundary(), ExcInternalError());
  
				   // Again, as in the previous
				   // function, we ask the ``FEValues''
				   // object for the shape values and
				   // the quadrature weights
  const FullMatrix<double> &v = fe_v.get_shape_values ();
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
					 // u,v)_{\partial K_+}$,
      if (beta_n>0)
	for (unsigned int i=0; i<fe_v.dofs_per_cell; ++i)
	  for (unsigned int j=0; j<fe_v.dofs_per_cell; ++j)
	    u_v_matrix(i,j) += beta_n *
			       v[j][point] *
			       v[i][point] *
			       JxW[point];
      else
					 // and the term $(\beta\cdot
					 // n g,v)_{\partial
					 // K_-\cap\partial\Omega}$,
	for (unsigned int i=0; i<fe_v.dofs_per_cell; ++i)
	  cell_vector(i) -= beta_n *
			    g[point] *
			    v[i][point] *
			    JxW[point];
    }
}


				 // @sect4{Function: assemble_face_term1}
				 //
				 // The ``assemble_face_term1''
				 // function assembles the face terms
				 // corresponding to the first version
				 // of the DG method, cf. above. For
				 // that case, the face terms are
				 // given as a sum of integrals over
				 // all cell boundaries.
				 //
				 // When this function is invoked,
				 // ``fe_v'' and ``fe_v_neighbor'' are
				 // already reinit'ed with the current
				 // cell and the neighoring cell,
				 // respectively, as well as with the
				 // current face. Hence they provide
				 // the inner and outer shape values
				 // on the face.
				 //
				 // In addition to the cell matrix
				 // ``u_v_matrix'' this function has
				 // got a new argument
				 // ``un_v_matrix'', that stores
				 // contributions to the system matrix
				 // that are based on outer values of
				 // u, see $\hat u_h$ in the
				 // introduction, and inner values of
				 // v, see $v_h$. Here we note that
				 // ``un'' is the short notation for
				 // ``u_neighbor'' and represents
				 // $\hat u_h$.
template <int dim>
void DGTransportEquation<dim>::assemble_face_term1(
  const FEFaceValuesBase<dim>& fe_v,
  const FEFaceValuesBase<dim>& fe_v_neighbor,      
  FullMatrix<double> &u_v_matrix,
  FullMatrix<double> &un_v_matrix) const
{
				   // First we check that the current
				   // face is not at the boundary by
				   // accident.
  Assert(!fe_v.get_face()->at_boundary(), ExcInternalError());
  
				   // Again, as in the previous
				   // function, we ask the FEValues
				   // objects for the shape values,
				   // the quadrature weights and the
				   // normals
  const FullMatrix<double> &v = fe_v.get_shape_values ();
  const FullMatrix<double> &v_neighbor = fe_v_neighbor.get_shape_values ();  
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
					 // u,v)_{\partial K_+}$,
      if (beta_n>0)
	for (unsigned int i=0; i<fe_v.dofs_per_cell; ++i)
	  for (unsigned int j=0; j<fe_v.dofs_per_cell; ++j)
	    u_v_matrix(i,j) += beta_n *
			       v[j][point] *
			       v[i][point] *
			       JxW[point];
      else
					 // and the
					 // term $(\beta\cdot n
					 // \hat u,v)_{\partial
					 // K_-}$.
	for (unsigned int i=0; i<fe_v.dofs_per_cell; ++i)
	  for (unsigned int k=0; k<fe_v_neighbor.dofs_per_cell; ++k)
	    un_v_matrix(i,k) += beta_n *
				v_neighbor[k][point] *
				v[i][point] *
				JxW[point];
    }
}


				 // @sect4{Function: assemble_face_term2}
				 //
				 // Now we look at the
				 // ``assemble_face_term2'' function
				 // that assembles the face terms
				 // corresponding to the second
				 // version of the DG method,
				 // cf. above. For that case the face
				 // terms are given as a sum of
				 // integrals over all faces.  Here we
				 // need two additional cell matrices
				 // ``u_vn_matrix'' and
				 // ``un_vn_matrix'' that will store
				 // contributions due to terms
				 // involving u and vn as well as un
				 // and vn.
template <int dim>
void DGTransportEquation<dim>::assemble_face_term2(
  const FEFaceValuesBase<dim>& fe_v,
  const FEFaceValuesBase<dim>& fe_v_neighbor,      
  FullMatrix<double> &u_v_matrix,
  FullMatrix<double> &un_v_matrix,
  FullMatrix<double> &u_vn_matrix,
  FullMatrix<double> &un_vn_matrix) const
{
				   // the first few lines are the same
  const FullMatrix<double> &v = fe_v.get_shape_values ();
  const FullMatrix<double> &v_neighbor = fe_v_neighbor.get_shape_values ();  
  const std::vector<double> &JxW = fe_v.get_JxW_values ();
  const std::vector<Point<dim> > &normals = fe_v.get_normal_vectors ();

  std::vector<Point<dim> > beta (fe_v.n_quadrature_points);
  
  beta_function.value_list (fe_v.get_quadrature_points(), beta);

  for (unsigned int point=0; point<fe_v.n_quadrature_points; ++point)
    {
      const double beta_n=beta[point] * normals[point];
      if (beta_n>0)
	{
					   // This terms we've already seen.
	  for (unsigned int i=0; i<fe_v.dofs_per_cell; ++i)
	    for (unsigned int j=0; j<fe_v.dofs_per_cell; ++j)
	      u_v_matrix(i,j) += beta_n *
				 v[j][point] *
				 v[i][point] *
				 JxW[point];

					   // We additionally assemble
					   // the term $(\beta\cdot n
					   // u,\hat v)_{\partial
					   // K_+},
	  for (unsigned int k=0; k<fe_v_neighbor.dofs_per_cell; ++k)
	    for (unsigned int j=0; j<fe_v.dofs_per_cell; ++j)
	      u_vn_matrix(k,j) -= beta_n *
				  v[j][point] *
				  v_neighbor[k][point] *
				  JxW[point];
	}
      else
	{
					   // This one we've already
					   // seen, too.
	  for (unsigned int i=0; i<fe_v.dofs_per_cell; ++i)
	    for (unsigned int l=0; l<fe_v_neighbor.dofs_per_cell; ++l)
	      un_v_matrix(i,l) += beta_n *
				  v_neighbor[l][point] *
				  v[i][point] *
				  JxW[point];

					   // And this is another new
					   // one: $(\beta\cdot n \hat
					   // u,\hat v)_{\partial
					   // K_-}$.
	  for (unsigned int k=0; k<fe_v_neighbor.dofs_per_cell; ++k)
	    for (unsigned int l=0; l<fe_v_neighbor.dofs_per_cell; ++l)
	      un_vn_matrix(k,l) -= beta_n *
				   v_neighbor[l][point] *
				   v_neighbor[k][point] *
				   JxW[point];
	}
    }
}


				 // @sect3{Class: DGMethod}
				 //
				 // After these preparations, we
				 // proceed with the main part of this
				 // program. The main class, here
				 // called ``DGMethod'' is basically
				 // the main class of step 6. One of
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
    const QGauss4<dim>   quadrature;
    const QGauss4<dim-1> face_quadrature;
    
				     // And there are two solution
				     // vectors, that store the
				     // solutions to the problems
				     // corresponding to the two
				     // different assembling routines
				     // ``assemble_system1'' and
				     // ``assemble_system2'';
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
DGMethod<dim>::DGMethod () :
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
};


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
			    *GeometryInfo<dim>::subfaces_per_face+1)*fe.dofs_per_cell);
  
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
};


				 // @sect4{Function: assemble_system1}
				 //
				 // We proceed with the
				 // ``assemble_system1'' function that
				 // implements the DG discretization
				 // in its first version. This
				 // function repeatedly calls the
				 // ``assemble_cell_term'',
				 // ``assemble_boundary_term'' and
				 // ``assemble_face_term1'' functions
				 // of the ``DGTransportEquation''
				 // object.  The
				 // ``assemble_boundary_term'' covers
				 // the first case mentioned in the
				 // introduction.
				 //
				 // 1. face is at boundary
				 //
				 // This function takes a
				 // ``FEFaceValues'' object as
				 // argument.  In contrast to that
				 // ``assemble_face_term1''
				 // takes two ``FEFaceValuesBase''
				 // objects; one for the shape
				 // functions on the current cell and
				 // the other for shape functions on
				 // the neighboring cell under
				 // consideration. Both objects are
				 // either of class ``FEFaceValues''
				 // or of class ``FESubfaceValues''
				 // (both derived from
				 // ``FEFaceValuesBase'') according to
				 // the remaining cases mentioned
				 // in the introduction:
				 //
				 // 2. neighboring cell is finer
				 // (current cell: ``FESubfaceValues'',
				 // neighboring cell: ``FEFaceValues'');
				 //
				 // 3. neighboring cell is of the same
				 // refinement level (both, current
				 // and neighboring cell:
				 // ``FEFaceValues'');
				 //
				 // 4. neighboring cell is coarser
				 // (current cell: ``FEFaceValues'',
				 // neighboring cell:
				 // ``FESubfaceValues'').
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
				   // ``UpdateFlags'' for the
				   // ``FEValues'' and the
				   // ``FEFaceValues'' objects.
  UpdateFlags update_flags = update_values
			     | update_gradients
			     | update_q_points
			     | update_JxW_values;

				   // Note, that on faces we do not
				   // need gradients but we need
				   // normal vectors.
  UpdateFlags face_update_flags = update_values
				  | update_q_points
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
  UpdateFlags neighbor_face_update_flags = update_values;
   
				   // Then we create the ``FEValues''
				   // object. Note, that since version
				   // 3.2.0 of deal.II the constructor
				   // of this class takes a
				   // ``Mapping'' object as first
				   // argument. Although the
				   // constructor without ``Mapping''
				   // argument is still supported it
				   // is recommended to use the new
				   // constructor. This reduces the
				   // effect of `hidden magic' (the
				   // old constructor implicitely
				   // assumes a ``MappingQ1'' mapping)
				   // and makes it easier to change
				   // the mapping object later.
  FEValues<dim> fe_v (
    mapping, fe, quadrature, update_flags);
  
				   // Similarly we create the
				   // ``FEFaceValues'' and
				   // ``FESubfaceValues'' objects for
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
				   // functions ``v'' (shape functions
				   // of the current cell). To be more
				   // precise, the first matrix will
				   // include the `u and v terms' and
				   // the second that will include the
				   // `un and v terms'. Here we recall
				   // the convention that `un' is
				   // the shortcut for `u_neighbor'
				   // and represents the $u_hat$, see
				   // introduction.
  FullMatrix<double> u_v_matrix (dofs_per_cell, dofs_per_cell);
  FullMatrix<double> un_v_matrix (dofs_per_cell, dofs_per_cell);

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
				       // ``assemble_face_term1''
				       // function contributions to
				       // the cell matrices and the
				       // cell vector are only
				       // ADDED. Therefore on each
				       // cell we need to reset the
				       // ``u_v_matrix'' and
				       // ``cell_vector'' to zero,
				       // before assembling the cell terms.
      u_v_matrix.clear ();
      cell_vector.clear ();
      
				       // Now we reinit the ``FEValues''
				       // object for the current cell
      fe_v.reinit (cell);

				       // and call the function
				       // that assembles the cell
				       // terms. The first argument is
				       // the ``FEValues'' that was
				       // previously reinit'ed on the
				       // current cell.
      dg.assemble_cell_term(fe_v,
			    u_v_matrix,
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
					   // ``un_v_matrix'' on each
					   // face.
	  un_v_matrix.clear();
		  
					   // Now we distinguish the
					   // four different cases in
					   // the ordering mentioned
					   // above. We start with
					   // faces belonging to the
					   // boundary of the domain.
	  if (face->at_boundary())
	    {
					       // We reinit the
					       // ``FEFaceValues''
					       // object to the
					       // current face
	      fe_v_face.reinit (cell, face_no);

					       // and assemble the
					       // corresponding face
					       // terms.
	      dg.assemble_boundary_term(fe_v_face,
					u_v_matrix,
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
					       // for ``dim==1''.
	      if (face->has_children())
		{
						   // First we store
						   // which number the
						   // current cell has
						   // in the list of
						   // neighbors of the
						   // neighboring
						   // cell. Hence,
						   // neighbor->neighbor(neighbor2)
						   // equals the
						   // current cell
						   // ``cell''.
		  const unsigned int neighbor2=
		    cell->neighbor_of_neighbor(face_no);
		  
		  
						   // We loop over
						   // subfaces
		  for (unsigned int subface_no=0;
		       subface_no<GeometryInfo<dim>::subfaces_per_face;
		       ++subface_no)
		    {
						       // and set the
						       // cell
						       // iterator
						       // ``neighbor_child''
						       // to the cell
						       // placed
						       // `behind' the
						       // current
						       // subface.
		      typename DoFHandler<dim>::active_cell_iterator neighbor_child=
			neighbor->child(GeometryInfo<dim>::
					child_cell_on_face(neighbor2,subface_no));
		      
						       // As these are
						       // quite
						       // complicated
						       // indirections
						       // which one
						       // does not
						       // usually get
						       // right at
						       // first
						       // attempt we
						       // check for
						       // the internal
						       // consistency.
		      Assert (neighbor_child->face(neighbor2) == face->child(subface_no),
			      ExcInternalError());
		      Assert (!neighbor_child->has_children(), ExcInternalError());

						       // We need to
						       // reset the
						       // ``un_v_matrix''
						       // on each
						       // subface
						       // because on
						       // each subface
						       // the ``un''
						       // belong to
						       // different
						       // neighboring
						       // cells.
		      un_v_matrix.clear();
		      
						       // As already
						       // mentioned
						       // above for
						       // the current
						       // case (case
						       // 2) we employ
						       // the
						       // ``FESubfaceValues''
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
					     u_v_matrix,
					     un_v_matrix);
		      
						       // Then we get
						       // the dof
						       // indices of
						       // the
						       // neighbor_child
						       // cell
		      neighbor_child->get_dof_indices (dofs_neighbor);
		      						
						       // and
						       // distribute
						       // ``un_v_matrix''
						       // to the
						       // system_matrix
		      for (unsigned int i=0; i<dofs_per_cell; ++i)
			for (unsigned int k=0; k<dofs_per_cell; ++k)
			  system_matrix.add(dofs[i], dofs_neighbor[k],
					    un_v_matrix(i,k));
		    }
						   // End of ``if
						   // (face->has_children())''
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
						       // ``FEFaceValues''
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
					     u_v_matrix,
					     un_v_matrix);
						       // End of ``if
						       // (neighbor->level()
						       // ==
						       // cell->level())''
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

		      Assert (neighbor->neighbor(neighbor_face_no)
			      ->child(GeometryInfo<dim>::child_cell_on_face(
				face_no,neighbor_subface_no)) == cell, ExcInternalError());

						       // Reinit the
						       // appropriate
						       // ``FEFaceValues''
						       // and assemble
						       // the face
						       // terms.
		      fe_v_face.reinit (cell, face_no);
		      fe_v_subface_neighbor.reinit (neighbor, neighbor_face_no,
						    neighbor_subface_no);
		      
		      dg.assemble_face_term1(fe_v_face,
					     fe_v_subface_neighbor,
					     u_v_matrix,
					     un_v_matrix);
		    }

						   // Now we get the
						   // dof indices of
						   // the
						   // ``neighbor_child''
						   // cell,
		  neighbor->get_dof_indices (dofs_neighbor);
		 						
						   // and distribute the
						   // ``un_v_matrix''.
		  for (unsigned int i=0; i<dofs_per_cell; ++i)
		    for (unsigned int k=0; k<dofs_per_cell; ++k)
		      system_matrix.add(dofs[i], dofs_neighbor[k],
					un_v_matrix(i,k));
		}
					       // End of ``face not at boundary'':
	    }
					   // End of loop over all faces:
	}
      
				       // Finally we distribute the
				       // ``u_v_matrix''
      for (unsigned int i=0; i<dofs_per_cell; ++i)
	for (unsigned int j=0; j<dofs_per_cell; ++j)
	  system_matrix.add(dofs[i], dofs[j], u_v_matrix(i,j));
      
				       // and the cell vector.
      for (unsigned int i=0; i<dofs_per_cell; ++i)
	right_hand_side(dofs[i]) += cell_vector(i);
    }
};


				 // @sect4{Function: assemble_system2}
				 //
				 // We proceed with the
				 // ``assemble_system2'' function that
				 // implements the DG discretization
				 // in its second version. This
				 // function is very similar to the
				 // ``assemble_system1''
				 // function. Therefore, here we only
				 // discuss the differences between
				 // the two functions. This function
				 // repeatedly calls the
				 // ``assemble_face_term2'' function
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

  UpdateFlags update_flags = update_values
			     | update_gradients
			     | update_q_points
			     | update_JxW_values;
  
  UpdateFlags face_update_flags = update_values
				  | update_q_points
				  | update_JxW_values
				  | update_normal_vectors;
  
  UpdateFlags neighbor_face_update_flags = update_values;

				   // Here we do not need
				   // ``fe_v_face_neighbor'' as case 4
				   // does not occur.
  FEValues<dim> fe_v (
    mapping, fe, quadrature, update_flags);
  FEFaceValues<dim> fe_v_face (
    mapping, fe, face_quadrature, face_update_flags);
  FESubfaceValues<dim> fe_v_subface (
    mapping, fe, face_quadrature, face_update_flags);
  FEFaceValues<dim> fe_v_face_neighbor (
    mapping, fe, face_quadrature, neighbor_face_update_flags);


  FullMatrix<double> u_v_matrix (dofs_per_cell, dofs_per_cell);
  FullMatrix<double> un_v_matrix (dofs_per_cell, dofs_per_cell);
  
				   // Additionally we need the
				   // following two cell matrices,
				   // both for face term that include
				   // test function ``vn'' (shape
				   // functions of the neighboring
				   // cell). To be more precise, the
				   // first matrix will include the `u
				   // and vn terms' and the second
				   // that will include the `un and vn
				   // terms'.
  FullMatrix<double> u_vn_matrix (dofs_per_cell, dofs_per_cell);
  FullMatrix<double> un_vn_matrix (dofs_per_cell, dofs_per_cell);
  
  Vector<double>  cell_vector (dofs_per_cell);

				   // The following lines are roughly
				   // the same as in the previous
				   // function.
  typename DoFHandler<dim>::active_cell_iterator
    cell = dof_handler.begin_active(),
    endc = dof_handler.end();
  for (;cell!=endc; ++cell) 
    {
      u_v_matrix.clear ();
      cell_vector.clear ();

      fe_v.reinit (cell);

      dg.assemble_cell_term(fe_v,
			    u_v_matrix,
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
					u_v_matrix,
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
		       subface_no<GeometryInfo<dim>::subfaces_per_face;
		       ++subface_no)
		    {
		      typename DoFHandler<dim>::cell_iterator neighbor_child=
			neighbor->child(GeometryInfo<dim>::child_cell_on_face(
			  neighbor2,subface_no));
		      Assert (neighbor_child->face(neighbor2) == face->child(subface_no),
			      ExcInternalError());
		      Assert (!neighbor_child->has_children(), ExcInternalError());
		      
		      un_v_matrix.clear();
		      u_vn_matrix.clear();
		      un_vn_matrix.clear();
		      
		      fe_v_subface.reinit (cell, face_no, subface_no);
		      fe_v_face_neighbor.reinit (neighbor_child, neighbor2);

		      dg.assemble_face_term2(fe_v_subface,
					     fe_v_face_neighbor,
					     u_v_matrix,
					     un_v_matrix,
					     u_vn_matrix,
					     un_vn_matrix);
		  
		      neighbor_child->get_dof_indices (dofs_neighbor);
		      						
		      for (unsigned int i=0; i<dofs_per_cell; ++i)
			for (unsigned int j=0; j<dofs_per_cell; ++j)
			  {
			    system_matrix.add(dofs[i], dofs_neighbor[j],
					      un_v_matrix(i,j));
			    system_matrix.add(dofs_neighbor[i], dofs[j],
					      u_vn_matrix(i,j));
			    system_matrix.add(dofs_neighbor[i], dofs_neighbor[j],
					      un_vn_matrix(i,j));
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
		      
		      un_v_matrix.clear();
		      u_vn_matrix.clear();
		      un_vn_matrix.clear();
		      
		      fe_v_face.reinit (cell, face_no);
		      fe_v_face_neighbor.reinit (neighbor, neighbor2);
		      
		      dg.assemble_face_term2(fe_v_face,
					     fe_v_face_neighbor,
					     u_v_matrix,
					     un_v_matrix,
					     u_vn_matrix,
					     un_vn_matrix);

		      neighbor->get_dof_indices (dofs_neighbor);

		      for (unsigned int i=0; i<dofs_per_cell; ++i)
			for (unsigned int j=0; j<dofs_per_cell; ++j)
			  {
			    system_matrix.add(dofs[i], dofs_neighbor[j],
					      un_v_matrix(i,j));
			    system_matrix.add(dofs_neighbor[i], dofs[j],
					      u_vn_matrix(i,j));
			    system_matrix.add(dofs_neighbor[i], dofs_neighbor[j],
					      un_vn_matrix(i,j));
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
	  system_matrix.add(dofs[i], dofs[j], u_v_matrix(i,j));
      
      for (unsigned int i=0; i<dofs_per_cell; ++i)
	right_hand_side(dofs[i]) += cell_vector(i);
    }
};


				 // @sect3{All the rest}
				 //
				 // For this simple problem we use the
				 // simplest possible solver, called
				 // Richardson iteration, that
				 // represents a simple defect
				 // correction. This, in combination
				 // with a block SSOR preconditioner,
				 // that uses the special block matrix
				 // structur of system matrices
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
  PrimitiveVectorMemory<> vector_memory;
  SolverRichardson<>      solver (solver_control, vector_memory);

				   // Here we create the
				   // preconditioner,
  PreconditionBlockSSOR<double> preconditioner;

				   // we asigned the matrix to it and
				   // set the right block size.
  preconditioner.initialize(system_matrix, fe.dofs_per_cell);

				   // As the inverses of the diagonal
				   // blocks are needed in each
				   // preconditioner step, it is wise
				   // to invert the diagonal blocks of
				   // the matrix before starting the
				   // solver. Otherwise, it takes less
				   // memory, but the diagonal blocks
				   // are inverted in each
				   // preconditioner step,
				   // significantly slowing down the
				   // linear solving process.
  preconditioner.invert_diagblocks();

				   // After these preparations we are
				   // ready to start the linear solver.
  solver.solve (system_matrix, solution, right_hand_side,
		preconditioner);
};


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
				 // ``DerivativeApproximation'' class
				 // that computes the approximate
				 // gradients in a way similar to the
				 // ``GradientEstimation'' described
				 // in Step 9 of this tutorial. In
				 // fact, the
				 // ``DerivativeApproximation'' class
				 // was developed following the
				 // ``GradientEstimation'' class of
				 // Step 9. Relating to the
				 // discussion in Step 9, here we
				 // consider $h^{1+d/2}|\nabla_h
				 // u_h|$. Futhermore we note that we
				 // do not consider approximate second
				 // derivatives because solutions to
				 // the linear advection equation are
				 // in general not in H^2 but in H^1
				 // (to be more precise, in H^1_\beta)
				 // only.
template <int dim>
void DGMethod<dim>::refine_grid ()
{
				   // The ``DerivativeApproximation''
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
};


				 // The following ``run'' function is
				 // similar to previous examples. The
				 // only difference is that the
				 // problem is assembled and solved
				 // twice on each refinement step;
				 // first by ``assemble_system1'' that
				 // implements the first version and
				 // then by ``assemble_system2'' that
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
      system_matrix.reinit();
      right_hand_side.clear();
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
      if (difference>1e-13)
	std::cout << "solution1 and solution2 differ!!" << std::endl;
      else
	std::cout << "solution1 and solution2 coincide." << std::endl;
	
				       // Finally we perform the
				       // output.
      output_results (cycle);
    }
}

				 // The following ``main'' function is
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
};


