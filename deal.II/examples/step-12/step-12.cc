/* $Id$ */
/* Author: Ralf Hartmann, University of Heidelberg, 2000 */

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
#include <dofs/dof_handler.h>
#include <grid/grid_generator.h>
#include <grid/tria_accessor.h>
#include <grid/tria_iterator.h>
#include <dofs/dof_accessor.h>
#include <dofs/dof_tools.h>
#include <fe/fe_values.h>
#include <numerics/data_out.h>
#include <grid/grid_out.h>
#include <grid/grid_refinement.h>

				 // This is the first new file. It
				 // declares the MappingQ1 class that
				 // gives the standard bilinear
				 // mapping. For bilinear mappings use
				 // an object of this class rather
				 // than an object of the MappingQ(1)
				 // class, as the MappingQ1 class is
				 // optimized due to the
				 // pre-knowledge of the actual
				 // polynomial degree 1.
#include <fe/mapping_q1.h>

				 // Here the discontinuous finite
				 // elements are defined. They are
				 // used as all other finite elements.
#include <fe/fe_dgq.h>

				 // We are going to use the simplest
				 // possible solver, called richardson
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
#include <fstream>


				 // First we define the class
				 // representing the equation-specific
				 // functions. Both classes, ``RHS''
				 // and ``BoundaryValues'', are
				 // derived from the Function
				 // class. Only the ``value_list''
				 // function are implemented because
				 // only lists of function values are
				 // computed rather than single
				 // values.
template <int dim>
class RHS:  public Function<dim>
{
  public:
    RHS() {};
    
    virtual void value_list (const std::vector<Point<dim> > &points,
			     std::vector<double> &values,
			     const unsigned int component=0) const;
};


template <int dim>
class BoundaryValues:  public Function<dim>
{
  public:
    BoundaryValues() {};
    
    virtual void value_list (const std::vector<Point<dim> > &points,
			     std::vector<double> &values,
			     const unsigned int component=0) const;
};


				 // The class ``Beta'' that represents
				 // the vector valued flow field of
				 // the linear transport equation is
				 // not derived from the Function
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
				 // set to be zero.
template <int dim>
void RHS<dim>::value_list(const std::vector<Point<dim> > &,
			  std::vector<double> &values,
			  const unsigned int) const
{
  for (unsigned int i=0; i<values.size(); ++i)
    values[i]=0;
}

				 // The flow field is chosen to be
				 // circular, anticlockwise, and with
				 // the origin as midpoint.
template <>
void Beta<2>::value_list(const std::vector<Point<2> > &points,
			 std::vector<Point<2> > &values) const
{
  Assert(values.size()==points.size(),
	 ExcDimensionMismatch(values.size(),points.size()));

  for (unsigned int i=0; i<points.size(); ++i)
    {
      const Point<2> &p=points[i];
      Point<2> &beta=values[i];

      beta(0)=-p(1);
      beta(1)=p(0);
      beta/=sqrt(beta*beta);
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

				 // Next we define the equation-
				 // dependent and DG-method-dependent
				 // class ``DGTransportEquation''. Its
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
			    Vector<double> &cell_vector);
    
    void assemble_face_term1(const FEFaceValuesBase<dim>& fe_v,
			     const FEFaceValuesBase<dim>& fe_v_neighbor,
			     FullMatrix<double> &u_v_matrix,
			     FullMatrix<double> &un_v_matrix,
			     Vector<double> &cell_vector);

    void assemble_face_term2(const FEFaceValuesBase<dim>& fe_v,
			     const FEFaceValuesBase<dim>& fe_v_neighbor,
			     FullMatrix<double> &u_v_matrix,
			     FullMatrix<double> &un_v_matrix,
			     FullMatrix<double> &u_vn_matrix,
			     FullMatrix<double> &un_vn_matrix,
			     Vector<double> &cell_vector);
  private:
    Beta<dim> beta_function;
    RHS<dim> rhs_function;
    BoundaryValues<dim> boundary_function;
};

				 // ``u_v_matrix'' is a cell matrix,
				 // i.e. for a DG method of degree 1,
				 // it is of size 4 times 4, and
				 // ``cell_vector'' is of size 4.
				 // When this function is invoked,
				 // ``fe_v'' was reinited with the
				 // current cell before and includes
				 // all shape values needed.
template <int dim>
void DGTransportEquation<dim>::assemble_cell_term(
  const FEValues<dim>& fe_v,
  FullMatrix<double> &u_v_matrix,
  Vector<double> &cell_vector)
{
				   // First we ask ``fe_v'' for the
				   // shape grads, shape values and
				   // quadrature weights,
  const vector<vector<Tensor<1,2> > > &grad_v = fe_v.get_shape_grads ();
  const FullMatrix<double> &v = fe_v.get_shape_values ();
  const vector<double> &JxW = fe_v.get_JxW_values ();

				   // Then the flow field beta and the
				   // ``rhs_function'' are evaluated at
				   // the quadrature points,
  vector<Point<dim> > beta (fe_v.n_quadrature_points);
  vector<double> rhs (fe_v.n_quadrature_points);
  
  beta_function.value_list (fe_v.get_quadrature_points(), beta);
  rhs_function.value_list (fe_v.get_quadrature_points(), rhs);
  
				   // and the cell matrix and cell
				   // vector are assembled as in
				   // previous tutorial steps.  Here,
				   // the terms $-(u,\beta\cdot\nabla
				   // v)_K$ and $(f,v)_K$ are
				   // assembled.
  for (unsigned int point=0; point<fe_v.n_quadrature_points; ++point)
    for (unsigned int i=0; i<fe_v.dofs_per_cell; ++i) 
      {
	for (unsigned int j=0; j<fe_v.dofs_per_cell; ++j)
	  u_v_matrix(i,j) -= beta[point]*grad_v[i][point]*
			      v(j,point) *
			      JxW[point];
	
	cell_vector(i) += rhs[point] *v(i,point) * JxW[point];
      }
}


				 // The ``assemble_face_term1''
				 // function assembles the face terms
				 // corresponding to the first version
				 // of the DG method, cf. above. Then,
				 // the face terms are given as a sum
				 // of integrals over all cell
				 // boundaries.
				 //
				 // When this function is invoked,
				 // ``fe_v'' and ``fe_v_neighbor'' are
				 // already reinited with the current
				 // cell and the neighoring cell,
				 // respectively, as well as with the
				 // current face. Hence they provide
				 // the inner and outer shape values
				 // on the face.
				 //
				 // In addition to the cell matrix
				 // ``u_v_matrix'' and the
				 // ``cell_vector'' this function has
				 // got a new argument
				 // ``un_v_matrix'', that stores
				 // contributions to the system matrix
				 // that are based on outer values of
				 // u, see $\hat u_h$ in the
				 // Introduction, and inner values of
				 // v, see $v_h$. Here we note that
				 // ``un'' is the short notation for
				 // ``u_neighbor'' and represents
				 // $\hat u_h$.
template <int dim>
void DGTransportEquation<dim>::assemble_face_term1(
  const FEFaceValuesBase<dim>& fe_v,
  const FEFaceValuesBase<dim>& fe_v_neighbor,      
  FullMatrix<double> &u_v_matrix,
  FullMatrix<double> &un_v_matrix,
  Vector<double> &cell_vector)
{
				   // Again, we ask the FEValues
				   // objects for the shape values and
				   // the quadrature weights
  const FullMatrix<double> &v = fe_v.get_shape_values ();
  const FullMatrix<double> &v_neighbor = fe_v_neighbor.get_shape_values ();  
  const vector<double> &JxW = fe_v.get_JxW_values ();
				   // but also for the normals.
  const vector<Point<dim> > &normals = fe_v.get_normal_vectors ();

				   // We also evaluate the flow field
				   // at the quadrature points
  vector<Point<dim> > beta (fe_v.n_quadrature_points);
  
  beta_function.value_list (fe_v.get_quadrature_points(), beta);

				   // and the boundary values if the
				   // current face belongs to the
				   // boundary.
  vector<double> g(fe_v.n_quadrature_points);
  DoFHandler<dim>::face_iterator face=fe_v.get_face();
  if (face->at_boundary())
    boundary_function.value_list (fe_v.get_quadrature_points(), g);

				   // Then we assemble the cell matrix
				   // and cell vector according to the
				   // DG method given in the
				   // introduction.
  for (unsigned int point=0; point<fe_v.n_quadrature_points; ++point)
    {
      double beta_n=beta[point] * normals[point];
      if (beta_n>0)
					 // The term $(\beta\cdot n
					 // u,v)_{\partial K_+}$,
	for (unsigned int i=0; i<fe_v.dofs_per_cell; ++i)
	  for (unsigned int j=0; j<fe_v.dofs_per_cell; ++j)
	    u_v_matrix(i,j) += beta_n *
			       v(j,point) *
			       v(i,point) *
			       JxW[point];
      else
	{
					   // at the boundary the term
					   // $(\beta\cdot n
					   // g,v)_{\partial
					   // K_-\cap\partial\Omega}$,
	  if (face->at_boundary())
	    for (unsigned int i=0; i<fe_v.dofs_per_cell; ++i)
	      cell_vector(i) -= beta_n *
				g[point] *
				v(i,point) *
				JxW[point];
	  else
					     // and on inner faces the
					     // term $(\beta\cdot n
					     // \hat u,v)_{\partial
					     // K_-}$
	    for (unsigned int i=0; i<fe_v.dofs_per_cell; ++i)
	      for (unsigned int k=0; k<fe_v_neighbor.dofs_per_cell; ++k)
		un_v_matrix(i,k) += beta_n *
				    v_neighbor(k,point) *
				    v(i,point) *
				    JxW[point];
	}
    }
}

				 // Now we look at the assembling
				 // function that assembles the face
				 // terms corresponding to the second
				 // version of the DG method,
				 // cf. above. Then, the face terms
				 // are given as a sum of integrals
				 // over all faces.  Here we need two
				 // additional cell matrices
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
  FullMatrix<double> &un_vn_matrix,
  Vector<double> &cell_vector)
{
				   // the first few lines are the same
  const FullMatrix<double> &v = fe_v.get_shape_values ();
  const FullMatrix<double> &v_neighbor = fe_v_neighbor.get_shape_values ();  
  const vector<double> &JxW = fe_v.get_JxW_values ();
  const vector<Point<dim> > &normals = fe_v.get_normal_vectors ();

  vector<Point<dim> > beta (fe_v.n_quadrature_points);
  
  beta_function.value_list (fe_v.get_quadrature_points(), beta);

  vector<double> g(fe_v.n_quadrature_points);
  DoFHandler<dim>::face_iterator face=fe_v.get_face();
  if (face->at_boundary())
    boundary_function.value_list (fe_v.get_quadrature_points(), g);

  for (unsigned int point=0; point<fe_v.n_quadrature_points; ++point)
    {
      double beta_n=beta[point] * normals[point];
      if (beta_n>0)
	{
					   // This terms we've already seen,
	  for (unsigned int i=0; i<fe_v.dofs_per_cell; ++i)
	    for (unsigned int j=0; j<fe_v.dofs_per_cell; ++j)
	      u_v_matrix(i,j) += beta_n *
				 v(j,point) *
				 v(i,point) *
				 JxW[point];

					   // on inner faces we
					   // additionally have the
					   // term $(\beta\cdot n
					   // u,\hat v)_{\partial K_+},
	  if (!face->at_boundary())
	    for (unsigned int k=0; k<fe_v_neighbor.dofs_per_cell; ++k)
	      for (unsigned int j=0; j<fe_v.dofs_per_cell; ++j)
		u_vn_matrix(k,j) -= beta_n *
				    v(j,point) *
				    v_neighbor(k,point) *
				    JxW[point];
	}
      else
	{
					   // this one we already know,
	  if (face->at_boundary())
	    for (unsigned int i=0; i<fe_v.dofs_per_cell; ++i)
	      cell_vector(i) -= beta_n *
				g[point] *
				v(i,point) *
				JxW[point];
	  else
	    {
					       // this one also,
	      for (unsigned int i=0; i<fe_v.dofs_per_cell; ++i)
		for (unsigned int l=0; l<fe_v_neighbor.dofs_per_cell; ++l)
		  un_v_matrix(i,l) += beta_n *
				      v_neighbor(l,point) *
				      v(i,point) *
				      JxW[point];

					       // and this is another
					       // new one:
					       // $(\beta\cdot n \hat
					       // u,\hat v)_{\partial
					       // K_-}$.
	      for (unsigned int k=0; k<fe_v_neighbor.dofs_per_cell; ++k)
		for (unsigned int l=0; l<fe_v_neighbor.dofs_per_cell; ++l)
		  un_vn_matrix(k,l) -= beta_n *
				       v_neighbor(l,point) *
				       v_neighbor(k,point) *
				       JxW[point];
	    }
	}
    }
}


				 // After these preparations, we
				 // proceed with the main part of this
				 // program. The main class, here
				 // called ``DGMethod'' is basically
				 // the main class of step 6. One of
				 // the differences is that there's no
				 // ConstraintMatrix object. This is,
				 // because there are no hanging nodes
				 // in DG discretizations.
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
    MappingQ1<dim>       mapping;
    
				     // Furthermore we want to
				     // use DG elements of degree 1
				     // (but this is only specified in
				     // the constructor):
    FE_DGQ<dim>          fe;
    DoFHandler<dim>      dof_handler;

    SparsityPattern      sparsity_pattern;
    SparseMatrix<double> system_matrix;

				     // We define the quadrature
				     // formulae for the cell and the
				     // face terms of the
				     // discretization.
    QGauss4<dim>   quadrature;
    QGauss4<dim-1> face_quadrature;
    
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
    DGTransportEquation<dim>     dg;
};



				 // Now for the implementation of the
				 // main class. Constructor and
				 // destructor follow the same
				 // pattern that was used previously,
				 // so we need not comment on these
				 // two functions:  
template <int dim>
DGMethod<dim>::DGMethod () :
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
				   // number of matrix entries is
				   // needed when all neighbors of a
				   // cell are once more refined than
				   // the cell under consideration.
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


				 // We proceed with the
				 // ``assemble_system1'' function that
				 // implements the DG discretization
				 // in its first version. This
				 // function repeatedly calls the
				 // ``assemble_cell_term'' and
				 // ``assemble_face_term1'' functions
				 // of the DGTransportEquation object.
				 // The ``assemble_face_term1''
				 // function takes two
				 // FEFaceValuesBase objects; one for
				 // the shape functions on the current
				 // cell and the other for shape
				 // functions on the neighboring cell
				 // under consideration. Both objects
				 // are either of class FEFaceValues
				 // or of class FESubfaceValues (both
				 // derived from FEFaceValuesBase)
				 // according to following cases
				 // already mentioned in the
				 // introduction:
				 //
				 // 1. face is at boundary (current
				 // cell: FEFaceValues, neighboring
				 // cell does not exist);
				 //
				 // 2. neighboring cell is finer
				 // (current cell: FESubfaceValues,
				 // neighboring cell: FEFaceValues);
				 //
				 // 3. neighboring cell is of the same
				 // refinement level (both, current
				 // and neighboring cell:
				 // FEFaceValues);
				 //
				 // 4. neighboring cell is coarser
				 // (current cell: FEFaceValues,
				 // neighboring cell:
				 // FESubfaceValues).
				 //
				 // If we considered globally refined
				 // meshes then only cases 1 and 3
				 // would occur. But as we consider
				 // also locally refined meshes we
				 // need to distinguish all four cases
				 // making the following assembling
				 // function a bit longish.
template <int dim>
void DGMethod<dim>::assemble_system1 () 
{
  const unsigned int dofs_per_cell = dof_handler.get_fe().dofs_per_cell;
  vector<unsigned int> dofs (dofs_per_cell);
  vector<unsigned int> dofs_neighbor (dofs_per_cell);

				   // First we create the Update flags
				   // for the FEValues and the
				   // FEFaceValues objects.
  UpdateFlags update_flags = UpdateFlags(update_values
					 | update_gradients
					 | update_q_points
					 | update_JxW_values);
  
				   // Note, that on faces we do not
				   // need gradients but we need
				   // normal vectors.
  UpdateFlags face_update_flags = UpdateFlags(update_values
					      | update_q_points
					      | update_JxW_values
					      | update_normal_vectors);
  
				   // On the neighboring cell we only
				   // need the shape values. Given a
				   // specific face, the quadrature
				   // points and `JxW values' are the
				   // same as for the current cells,
				   // the normal vectors are known to
				   // be the negative of the normal
				   // vectors of the current cell.
  UpdateFlags neighbor_face_update_flags = UpdateFlags(update_values);
   
				   // Then we create the FEValues
				   // object. Note, that since version
				   // 3.2.0 the constructor of this
				   // class takes a Mapping object as
				   // first argument. Although the
				   // constructor without Mapping
				   // argument is still supported it
				   // is recommended to use the new
				   // constructor. This reduces the
				   // effect of `hidden magic' (the
				   // old constructor implicitely
				   // assumes a MappingQ1 mapping) and
				   // makes it easier to change the
				   // Mapping object later.
  FEValues<dim> fe_v (
    mapping, fe, quadrature, update_flags);
  
				   // Similarly we create the
				   // FEFaceValues and FESubfaceValues
				   // objects for both, the current
				   // and the neighboring cell. Within
				   // the following nested loop over
				   // all cells and all faces of the
				   // cell they will be reinited to
				   // the current cell and the face
				   // (and subface) number.
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
				   // our the convention that `un' is
				   // the short cut for `u_neighbor'
				   // and represents the $u_hat$, see
				   // introduction.
  FullMatrix<double> u_v_matrix (dofs_per_cell, dofs_per_cell);
  FullMatrix<double> un_v_matrix (dofs_per_cell, dofs_per_cell);

  Vector<double>  cell_vector (dofs_per_cell);

				   // Furthermore we need some cell
				   // and face iterators
  DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(),
					endc = dof_handler.end();
  DoFHandler<dim>::face_iterator face;
  DoFHandler<dim>::cell_iterator neighbor;
  DoFHandler<dim>::active_cell_iterator neighbor_child;

				   // Now we start the loop over all
				   // active cells
  for (;cell!=endc; ++cell) 
    {
				       // and reinit the FEValues
				       // object for the current cell,
      fe_v.reinit (cell);

				       // Call the function that
				       // assembles the cell
				       // terms. The first argument is
				       // the FEValues that was
				       // already reinited on current
				       // the cell.
      dg.assemble_cell_term(fe_v,
			    u_v_matrix,
			    cell_vector);

				       // As in previous example steps
				       // the vector `dofs' includes
				       // the dof_indices.
      cell->get_dof_indices (dofs);

				       // This is the start of the
				       // nested loop over all faces.
      for (unsigned int face_no=0; face_no<GeometryInfo<dim>::faces_per_cell; ++face_no)
	{
					   // First we set the face
					   // iterator.
	  face = cell->face(face_no);

					   // Now we distinguish the
					   // four different cases in
					   // the ordering mentioned
					   // above. We start with
					   // faces belonging to the
					   // boundary of the domain.
	  if (face->at_boundary())
	    {
					       // We reinit the
					       // FEFaceValues object
					       // to the current face
	      fe_v_face.reinit (cell, face_no);

					       // and assemble the
					       // corresponding face
					       // terms. Here, the
					       // second and fourth
					       // arguments are only
					       // dummy arguments. On
					       // the boundary of the
					       // domain the
					       // ``assemble_face_term1''
					       // function will not
					       // access to shape
					       // values on the
					       // non-existent
					       // neighboring
					       // cell. Also,
					       // ``un_v_matrix'' will
					       // be unchanged.
	      dg.assemble_face_term1(fe_v_face,
				     fe_v_face,
				     u_v_matrix,
				     un_v_matrix,
				     cell_vector);
	    }
	  else
	    {
					       // When we are not on the
					       // boundary of the
					       // domain then there
					       // must exist a
					       // neighboring cell.
	      neighbor = cell->neighbor(face_no);

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
					       // known to be only
					       // ONCE more refined
					       // than the current
					       // cell. Furthermore
					       // also the face is
					       // once more refined,
					       // i.e. it has
					       // children.
	      if (face->has_children())
		{
						   // first we store
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
		      neighbor_child = neighbor->child(GeometryInfo<dim>::
						       child_cell_on_face(neighbor2,subface_no));

						       // As these are
						       // quite
						       // complicated
						       // indirections
						       // we check for
						       // the internal
						       // consistency.
		      Assert (neighbor_child->face(neighbor2) == face->child(subface_no),
			      ExcInternalError());
		      Assert (!neighbor_child->has_children(), ExcInternalError());

						       // As already
						       // mentioned
						       // above for
						       // this case
						       // (case 2) we
						       // employ the
						       // FESubfaceValues
						       // of the
						       // current
						       // cell, here
						       // reinited for
						       // the current
						       // cell, face
						       // and subface,
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
					     un_v_matrix,
					     cell_vector);
		      
						       // get dof
						       // indices of
						       // the
						       // neighbor_child
						       // cell
		      neighbor_child->get_dof_indices (dofs_neighbor);
		      						
						       // distribute
						       // cell matrix
						       // to the
						       // system_matrix
		      for (unsigned int i=0; i<dofs_per_cell; ++i)
			for (unsigned int k=0; k<dofs_per_cell; ++k)
			  system_matrix.add(dofs[i], dofs_neighbor[k],
					    un_v_matrix(i,k));

		      				       // In the
						       // ``assemble_face_term1''
						       // function contributions to
						       // the cell matrices and the
						       // cell vector are only
						       // ADDED. Therefore on each
						       // subface we need to reset the
						       // un_v_matrix
						       // to zero, before assembling
						       // the face terms corresponding
						       // to the following neighbor_child cell.
		      un_v_matrix.clear();
		    }
		}
					       // End of ``if
					       // (face->has_children())''
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
						   // Like before we
						   // store which
						   // number the
						   // current cell has
						   // in the list of
						   // neighbors of the
						   // neighboring
						   // cell.
		      const unsigned int neighbor2=cell->neighbor_of_neighbor(face_no);

						       // We reinit
						       // the
						       // FEFaceValues
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
					     un_v_matrix,
					     cell_vector);
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
						       // FEFaceValues
						       // and assemble
						       // the face
						       // terms.
		      fe_v_face.reinit (cell, face_no);
		      fe_v_subface_neighbor.reinit (neighbor, neighbor_face_no,
						    neighbor_subface_no);
		      
		      dg.assemble_face_term1(fe_v_face,
					     fe_v_subface_neighbor,
					     u_v_matrix,
					     un_v_matrix,
					     cell_vector);
		    }

						   // Get dof indices
						   // of the
						   // neighbor_child
						   // cell,
		  neighbor->get_dof_indices (dofs_neighbor);
		 						
						   // distribute the
						   // un_v_matrix,
		  for (unsigned int i=0; i<dofs_per_cell; ++i)
		    for (unsigned int k=0; k<dofs_per_cell; ++k)
		      system_matrix.add(dofs[i], dofs_neighbor[k],
					un_v_matrix(i,k));
		  
						   // and clear the
						   // ``un_v_matrix''
						   // on each face.
		  un_v_matrix.clear();
		}
					       // End of ``face not at boundary'':
	    }
					   // End of loop over all faces:
	}
      
				       // Finally we distribute the
				       // u_v_matrix,
      for (unsigned int i=0; i<dofs_per_cell; ++i)
	for (unsigned int j=0; j<dofs_per_cell; ++j)
	  system_matrix.add(dofs[i], dofs[j], u_v_matrix(i,j));
      
				       // the cell vector
      for (unsigned int i=0; i<dofs_per_cell; ++i)
	right_hand_side(dofs[i]) += cell_vector(i);

				       // and clear them both.
      u_v_matrix.clear ();
      cell_vector.clear ();
    }
};



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
  vector<unsigned int> dofs (dofs_per_cell);
  vector<unsigned int> dofs_neighbor (dofs_per_cell);

  UpdateFlags update_flags = UpdateFlags(update_values
					 | update_gradients
					 | update_q_points
					 | update_JxW_values);
  
  UpdateFlags face_update_flags = UpdateFlags(update_values
					      | update_q_points
					      | update_JxW_values
					      | update_normal_vectors);
   
  UpdateFlags neighbor_face_update_flags = UpdateFlags(update_values);

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
  
				   // Additionally we need following
				   // two cell matrices, both for face
				   // term that include test function
				   // ``vn'' (shape functions of the
				   // neighboring cell). To be more
				   // precise, the first matrix will
				   // include the `u and vn terms' and
				   // the second that will include the
				   // `un and vn terms'.
  FullMatrix<double> u_vn_matrix (dofs_per_cell, dofs_per_cell);
  FullMatrix<double> un_vn_matrix (dofs_per_cell, dofs_per_cell);
  
  Vector<double>  cell_vector (dofs_per_cell);

				   // Furthermore, here we define a
				   // dummy matrix and rhs to
				   // emphasize when arguments of the
				   // ``assemble_face_term2''
				   // functions will not be access.
  FullMatrix<double> dummy_matrix;
  Vector<double>     dummy_rhs;

				   // The following lines are roughly
				   // the same as in the previous
				   // function.
  DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(),
					endc = dof_handler.end();
  DoFHandler<dim>::face_iterator face;
  DoFHandler<dim>::cell_iterator neighbor;
  DoFHandler<dim>::cell_iterator neighbor_child;

  for (;cell!=endc; ++cell) 
    {
      fe_v.reinit (cell);

      dg.assemble_cell_term(fe_v,
			    u_v_matrix,
			    cell_vector);
      
      cell->get_dof_indices (dofs);

      for (unsigned int face_no=0; face_no<GeometryInfo<dim>::faces_per_cell; ++face_no)
	{
	  face = cell->face(face_no);

					   // Case 1:
	  if (face->at_boundary())
	    {
	      fe_v_face.reinit (cell, face_no);

	      dg.assemble_face_term2(fe_v_face,
				     fe_v_face,
				     u_v_matrix,
				     dummy_matrix,
				     dummy_matrix,
				     dummy_matrix,
				     cell_vector);
	    }
	  else
	    {
	      Assert (cell->neighbor(face_no).state() == valid, ExcInternalError());
	      neighbor = cell->neighbor(face_no);

					       // Case 2:
	      if (face->has_children())
		{
		  const unsigned int neighbor2=
		    cell->neighbor_of_neighbor(face_no);
		  
		  for (unsigned int subface_no=0;
		       subface_no<GeometryInfo<dim>::subfaces_per_face;
		       ++subface_no)
		    {
		      neighbor_child = neighbor->child(
			GeometryInfo<dim>::child_cell_on_face(neighbor2,subface_no));
		      Assert (neighbor_child->face(neighbor2) == face->child(subface_no),
			      ExcInternalError());
		      Assert (!neighbor_child->has_children(), ExcInternalError());

		      fe_v_subface.reinit (cell, face_no, subface_no);
		      fe_v_face_neighbor.reinit (neighbor_child, neighbor2);

		      dg.assemble_face_term2(fe_v_subface,
					     fe_v_face_neighbor,
					     u_v_matrix,
					     un_v_matrix,
					     u_vn_matrix,
					     un_vn_matrix,
					     dummy_rhs);
		  
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
		      
		      un_v_matrix.clear();
		      u_vn_matrix.clear();
		      un_vn_matrix.clear();
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

		      fe_v_face.reinit (cell, face_no);
		      fe_v_face_neighbor.reinit (neighbor, neighbor2);
		      
		      dg.assemble_face_term2(fe_v_face,
					     fe_v_face_neighbor,
					     u_v_matrix,
					     un_v_matrix,
					     u_vn_matrix,
					     un_vn_matrix,
					     dummy_rhs);

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
		      
		      un_v_matrix.clear();
		      u_vn_matrix.clear();
		      un_vn_matrix.clear();
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
      
      u_v_matrix.clear ();
      cell_vector.clear ();
    }
};

				 // For this simple solver we use the
				 // simplest possible solver, called
				 // richardson iteration, that
				 // represents a simple defect
				 // correction. This, in combination
				 // with a block SSOR preconditioner,
				 // that uses the special block matrix
				 // structur of system matrices
				 // arising from DG
				 // discretizations. The size of these
				 // blocks are the number of DoFs
				 // per cell. Here, we use a SSOR
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
				   // solver. Otherwise, the diagonal
				   // blocks are inverted in each
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
				 // namely the gradients of the
				 // solution. As here we consider the
				 // DG(1) method (i.e. we use
				 // piecewise bilinear shape
				 // functions) we could simply compute
				 // the gradients on each cell. But we
				 // do not want to base our refinement
				 // indicator on the gradients on each
				 // cell only, but want to base them
				 // also on jumps of the discontinuous
				 // solution function over faces
				 // between neighboring cells. The
				 // simpliest way of doing that is to
				 // compute approximative gradients by
				 // difference quotients including the
				 // cell under consideration and its
				 // neighbors. This is done by the
				 // DerivativeApproximation class that
				 // computes the approximate
				 // gradients in a way similar to the
				 // GradientEstimation described in
				 // Step 9 of this tutorial. According
				 // to the argumentation in Step 9,
				 // here we consider
				 // $h^{1+d/2}|\nabla_h
				 // u_h|$. Futhermore we note that we
				 // do not consider approximate
				 // second derivatives because
				 // solutions to the linear advection
				 // equation are in general not in H^2
				 // but in H^1 (to be more precise, in
				 // H^1_\beta) only.
template <int dim>
void DGMethod<dim>::refine_grid ()
{
				   // The DerivativeApproximation
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
  DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(),
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
  cout << "Writing grid to <" << filename << ">..." << endl;
  std::ofstream eps_output (filename.c_str());

  GridOut grid_out;
  grid_out.write_eps (triangulation, eps_output);
  
				   // Output of the solution in
				   // gnuplot format.
  filename = "sol-";
  filename += ('0' + cycle);
  Assert (cycle < 10, ExcInternalError());
  
  filename += ".gnuplot";
  cout << "Writing solution to <" << filename << ">..." << endl;
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
      cout << "Time of assemble_system1: " << assemble_timer() << endl;
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
      cout << "Time of assemble_system2: " << assemble_timer() << endl;
      solve (solution2);

				       // To make sure that both
				       // versions of the DG method
				       // yield the same
				       // discretization and hence the
				       // same solution we check the
				       // two solutions for equality.
      solution1-=solution2;
      const double difference=solution1.linfty_norm();
      if (difference<1e-13)
	cout << "solution1 and solution2 do not differ." << endl;

				       // Finally we perform the
				       // output.
      output_results (cycle);
    }
}



int main () 
{
  DGMethod<2> dgmethod_2d;
  dgmethod_2d.run ();

  return 0;
};


