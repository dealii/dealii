/*----------------------------   fe_values.h     ---------------------------*/
/*      $Id$                 */
#ifndef __fe_values_H
#define __fe_values_H
/*----------------------------   fe_values.h     ---------------------------*/


#include <lac/dfmatrix.h>
#include <grid/point.h>
#include <base/exceptions.h>
#include <grid/tria.h>
#include <fe/fe_update_flags.h>


// forward declarations
template <int dim> class Boundary;
template <int dim> class FiniteElement;
template <int dim> class Quadrature;



/**
  This class offers a multitude of arrays and other fields which are used by
   the derived classes #FEValues# and #FEFaceValues#. In principle, it is the
   back end of the front end for the unification of a certain finite element
   and a quadrature formula which evaluates certain aspects of the finite
   element at quadrature points.
   
   This class is an optimization which avoids evaluating the shape functions
   at the quadrature points each time a quadrature takes place. Rather, the
   values and gradients (and possibly higher order derivatives in future
   versions of this library) are evaluated once and for all on the unit
   cell or face before doing the quadrature itself. Only the Jacobian matrix of
   the transformation from the unit cell or face to the real cell or face and
   the integration points in real space are calculated each time we move on
   to a new face.

   Actually, this class does none of the evaluations at startup itself; this is
   all done by the derived classes. It only offers the basic functionality,
   like providing those fields that are common to the derived classes and
   access to these fields. Any computations are in the derived classes. See there
   for more information.

   It has support for the restriction of finite elements to faces of cells or
   even to subfaces (i.e. refined faces). For this purpose, it offers an array
   of matrices of ansatz function values, rather than one. Since the value of
   a function at a quadrature point is an invariant under the transformation
   from the unit cell to the real cell, it is only evaluated once upon startup.
   However, when considering the restriction of a finite element to a face of
   a cell (using a given quadrature rule), we may be tempted to compute the
   restriction to all faces at startup (thus ending in four array of ansatz
   function values in two dimensions, one per face, and even more in higher
   dimensions) and let the respective #reinit# function of the derived classes
   set a number which of the fields is to be taken when the user requests the
   function values. This is done through the #selected_dataset# variable. See
   the derived classes and the #get_values# function for the exact usage of
   this variable.

   
   {\bf Member functions}

   The functions of this class fall into different cathegories:
   \begin{itemize}
   \item #shape_value#, #shape_grad#, etc: return one of the values
     of this object at a time. In many cases you will want to get
     a whole bunch at a time for performance or convenience reasons,
     then use the #get_*# functions.
    
   \item #get_shape_values#, #get_shape_grads#, etc: these return
     a reference to a whole field. Usually these fields contain
     the values of all ansatz functions at all quadrature points.

   \item #get_function_values#, #get_function_gradients#: these
     two functions offer a simple way to avoid the detour of the
     ansatz functions, if you have a finite solution (resp. the
     vector of values associated with the different ansatz functions.)
     Then you may want to get information from the restriction of
     the finite element function to a certain cell, e.g. the values
     of the function at the quadrature points or the values of its
     gradient. These two functions provide the information needed:
     you pass it a vector holding the finite element solution and the
     functions return the values or gradients of the finite element
     function restricted to the cell which was given last time the
     #reinit# function was given.
    
     Though possible in principle, these functions do not call the
     #reinit# function, you have to do so yourself beforehand. On the
     other hand, a copy of the cell iterator is stored which was used
     last time the #reinit# function was called. This frees us from
     the need to pass the cell iterator again to these two functions,
     which guarantees that the cell used here is in sync with that used
     for the #reinit# function. You should, however, make sure that
     nothing substantial happens to the #DoFHandler# object or any
     other involved instance between the #reinit# and the #get_function_*#
     functions are called.

   \item #reinit#: initialize the #FEValues# object for a certain cell.
     This function is not in the present class but only in the derived
     classes and has a variable call syntax. 
     See the docs for the derived classes for more information.
  \end{itemize}

   @author Wolfgang Bangerth, 1998
*/
template <int dim>
class FEValuesBase {
  public:
        
				     /**
				      * Number of quadrature points.
				      */
    const unsigned int n_quadrature_points;

				     /**
				      * Total number of shape functions
				      * per cell. If we use this base class
				      * to evaluate a finite element on
				      * faces of cells, this is still the
				      * number of degrees of freedom per
				      * cell, not per face.
				      */
    const unsigned int total_dofs;

				     /**
				      * Constructor. Set up the array sizes
				      * with #n_q_points# quadrature points,
				      * #n_ansatz_points# ansatz points (on
				      * the cell or face), #n_dof# ansatz
				      * functions per cell and with the
				      * given pattern to update the fields
				      * when the #reinit# function of the
				      * derived classes is called. The
				      * fields themselves are not set up,
				      * this must happen in the derived
				      * class's constructor, only the sizes
				      * are set correctly.
				      */
    FEValuesBase (const unsigned int n_q_points,
		  const unsigned int n_ansatz_points,
		  const unsigned int n_dofs,
		  const unsigned int n_values_array,
		  const UpdateFlags update_flags);
    

				     /**
				      * Return the value of the #i#th shape
				      * function at the #j# quadrature point
				      * on the cell, face or subface selected
				      * the last time the #reinit# function
				      * of the derived class was called.
				      */
    double shape_value (const unsigned int i,
			const unsigned int j) const;

				     /**
				      * Return a pointer to the matrix holding
				      * all values of shape functions at all
				      * integration points, on the present cell,
				      * face or subface selected
				      * the last time the #reinit# function
				      * of the derived class was called.
				      * For the format of this matrix, see the
				      * documentation for the matrix itself.
				      */
    const dFMatrix & get_shape_values () const;

				     /**
				      * Return the values of the finite
				      * element function characterized
				      * by #fe_function# restricted to
				      * the cell, face or subface selected
				      * the last time the #reinit# function
				      * of the derived class was called,
				      * at the quadrature points.
				      *
				      * The function assumes that the
				      * #values# object already has the
				      * right size.
				      */
    void get_function_values (const dVector      &fe_function,
			      vector<double>     &values) const;

    				     /**
				      * Return the gradient of the #i#th shape
				      * function at the #j# quadrature point.
				      * If you want to get the derivative in
				      * one of the coordinate directions, use
				      * the appropriate function of the #Point#
				      * class to extract one component. Since
				      * only a reference to the gradient's value
				      * is returned, there should be no major
				      * performance drawback.
				      * The function returns the gradient on the
				      * real element, not the reference element.
				      */
    const Point<dim> & shape_grad (const unsigned int i,
				   const unsigned int j) const;

				     /** 
				      * Return a pointer to the matrix holding
				      * all gradients of shape functions at all
				      * integration points, on the present cell.
				      * For the format of this matrix, see the
				      * documentation for the matrix itself.
				      */
    const vector<vector<Point<dim> > > & get_shape_grads () const;
    
				     /**
				      * Return the gradients of the finite
				      * element function characterized
				      * by #fe_function# restricted to
				      * #cell# at the quadrature points.
				      *
				      * The function assumes that the
				      * #gradients# object already has the
				      * right size.
				      */
    void get_function_grads (const dVector       &fe_function,
			     vector<Point<dim> > &gradients) const;

				     /**
				      * Return the position of the #i#th
				      * quadrature point in real space.
				      *
				      * If this object is used to evaluate
				      * finite elements on faces of cells,
				      * and for curved boundary cells, using
				      * biquadratic or higher mappings
				      * of the unit cell to the real cell,
				      * these points may not be on the
				      * plane submannifold on which the
				      * vertices of the face lie.
				      */
    const Point<dim> & quadrature_point (const unsigned int i) const;

				     /**
				      * Return a pointer to the vector of
				      * quadrature points.
				      */
    const vector<Point<dim> > & get_quadrature_points () const;

				     /**
				      * Return the point in real space where
				      * the #i#th ansatz function is located
				      * (location is in the sense of where it
				      * assumes its nominal properties, e.g. at
				      * the vertex of a cell, at the center of
				      * a line, etc).
				      *
				      * This function is needed for the
				      * interpolation problem: if we want to
				      * transfer a continuous function to a
				      * finite element function by interpolation
				      * we have to take the continuous
				      * function's value at the ansatz function
				      * locations.
				      *
				      * For the evaluation of finite elements on
				      * faces of cells, #i# is the number
				      * of the ansatz function on the face, not
				      * on the cell.
				      */
    const Point<dim> & ansatz_point (const unsigned int i) const;

				     /**
				      * Return a pointer to the vector of points
				      * denoting the location of the ansatz
				      * functions.
				      */
    const vector<Point<dim> > & get_ansatz_points () const;
    
				     /**
				      * Return the Jacobi determinant times
				      * the weight of the #i#th quadrature
				      * point.
				      *
				      * If faces of cells are concerned,
				      * the jacobi determinant is that of the
				      * transformation of the unit face to
				      * the present face, not of the unit
				      * cell to the real cell (unlike for
				      * the #jacobi_matrix# array of the
				      * derived classes which store the cell
				      * transformation's Jacobi matrix in
				      * all cases).
				      */
    double JxW (const unsigned int i) const;

				     /**
				      * Return a pointer to the array holding
				      * the JxW values at the different
				      * quadrature points.
				      */
    const vector<double> & get_JxW_values () const;

				     /**
				      * Exception
				      */
    DeclException2 (ExcInvalidIndex,
		    int, int,
		    << "The index " << arg1
		    << " is out of range, it should be less than " << arg2);
				     /**
				      * Exception
				      */
    DeclException0 (ExcAccessToUninitializedField);
				     /**
				      * Exception
				      */
    DeclException0 (ExcCannotInitializeField);
				     /**
				      * Exception
				      */
    DeclException2 (ExcWrongVectorSize,
		    int, int,
		    << "Vector has wrong size " << arg1
		    << ", expected size " << arg2);
				     /**
				      * Exception
				      */
    DeclException0 (ExcInternalError);
				     /**
				      * Exception
				      */
    DeclException0 (ExcNotImplemented);
    
  protected:
				     /**
				      * Store the values of the shape functions
				      * at the quadrature points. Rows in the
				      * matrices denote the values of a single
				      * shape function at the different points,
				      * columns are for a single point with the
				      * different shape functions.
				      *
				      * For cell values, the vector contains
				      * only one entry, representing the
				      * restriction of the finite element ansatz
				      * space to a cell. For face values, the
				      * vector contains as many elements as
				      * there are faces, for subfaces the same
				      * applies. Which of the matrices is active
				      * is determined by the #selected_dataset#
				      * variable.
				      */
    vector<dFMatrix>     shape_values;

    				     /**
				      * Store the gradients of the shape
				      * functions at the quadrature points.
				      * Since unfortunately the full matrix
				      * classes of DEAL are not templated,
				      * we have to store them in an
				      * archetypic style.
				      *
				      * This field is reset each time
				      * #reinit# is called and contains the
				      * gradients on the real element, rather
				      * than on the reference element.
				      */
    vector<vector<Point<dim> > >  shape_gradients;

				     /**
				      * Store an array of the weights of the
				      * quadrature points. This array is
				      * set up upon construction.
				      *
				      * If faces rather than cells are
				      * considered, the weights are stored
				      * only once still, since they are
				      * not transformed and are thus the same
				      * for all faces.
				      */
    vector<double>       weights;

				     /**
				      * Store an array of weights times the
				      * Jacobi determinant at the quadrature
				      * points. This function is reset each time
				      * #reinit# is called. The Jacobi determinant
				      * is actually the reciprocal value of the
				      * Jacobi matrices stored in this class,
				      * see the general documentation of this
				      * class for more information.
				      */
    vector<double>       JxW_values;

				     /**
				      * Array of quadrature points. This array
				      * is set up upon calling #reinit# and
				      * contains the quadrature points on the
				      * real element, rather than on the
				      * reference element.
				      */
    vector<Point<dim> >  quadrature_points;

    				     /**
				      * Array of points denoting the off-point
				      * of the ansatz functions. In real space
				      * (no-one seems to need the off-point
				      * on the unit cell, so no function is
				      * provided for this).
				      */
    vector<Point<dim> >  ansatz_points;
    
				     /**
				      * Store the jacobi matrices at the
				      * different quadrature points. This field
				      * is set each time #reinit# is called.
				      *
				      * If faces rather than cells are considered
				      * this is the Jacobi matrix of the
				      * transformation of the unit cell to the
				      * real cell, not of the unit face to the
				      * face. We need this full matrix for the
				      * transformation of the gradients to the
				      * real cell.
				      */
    vector<dFMatrix>     jacobi_matrices;

				     /**
				      * Store which of the data sets in the
				      * #shape_values# array is presently
				      * active. This variable is set by the
				      * #reinit# functions of the derived
				      * classes. For the exact meaning see
				      * there and in the doc for this class.
				      */
    unsigned int         selected_dataset;
    
				     /**
				      * Store which fields are to be updated by
				      * the reinit function.
				      */
    UpdateFlags         update_flags;

				     /**
				      * Store the cell selected last time
				      * the #reinit# function was called
				      * to make access
				      * to the #get_function_*# functions
				      * safer.
				      */
    DoFHandler<dim>::cell_iterator present_cell;
};



/**
  Represent a finite element evaluated with a specific quadrature rule on
  a cell.

  The unit cell is defined to be the tensor product of the interval $[0,1]$
  in the present number of dimensions. In part of the literature, the convention
  is used that the unit cell be the tensor product of the interval $[-1,1]$,
  which is to distinguished properly.

  Objects of this class store a multitude of different values needed to
  do the assemblage steps on real cells rather than on the unit cell. Among
  these values are the values and gradients of the shape functions at the
  quadrature points on the real and the unit cell, the location of the
  quadrature points on the real and on the unit cell, the weights of the
  quadrature points, the Jacobian matrices of the mapping from the unit to
  the real cell at the quadrature points and so on.

  The Jacobian matrix is defined to be
  $$ J_{ij} = {d\xi_i \over dx_j} $$
  where the $\xi_i$ are the coordinates on the unit cell and the $x_i$ are
  the coordinates on the real cell.
  This is the form needed to compute the gradient on the real cell from
  the gradient on the unit cell. If we want to transform the area element
  $dx dy$ from the real to the unit cell, we have to take the determinant of
  the inverse matrix, which is the reciprocal value of the determinant of the
  matrix defined above.

  The #FEValues# object keeps track of those fields which really need to
  be computed, since the computation of the gradients of the ansatz functions
  on each real cell can be quite an expensive thing if it is not needed. The
  object knows about which fields are needed by the #UpdateFlags# object
  passed through the constructor. In debug mode, the accessor functions, which
  return values from the different fields, check whether the required field
  was initialized, thus avoiding use of unitialized data.

  @author Wolfgang Bangerth, 1998  
  */
template <int dim>
class FEValues : public FEValuesBase<dim> {
  public:

    
    
				     /**
				      * Constructor. Fill all arrays with the
				      * values of the shape functions of the
				      * specified finite element using the
				      * quadrature points of the given
				      * quadrature rule.
				      *
				      * This function actually only fills
				      * the fields related to the unit face,
				      * the fields related to a real face (like
				      * gradients, true quadrature points, etc.)
				      * need to be initialized using the
				      * #reinit# function.
				      */
    FEValues (const FiniteElement<dim> &,
	      const Quadrature<dim> &,
	      const UpdateFlags);
    
				     /**
				      * Reinitialize the gradients, Jacobi
				      * determinants, etc for the given cell
				      * and the given finite element.
				      *
				      * This function needs a boundary object
				      * passed, since this class needs to know
				      * how to handle faces which are located
				      * on the boundary of the domain. In that
				      * case, faces may be curved and the
				      * calculation of quadrature points,
				      * gradients and the like may need
				      * additional effort, depending on the
				      * mapping from the unit to the real cell
				      * (linear mappings use straight boundary
				      * segments, but higher order elements
				      * may use other ways.)
				      */
    void reinit (const typename DoFHandler<dim>::cell_iterator &,
		 const FiniteElement<dim> &,
		 const Boundary<dim> &);

  private:
				     /**
				      * Store the gradients of the shape
				      * functions at the quadrature points on
				      * the unit cell.
				      * This field is set up upon construction
				      * of the object and contains the gradients
				      * on the reference element.
				      */
    vector<vector<Point<dim> > >   unit_shape_gradients;
    

				     /**
				      * Array of quadrature points in the unit
				      * cell. This array is set up upon
				      * construction and contains the quadrature
				      * points on the reference element.
				      */
    vector<Point<dim> >  unit_quadrature_points;
    
};




/**
  Represent a finite element evaluated with a specific quadrature rule on
  the face of a cell.

  The unit face is defined to be the tensor product of the interval $[0,1]$
  in the present number of dimensions minus one. In part of the literature,
  the convention is used that the unit cell be the tensor product of the
  interval $[-1,1]$, which is to distinguished properly.

  This class is very similar to the #FEValues# class; see there for more
  documentation. It is, however, a bit more involved: since we want to
  compute the restriction of finite element functions (here: the basis
  functions, but a finite element function is obtained by multiplication
  with the nodal values and summation) to the face of a cell and since
  finite element functions and especially their gradients need not be
  continuous at faces, we can not compute the wanted information from
  the face and a finite element class on the unit cell alone, but we
  need the real cell as well. In addition, we need to know what number
  the face is in the set of faces of the cell we want to restrict.
  Finally, since we may want to use higher order elements with unit cell
  to real cell mappings of higher than first order, thus applying curved
  boundaries, we need to know an object describing the boundary of the
  domain.


  {\bf Technical issues}

  Just like in the #FEValues# class, function values and gradients on the unit
  cell are evaluated at the quadrature points only once, in the constructor.
  Being a tensor of rank zero, the function values remain the same when we
  want them at the quadrature points on the real cell, while we get the
  gradients (a tensor of rank one) by multiplication with the Jacobi matrix
  of the transformation, which we need to compute for each cell and each
  quadrature point.

  However, while in the #FEValues# class the quadrature points are always the
  same, here we deal with more than one face. We therefore store the values
  and gradients of the ansatz functions on the unit cell in an array with as
  many elements as there are faces on a cell. The same applies for the
  quadrature points on the faces: for each face we store the position on the
  cell. This way we still need to evaluate unit gradients and function values
  only once.

  When the reinit function is called, only those gradients, quadrature points
  etc are transformed to the real cell which belong to the selected face. The
  number of the selected face is stored such that the #shape_value# function
  can return the shape function's values on the face which was last selected
  by a call to the #reinit# function.

  In addition to the complications described above, we need two different
  Jacobi matrices and determinant in this context: one for the transformation
  of the unit cell to the real cell (this Jacobi matrix is needed to
  compute the restriction of the real gradient to the given face) and one
  for the transformation of the unit face to the real face (needed to
  compute the weight factors for integration along faces). These two
  concepts have to be carefully separated.

  Finally, we will often need the outward normal to a cell at the quadrature
  points. While this could in principle be easily done using the Jacobi
  matrices at the quadrature points and the normal vectors to the unit cell
  (also easily derived, since they have an appealingly easy form for the unit
  cell ;-), it is more efficiently done by the finite element class itself.
  For example for (bi-, tri-)linear mappings the normal vector is readily
  available without complicated matrix-vector-multiplications.
  */
template <int dim>
class FEFaceValues : public FEValuesBase<dim> {
  public:
				     /**
				      * Constructor. Fill all arrays with the
				      * values of the shape functions of the
				      * specified finite element using the
				      * quadrature points of the given
				      * quadrature rule for the face, which
				      * has a dimension one less than the
				      * cell.
				      *
				      * This function actually only fills
				      * the fields related to the unit face,
				      * the fields related to a real face (like
				      * gradients, true quadrature points, etc.)
				      * need to be initialized using the
				      * #reinit# function.
				      */
    FEFaceValues (const FiniteElement<dim> &,
		  const Quadrature<dim-1> &,
		  const UpdateFlags);

				     /**
				      * Return the outward normal vector to
				      * the cell at the #i#th quadrature
				      * point. The length of the vector
				      * is normalized to one.
				      */
    const Point<dim> & normal_vector (const unsigned int i) const;
    
				     /**
				      * Return the list of outward normal
				      * vectors to the cell at the
				      * quadrature points.
				      */
    const vector<Point<dim> > & get_normal_vectors () const;
    
				     /**
				      * Reinitialize the gradients, Jacobi
				      * determinants, etc for the face with
				      * number #face_no# of #cell#
				      * and the given finite element.
				      *
				      * The constructor needs a boundary object
				      * passed, since this class needs to know
				      * how to handle faces which are located
				      * on the boundary of the domain. In that
				      * case, faces may be curved and the
				      * calculation of quadrature points,
				      * gradients and the like may need
				      * additional effort, depending on the
				      * mapping from the unit to the real cell
				      * (linear mappings use straight boundary
				      * segments, but higher order elements
				      * may use other ways.)
				      */
    void reinit (const typename DoFHandler<dim>::cell_iterator &cell,
		 const unsigned int                    face_no,
		 const FiniteElement<dim>             &fe,
		 const Boundary<dim>                  &boundary);

  private:
				     /**
				      * Store the gradients of the shape
				      * functions at the quadrature points on
				      * the unit cell.
				      * This field is set up upon construction
				      * of the object and contains the gradients
				      * on the reference element.
				      *
				      * There is one element for each face.
				      */
    vector<vector<Point<dim> > >   unit_shape_gradients[2*dim];
    

				     /**
				      * Array of quadrature points on the
				      * unit face. This is a copy of the
				      * alike field of the quadrature formula
				      * passed upon construction.
				      */
    vector<Point<dim-1> > unit_quadrature_points;
    
				     /**
				      * Array of quadrature points in the unit
				      * cell. This array is set up upon
				      * construction and contains the quadrature
				      * points on the reference element.
				      *
				      * There is one element for each face. The
				      * points are computed from those on the
				      * unit face, but are stored as coordinates
				      * on the unit cell.
				      */
    vector<Point<dim> >  global_unit_quadrature_points[2*dim];
    
				     /**
				      * List of values denoting the determinant
				      * of the transformation from the unit face
				      * to the real face. Needed to actually
				      * compute the JxW values.
				      */
    vector<double>       face_jacobi_determinants;

				     /**
				      * List of outward normal vectors at the
				      * quadrature points. This field is filled
				      * in by the finite element class.
				      */
    vector<Point<dim> >  normal_vectors;
};





/*------------------------ Inline functions: FEValuesBase ------------------------*/




template <int dim>
inline
const dFMatrix & FEValuesBase<dim>::get_shape_values () const {
  Assert (selected_dataset<shape_values.size(),
	  ExcInvalidIndex (selected_dataset, shape_values.size()));
  return shape_values[selected_dataset];
};



template <int dim>
inline
const vector<vector<Point<dim> > > &
FEValuesBase<dim>::get_shape_grads () const {
  Assert (update_flags & update_gradients, ExcAccessToUninitializedField());
  return shape_gradients;
};



template <int dim>
inline
const vector<Point<dim> > &
FEValuesBase<dim>::get_quadrature_points () const {
  Assert (update_flags & update_q_points, ExcAccessToUninitializedField());
  return quadrature_points;
};



template <int dim>
inline
const vector<Point<dim> > &
FEValuesBase<dim>::get_ansatz_points () const {
  Assert (update_flags & update_ansatz_points, ExcAccessToUninitializedField());
  return ansatz_points;
};



template <int dim>
inline
const vector<double> &
FEValuesBase<dim>::get_JxW_values () const {
  Assert (update_flags & update_JxW_values, ExcAccessToUninitializedField());
  return JxW_values;
};




/*------------------------ Inline functions: FEFaceValues ------------------------*/


template <int dim>
inline
const vector<Point<dim> > &
FEFaceValues<dim>::get_normal_vectors () const {
  Assert (update_flags & update_normal_vectors, ExcAccessToUninitializedField());
  return normal_vectors;
};




/*----------------------------   fe_values.h     ---------------------------*/
/* end of #ifndef __fe_values_H */
#endif
/*----------------------------   fe_values.h     ---------------------------*/
