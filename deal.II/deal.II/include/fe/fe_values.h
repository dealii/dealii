/*----------------------------   fe_values.h     ---------------------------*/
/*      $Id$                 */
/*      Copyright W. Bangerth, University of Heidelberg, 1998 */
#ifndef __fe_values_H
#define __fe_values_H
/*---------------------------   fe_values.h     ---------------------------*/


#include <base/exceptions.h>
#include <base/subscriptor.h>
#include <lac/full_matrix.h>
#include <dofs/dof_handler.h>
#include <base/point.h>
#include <base/tensor.h>
#include <grid/tria.h>
#include <fe/fe_update_flags.h>




/**
 *  This class offers a multitude of arrays and other fields which are used by
 *  the derived classes #FEValues#, #FEFaceValues# and #FESubfaceValues#.
 *  In principle, it is the
 *  back end of the front end for the unification of a certain finite element
 *  and a quadrature formula which evaluates certain aspects of the finite
 *  element at quadrature points.
 *  
 *  This class is an optimization which avoids evaluating the shape functions
 *  at the quadrature points each time a quadrature takes place. Rather, the
 *  values and gradients (and possibly higher order derivatives in future
 *  versions of this library) are evaluated once and for all on the unit
 *  cell or face before doing the quadrature itself. Only the Jacobian matrix of
 *  the transformation from the unit cell or face to the real cell or face and
 *  the integration points in real space are calculated each time we move on
 *  to a new face.
 *
 *  Actually, this class does none of the evaluations at startup itself; this is
 *  all done by the derived classes. It only offers the basic functionality,
 *  like providing those fields that are common to the derived classes and
 *  access to these fields. Any computations are in the derived classes. See there
 *  for more information.
 *
 *  It has support for the restriction of finite elements to faces of cells or
 *  even to subfaces (i.e. refined faces). For this purpose, it offers an array
 *  of matrices of trial function values, rather than one. Since the value of
 *  a function at a quadrature point is an invariant under the transformation
 *  from the unit cell to the real cell, it is only evaluated once upon startup.
 *  However, when considering the restriction of a finite element to a face of
 *  a cell (using a given quadrature rule), we may be tempted to compute the
 *  restriction to all faces at startup (thus ending in four array of trial
 *  function values in two dimensions, one per face, and even more in higher
 *  dimensions) and let the respective #reinit# function of the derived classes
 *  set a number which of the fields is to be taken when the user requests the
 *  function values. This is done through the #selected_dataset# variable. See
 *  the derived classes and the #get_values# function for the exact usage of
 *  this variable.
 *
 *  For many of the actual computations done by the #fill_fe_*# functions of
 *  the #FiniteElement# class and its decendants, the values and gradients of
 *  the transformation functions are needed. For example, for the computation
 *  of the real location of a quadrature point from the location on the unit
 *  cell, the values are needed, while for the computation of the Jacobian
 *  matrix the gradient is needed. While for linear elements the transformation
 *  functions coincide with the trial functions, this does not hold for higher
 *  order elements with subparametric mappings and for other types of elements
 *  such as non-conforming ones, etc, such that the precomputed values and
 *  gradients of the trial functions (#unit_shape_values# and
 *  #unit_shape_grads#) cannot be used for the present purpose.
 *  In principle, these values could be computed each time the #fill_fe_*#
 *  function is called; however, this computation is highly redundant, since
 *  only the values on the unit cell and only at the quadrature points are
 *  needed, i.e. they are the same for each cell that #fill_fe_*# is called.
 *  Therefore, two additional arrays, #unit_shape_values_transform# and
 *  #unit_shape_grads_transform# are provided, which are filled upon construction
 *  of an object of this type, which the actual finite element may or may not
 *  use. Later on, the #fill_fe_*# functions are passed pointers to these
 *  arrays, which they may use to extract the values and gradients of the
 *  transform functions. If a concrete finite element choses not to use this
 *  field, it shall set its field #transform_functions# to zero.
 *
 *  The #unit_shape_grads_transform# array is provided by the derived classes
 *  to allow for the inclusion of multiple faces, etc.
 *
 *
 *  \subsection{Definitions}
 *
 *  The Jacobian matrix is defined to be
 *  $$ J_{ij} = {d\xi_i \over dx_j} $$
 *  where the $\xi_i$ are the coordinates on the unit cell and the $x_i$ are
 *  the coordinates on the real cell.
 *  This is the form needed to compute the gradient on the real cell from
 *  the gradient on the unit cell. If we want to transform the area element
 *  $dx dy$ from the real to the unit cell, we have to take the determinant of
 *  the inverse matrix, which is the reciprocal value of the determinant of the
 *  matrix defined above.
 *
 *  The Jacobi matrix is always that of the transformation of unit to real cell.
 *  This applies also to the case where the derived class handles faces or
 *  subfaces, in which case also the transformation of unit to real cell is
 *  needed. However, the Jacobi matrix of the full transformation is always
 *  needed if we want to get the values of the gradients, which need to be
 *  transformed with the full Jacobi matrix, while we only need the
 *  transformation from unit to real face to compute the determinant of the
 *  Jacobi matrix to get the scaling of the surface element $do$.
 *
 *  The question whether to compute the Jacobi matrix as the inverse of another
 *  matrix M (which we can compute from the transformation, while we can't do
 *  so for the Jacobi matrix itself) or its transpose is a bit delicate. It
 *  should be kept in mind that when we compute the gradients in real space
 *  from those on the unit cell, we multiply with the Jacobi matrix
 *  \textit{from the right}; the whole situation is a bit confusing and it
 *  either takes deep though or trial-and-error to do it right. Some more
 *  information on this can be found in the source code documentation for the
 *  #FEQ1Mapping<dim>::fill_fe_values# function, where also a small test
 *  program is presented.
 *
 *  The derivatives of the Jacobi matrices at the quadrature points with respect
 *  to unit cell coordinates is stored in the field names
 *  #jacobi_matrices_grad#. Since the gradient of a shape function is given by
 *  $\partial_i \phi = \sum_k  \hat\partial_k \hat\phi  J_{ki}$, where
 *  $\hat\partial$ denotes differentiation on the unit cell, the second
 *  derivative of a function is given by
 *  $\partial_j \partial i \phi
 *   =
 *   \hat\partial_l [ (\hat \partial_k \hat\phi) J_{ki} ] J_{lj}
 *   =
 *   (\hat\partial_k \hat\partial_l \hat\phi) J_{ki} J_{lj}
 *   +
 *   (\hat \partial_l \hat\phi) (\hat\partial_l J_{ki}) J_{lj}$.
 *  While we already have access to the Jacobian matrix, the derivatives are
 *  stored in the named field.
 *
 *  
 *  \subsection{Member functions}
 *
 *  The functions of this class fall into different cathegories:
 *  \begin{itemize}
 *  \item #shape_value#, #shape_grad#, etc: return one of the values
 *    of this object at a time. In many cases you will want to get
 *    a whole bunch at a time for performance or convenience reasons,
 *    then use the #get_*# functions.
 *   
 *  \item #get_shape_values#, #get_shape_grads#, etc: these return
 *    a reference to a whole field. Usually these fields contain
 *    the values of all trial functions at all quadrature points.
 *
 *  \item #get_function_values#, #get_function_grads#, #...#: these
 *    functions offer a simple way to avoid the detour of the
 *    trial functions, if you have a finite element solution (resp. the
 *    vector of values associated with the different trial functions.)
 *    Then you may want to get information from the restriction of
 *    the finite element function to a certain cell, e.g. the values
 *    of the function at the quadrature points or the values of its
 *    gradient. These two functions provide the information needed:
 *    you pass it a vector holding the finite element solution and the
 *    functions return the values or gradients of the finite element
 *    function restricted to the cell which was given last time the
 *    #reinit# function was given. The same applies for the functions
 *    returning higher derivatives of the solution.
 *   
 *    Though possible in principle, these functions do not call the
 *    #reinit# function, you have to do so yourself beforehand. On the
 *    other hand, a copy of the cell iterator is stored which was used
 *    last time the #reinit# function was called. This frees us from
 *    the need to pass the cell iterator again to these two functions,
 *    which guarantees that the cell used here is in sync with that used
 *    for the #reinit# function. You should, however, make sure that
 *    nothing substantial happens to the #DoFHandler# object or any
 *    other involved instance between the #reinit# and the #get_function_*#
 *    functions are called.
 *
 *  \item #reinit#: initialize the #FEValues# object for a certain cell.
 *    This function is not in the present class but only in the derived
 *    classes and has a variable call syntax. 
 *    See the docs for the derived classes for more information.
 * \end{itemize}
 *
 *
 * \subsection{Implementational issues}
 *
 * The #FEValues# object keeps track of those fields which really need to
 * be computed, since the computation of the gradients of the trial functions
 * and of other values on each real cell can be quite an expensive thing
 * if it is not needed. The
 * object knows about which fields are needed by the #UpdateFlags# object
 * passed through the constructor. In debug mode, the accessor functions, which
 * return values from the different fields, check whether the required field
 * was initialized, thus avoiding use of unitialized data.
 *
 * Functions should not assume that one flag is needed for another object as
 * well. For example, the computation of the Jacobi determinant usually
 * requires the computation of the Jacobi matrix. However, functions shall
 * not assume that someone who wants to get the #JxW_values# must set the
 * #update_jacobians# flag besides the #update_JxW_values# flag.
 *
 * It is also forbidden that the constructor of a class set the
 * #update_jacobians# flag if the user specifies #update_JxW_values#. This is
 * since derived classes may be able to compute the #JxW_values# field without
 * the Jacobian matrices, but need to do the latter since they can't know who
 * set the #update_jacobians# flag.
 *
 * @author Wolfgang Bangerth, 1998
 */
template <int dim>
class FEValuesBase
{
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
    const unsigned int dofs_per_cell;

				     /**
				      * Number of basis functions for the
				      * transformation from the unit cell
				      * to the real cell. See the docs for
				      * more information on this field.
				      */
    const unsigned int n_transform_functions;
    
				     /**
				      * Constructor. Set up the array sizes
				      * with #n_q_points# quadrature points,
				      * #n_support_points# support points (on
				      * the cell or face), #n_dof# trial
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
		  const unsigned int n_support_points,
		  const unsigned int dofs_per_cell,
		  const unsigned int n_transform_functions,
		  const unsigned int n_values_array,
		  const UpdateFlags         update_flags,
		  const FiniteElement<dim> &fe);
    

				     /**
				      * Return the value of the #i#th shape
				      * function at the #j# quadrature point
				      * on the cell, face or subface selected
				      * the last time the #reinit# function
				      * of the derived class was called.
				      */
    double shape_value (const unsigned int function,
			const unsigned int quadrature_point) const;

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
    const FullMatrix<double> & get_shape_values () const;

				     /**
				      * Return the values of the finite
				      * element function characterized
				      * by #fe_function# restricted to
				      * the cell, face or subface selected
				      * the last time the #reinit# function
				      * of the derived class was called,
				      * at the quadrature points.
				      *
				      * If the present cell is not an active
				      * one the interpolated function values
				      * are returned.
				      *
				      * To get values of
				      * multi-component elements,
				      * there is another
				      * #get_function_values#
				      * returning a vector of vectors of
				      * results.
				      *
				      * The function assumes that the
				      * #values# object already has the
				      * right size. 
				      */
    void get_function_values (const Vector<double> &fe_function,
			      vector<double>       &values) const;

				     /**
				      * Access to vector valued finite element functions.
				      *
				      * This function does the same as
				      * the other #get_function_values#,
				      * but applied to multi-component
				      * elements.
				      */
    void get_function_values (const Vector<double>    &fe_function,
			      vector<Vector<double> > &values) const;

    				     /**
				      * Return the gradient of the #i#th shape
				      * function at the #j# quadrature point.
				      * If you want to get the derivative in
				      * one of the coordinate directions, use
				      * the appropriate function of the #Tensor#
				      * class to extract one component. Since
				      * only a reference to the gradient's value
				      * is returned, there should be no major
				      * performance drawback.
				      * The function returns the gradient on the
				      * real element, not the reference element.
				      */
    const Tensor<1,dim> & shape_grad (const unsigned int function,
				      const unsigned int quadrature_point) const;

				     /** 
				      * Return a pointer to the matrix holding
				      * all gradients of shape functions at all
				      * integration points, on the present cell.
				      * For the format of this matrix, see the
				      * documentation for the matrix itself.
				      */
    const vector<vector<Tensor<1,dim> > > & get_shape_grads () const;

				     /**
				      * Return the gradients of the finite
				      * element function characterized
				      * by #fe_function# restricted to
				      * #cell# at the quadrature points.
				      *
				      * If the present cell is not an active
				      * one the interpolated function values
				      * are returned.
				      *
				      * The function assumes that the
				      * #gradients# object already has the
				      * right size.
				      */
    void get_function_grads (const Vector<double>   &fe_function,
			     vector<Tensor<1,dim> > &gradients) const;

				     /**
				      * Return the gradients of the finite
				      * element function characterized
				      * by #fe_function# restricted to
				      * #cell# at the quadrature points.
				      *
				      * If the present cell is not an active
				      * one the interpolated function values
				      * are returned.
				      *
				      * The function assumes that the
				      * #gradients# object already has the
				      * right size.
				      *
				      * This function does the same as
				      * the other #get_function_values#,
				      * but applied to multi-component
				      * elements.
				      */
    void get_function_grads (const Vector<double>   &fe_function,
			     vector<vector<Tensor<1,dim> > > &gradients) const;

    				     /**
				      * Return the 2nd derivatives of
				      * the #i#th shape function at
				      * the #j# quadrature point. If
				      * you want to get the derivatives
				      * in one of the coordinate
				      * directions, use the
				      * appropriate function of the
				      * #Tensor# class to extract one
				      * component. Since only a
				      * reference to the derivatives'
				      * values is returned, there
				      * should be no major performance
				      * drawback.  The function
				      * returns the derivatives on the
				      * real element, not the
				      * reference element.
				      */
    const Tensor<2,dim> & shape_2nd_derivative (const unsigned int function,
						const unsigned int quadrature_point) const;

				     /**
				      * Return a pointer to the
				      * matrix holding all 2nd
				      * derivatives of shape functions
				      * at all integration points, on
				      * the present cell.  For the
				      * format of this matrix, see the
				      * documentation for the matrix
				      * itself.
				      */
    const vector<vector<Tensor<2,dim> > > & get_shape_2nd_derivatives () const;
    
				     /**
				      * Return the tensor of second
				      * derivatives of the finite
				      * element function characterized
				      * by #fe_function# restricted to
				      * #cell# at the quadrature points.
				      *
				      * The function assumes that the
				      * #second_derivatives# object already has
				      * the right size.
				      */
    void get_function_2nd_derivatives (const Vector<double>   &fe_function,
				       vector<Tensor<2,dim> > &second_derivatives) const;

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
				      * the #i#th trial function is located
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
				      * function's value at the trial function
				      * locations.
				      *
				      * For the evaluation of finite elements on
				      * faces of cells, #i# is the number
				      * of the trial function on the face, not
				      * on the cell.
				      */
    const Point<dim> & support_point (const unsigned int i) const;

				     /**
				      * Return a pointer to the vector of points
				      * denoting the location of the trial
				      * functions.
				      */
    const vector<Point<dim> > & get_support_points () const;
    
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
    double JxW (const unsigned int quadrature_point) const;

				     /**
				      * Return a pointer to the array holding
				      * the Jacobi determinant times the
				      * quadrature weight at the different
				      * quadrature points.
				      */
    const vector<double> & get_JxW_values () const;

				     /**
				      * Return a constant reference to the
				      * selected finite element object. This
				      * function is inline, so it should
				      * be reasonably fast.
				      */
    const FiniteElement<dim> & get_fe () const;

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
    DeclException0 (ExcWrongNoOfComponents);
				     /**
				      * Exception.
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
    DeclException0 (ExcInvalidUpdateFlag);
    
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
				      * restriction of the finite element trial
				      * space to a cell. For face values, the
				      * vector contains as many elements as
				      * there are faces, for subfaces the same
				      * applies. Which of the matrices is active
				      * is determined by the #selected_dataset#
				      * variable.
				      */
    vector<FullMatrix<double> >     shape_values;

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
    vector<vector<Tensor<1,dim> > >  shape_gradients;

				     /**
				      * Store the 2nd derivatives of the shape
				      * functions at the quadrature points.
				      *
				      * This field is reset each time
				      * #reinit# is called and contains the
				      * gradients on the real element, rather
				      * than on the reference element.
				      */
    vector<vector<Tensor<2,dim> > >  shape_2nd_derivatives;

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
				      * of the trial functions. In real space
				      * (no-one seems to need the off-point
				      * on the unit cell, so no function is
				      * provided for this).
				      */
    vector<Point<dim> >  support_points;
    
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
    vector<Tensor<2,dim> > jacobi_matrices;

				     /**
				      * Store the derivatives of the jacobi
				      * matrices. If #J[j][k]# is the jacobi
				      * matrix, then the index #[i][j][k]#
				      * of this field denotes the derivation
				      * of #J[j][k]# with respect to the
				      * #i#th variable.
				      *
				      * The same general remarks apply as for
				      * #jacobi_matrices#.
				      */
    vector<Tensor<3,dim> > jacobi_matrices_grad;
    
				     /**
				      * Store the values of the basis functions
				      * of the transformation from unit cell
				      * to real cell at the quadrature points.
				      *
				      * This field stores some data which is not
				      * really needed for the assemblage of
				      * matrices and vectors but makes that
				      * operation much faster. Each time the
				      * #FEValues::reinit# function calls
				      * the #FiniteElemenet::fill_fe_values#
				      * function, this and the next array are
				      * passed. The #fill_fe_values# function
				      * may or may not use this field.
				      *
				      * The element #(i,j)# denotes the value
				      * of the #i#th transfer basis function
				      * at the #j#th quadrature point.
				      */
    vector<FullMatrix<double> >     shape_values_transform;

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
    UpdateFlags          update_flags;

				     /**
				      * Store the cell selected last time
				      * the #reinit# function was called
				      * to make access
				      * to the #get_function_*# functions
				      * safer.
				      */
    DoFHandler<dim>::cell_iterator present_cell;

				     /**
				      * Store the finite element for later use.
				      */
    const SmartPointer<const FiniteElement<dim> > fe;
};




/**
 * Represent a finite element evaluated with a specific quadrature rule on
 * a cell.
 *
 * The unit cell is defined to be the tensor product of the interval $[0,1]$
 * in the present number of dimensions. In part of the literature, the convention
 * is used that the unit cell be the tensor product of the interval $[-1,1]$,
 * which is to distinguished properly.
 *
 * Objects of this class store a multitude of different values needed to
 * do the assemblage steps on real cells rather than on the unit cell. Among
 * these values are the values and gradients of the shape functions at the
 * quadrature points on the real and the unit cell, the location of the
 * quadrature points on the real and on the unit cell, the weights of the
 * quadrature points, the Jacobian matrices of the mapping from the unit to
 * the real cell at the quadrature points and so on.
 *
 * @author Wolfgang Bangerth, 1998  
 */
template <int dim>
class FEValues : public FEValuesBase<dim>
{
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
				      *
				      * UPDATE!
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
    FEValues (const FiniteElement<dim> &fe,
	      const Quadrature<dim>    &quadrature,
	      const UpdateFlags         update_flags);
    
				     /**
				      * Reinitialize the gradients, Jacobi
				      * determinants, etc for the given cell
				      * and the given finite element.
				      */
    void reinit (const typename DoFHandler<dim>::cell_iterator &);

  private:
				     /**
				      * Store the gradients of the shape
				      * functions at the quadrature points on
				      * the unit cell.
				      * This field is set up upon construction
				      * of the object and contains the gradients
				      * on the reference element.
				      */
    vector<vector<Tensor<1,dim> > > unit_shape_gradients;

				     /**
				      * Store the 2nd derivatives of the shape
				      * functions at the quadrature points on
				      * the unit cell.
				      * This field is set up upon construction
				      * of the object and contains the
				      * derivatives on the reference element.
				      */
    vector<vector<Tensor<2,dim> > > unit_shape_2nd_derivatives;
    
				     /**
				      * Gradients of the basis
				      * functions of the transformation.
				      * Analogous to the #shape_values_transform#
				      * array of the base class.
				      */
    vector<vector<Tensor<1,dim> > > unit_shape_gradients_transform;
    
				     /**
				      * Array of quadrature points in the unit
				      * cell. This array is set up upon
				      * construction and contains the quadrature
				      * points on the reference element.
				      */
    vector<Point<dim> >  unit_quadrature_points;    
};





/**
 *  This class provides for the data elements needed for the restriction of
 *  finite elements to faces or subfaces. It does no real computations, apart
 *  from initialization of the fields with the right size. It more or
 *  less is only a base class to the #FEFaceValues# and #FESubfaceValues#
 *  classes which do the real computations. See there for descriptions of
 *  what is really going on.
 *
 *  Since many of the concepts are the same whether we restrict a finite element
 *  to a face or a subface (i.e. the child of the face of a cell), we describe
 *  those common concepts here, rather than in the derived classes.
 *
 *
 *  \subsection{Technical issues}
 *
 *  The unit face is defined to be the tensor product of the interval $[0,1]$
 *  in the present number of dimensions minus one. In part of the literature,
 *  the convention is used that the unit cell/face be the tensor product of the
 *  interval $[-1,1]$, which is to distinguished properly. A subface is the
 *  child of a face; they are numbered in the way laid down in the
 *  #Triangulation# class.
 *
 *  Just like in the #FEValues# class, function values and gradients on the unit
 *  face or subface are evaluated at the quadrature points only once, and stored
 *  by the common base class. Being a tensor of rank zero, the function values
 *  remain the same when we want them at the quadrature points on the real cell,
 *  while we get the gradients (a tensor of rank one) by multiplication with the
 *  Jacobi matrix of the transformation, which we need to compute for each cell
 *  and each quadrature point.
 *
 *  However, while in the #FEValues# class the quadrature points are always the
 *  same, here we deal with more than one (sub)face. We therefore store the values
 *  and gradients of the trial functions on the unit cell in an array with as
 *  many elements as there are (sub)faces on a cell. The same applies for the
 *  quadrature points on the (sub)faces: for each (sub)face we store the position
 *  on the cell. This way we still need to evaluate unit gradients and function
 *  values only once and only recompute the gradients on the real (sub)face by
 *  multiplication of the unit gradients on the presently selected (sub)face
 *  with the Jacobi matrix.
 *
 *  
 *  When the #reinit# function of a derived class is called, only those
 *  gradients, quadrature points etc are transformed to the real cell which
 *  belong to the selected face or subface. The number of the selected face
 *  or subface is stored in the #selected_dataset# variable of the base class
 *  such that the #shape_value# function can return the shape function's
 *  values on the (sub)face which was last selected by a call to the #reinit#
 *  function.
 *
 *  In addition to the complications described above, we need two different
 *  Jacobi matrices and determinants in this context: one for the transformation
 *  of the unit cell to the real cell (this Jacobi matrix is needed to
 *  compute the restriction of the real gradient to the given face) and one
 *  for the transformation of the unit face to the real face or subface
 *  (needed to compute the weight factors for integration along faces). These two
 *  concepts have to be carefully separated.
 *
 *  Finally, we will often need the outward normal to a cell at the quadrature
 *  points. While this could in principle be easily done using the Jacobi
 *  matrices at the quadrature points and the normal vectors to the unit cell
 *  (also easily derived, since they have an appealingly simple form for the unit
 *  cell ;-), it is more efficiently done by the finite element class itself.
 *  For example for (bi-, tri-)linear mappings the normal vector is readily
 *  available without complicated matrix-vector-multiplications.
 *    
 *  @author Wolfgang Bangerth, 1998
*/
template <int dim>
class FEFaceValuesBase : public FEValuesBase<dim>
{
  public:
				     /**
				      * Constructor. Call the constructor of
				      * the base class and set up the arrays
				      * of this class with the right sizes.
				      * Actually filling these arrays is a
				      * duty of the derived class's
				      * constructors.
				      *
				      * #n_faces_or_subfaces# is the number
				      * of faces or subfaces that this object
				      * is to store. The actual number depends
				      * on the derived class, for
				      * #FEFaceValues# it is #2*dim#, while for
				      * the #FESubfaceValues# class it is
				      * #2*dim*(1<<(dim-1))#, i.e. the number
				      * of faces times the number of subfaces
				      * per face.
				      */
    FEFaceValuesBase (const unsigned int n_q_points,
		      const unsigned int n_support_points,
		      const unsigned int dofs_per_cell,
		      const unsigned int n_transform_functions,
		      const unsigned int n_faces_or_subfaces,
		      const UpdateFlags         update_flags,
		      const FiniteElement<dim> &fe);

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

  protected:
				     /**
				      * Store the gradients of the shape
				      * functions at the quadrature points on
				      * the unit cell.
				      * This field is set up upon construction
				      * of the object and contains the gradients
				      * on the reference element.
				      *
				      * There is one element for each face or
				      * subface, with indices like that:
				      * #unit_shape_gradients[face][dof][q_point]#
				      */
    vector<vector<vector<Tensor<1,dim> > > > unit_shape_gradients;
    
				     /**
				      * Store the 2nd derivatives of the shape
				      * functions at the quadrature points on
				      * the unit cell for each face.
				      * This field is set up upon construction
				      * of the object and contains the
				      * derivatives on the reference element.
				      */
    vector<vector<vector<Tensor<2,dim> > > > unit_shape_2nd_derivatives;

				     /**
				      * Gradients of the basis
				      * functions of the transformation.
				      * Analogous to the #shape_values_transform#
				      * array of the base class.
				      */
    vector<vector<vector<Tensor<1,dim> > > > unit_shape_gradients_transform;

    				     /**
				      * Array of quadrature points on the
				      * unit face. This is a copy of the
				      * alike field of the quadrature formula
				      * passed upon construction.
				      */
    vector<Point<dim-1> > unit_face_quadrature_points;
    
				     /**
				      * Array of quadrature points in the unit
				      * cell. This array is set up upon
				      * construction and contains the quadrature
				      * points on the reference element.
				      *
				      * There is one element for each face or
				      * subface. The points are computed from
				      * those on the unit face, but are stored
				      * as coordinates on the unit cell.
				      */
    vector<vector<Point<dim> > > unit_quadrature_points;
    
				     /**
				      * List of values denoting the determinant
				      * of the transformation from the unit face
				      * to the real face or subface. Needed to
				      * actually compute the JxW values.
				      */
    vector<double>       face_jacobi_determinants;

				     /**
				      * List of outward normal vectors at the
				      * quadrature points. This field is filled
				      * in by the finite element class.
				      */
    vector<Point<dim> >  normal_vectors;
};



/**
 * Represent a finite element evaluated with a specific quadrature rule on
 * the face of a cell.
 *
 * This class is very similar to the #FEValues# class; see there for more
 * documentation. It is, however, a bit more involved: since we want to
 * compute the restriction of finite element functions (here: the basis
 * functions, but a finite element function is obtained by multiplication
 * with the nodal values and summation) to the face of a cell and since
 * finite element functions and especially their gradients need not be
 * continuous at faces, we can not compute the wanted information from
 * the face and a finite element class on the unit cell alone, but we
 * need the real cell as well. In addition, we need to know what number
 * the face is in the set of faces of the cell we want to restrict.
 * Finally, since we may want to use higher order elements with unit cell
 * to real cell mappings of higher than first order, thus applying curved
 * boundaries, we need to know an object describing the boundary of the
 * domain.
 *    
 * @author Wolfgang Bangerth, 1998
 */
template <int dim>
class FEFaceValues : public FEFaceValuesBase<dim>
{
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
				      *
				      * UPDATE!
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
    FEFaceValues (const FiniteElement<dim> &,
		  const Quadrature<dim-1> &,
		  const UpdateFlags);

				     /**
				      * Reinitialize the gradients, Jacobi
				      * determinants, etc for the face with
				      * number #face_no# of #cell#
				      * and the given finite element.
				      */
    void reinit (const typename DoFHandler<dim>::cell_iterator &cell,
		 const unsigned int                    face_no);
};




/**
 * Represent a finite element evaluated with a specific quadrature rule on
 * the child of the face of a cell.
 *
 * This class is very similar to the #FEFaceValues# class; see there for
 * more documentation. It serves the computation of interface integrals
 * where the cells on both sides of the face have different refinement
 * levels. This is useful for example when we want to integrate the jump
 * of the gradient of the finite element solution along the boundary of
 * a cell to estimate the error. Now, this is not so much of a problem
 * if all neighbors of the cell have the same refinement level, then we
 * will use the #FEFaceValues# class, but it gets trickier if one of the
 * cells is more refined than the other.
 *
 * To this end, there seem to be two ways which may be applicable:
 * \begin{itemize}
 * \item Prolong the coarser cell to the finer refinement level: we could
 *   compute the prolongation of the finite element functions to the
 *   child cells and consider the subface a face of one of the child cells.
 *   This approach seems clear and rather simple to implement, however it
 *   has two major drawbacks: first, the finite element space on the
 *   refined (child) cells may not be included in the space of the unrefined
 *   cell, in which case the prolongation would alter information and thus
 *   make computations worthless in the worst case. The second reason is
 *   a practical one, namely that by refining the cell virtually, we would
 *   end up with child cells which do not exist in real and can thus not be
 *   represented in terms of iterators. This would mean that we had to change
 *   the whole interface to the #FE*Values# classes to accept cell corner
 *   points by value, etc, instead of relying on appropriate iterators. This
 *   seems to be clumsy and not very suitable to maintain an orthogonal
 *   programming style. Apart from that, we already have iterators, why
 *   shouldn't we use them?
 *   
 * \item Use 'different' quadrature formulae: this second approach is the
 *   way we chose here. The idea is to evaluate the finite element trial
 *   functions on the two cells restricted to the face in question separately,
 *   by restricting the trial functions on the less refined cell to its
 *   face and the functions on the more refined cell to its face as well,
 *   the second face being a child to the first one. Now, if we would use
 *   the same quadrature formula for both restrictions, we would end up with
 *   the same number of quadrature points, but at different locations since
 *   they were evaluated on faces of different size. We therefore use the
 *   original quadrature formula for the refined cell and a modified one for
 *   the coarse cell, the latter being modified in such a way that the
 *   locations of the quadrature points match each other.
 *
 *   An example may shed more light onto this: assume we are in two dimension,
 *   we have a cell of which we want to evaluate a finite element function on
 *   face zero, and neighbor zero is refined (then so is face zero). The
 *   quadrature formula shall be the Simpson rule with quadrature points
 *   $0$, $0.5$ and $1$. The present cell shall be the unit cell, without
 *   loss of generality. Then the face in question is the line $(0,0)$ to
 *   $(1,0)$, subdivided into two subfaces. We will then compute the
 *   restriction of the present cell to the common subface $(0,0)$ to
 *   $(0.5,5)$ by using a modified quadrature formulae with quadrature
 *   points $(0,0)$, $(0.25,0)$ and $(0.5,0)$ (coordinates on the cell)
 *   which is not symmetric as was the original quadrature rule for a line.
 *   This modified quadrature rule is computed by projection onto the subface
 *   using the #QProjector<dim>::project_to_subface()# function. The neighboring
 *   cell, being refined once more than the present is evaluated with the
 *   quadrature formula projected to the common face, but using the original
 *   quadrature formula. This way, the locations of the quadrature points
 *   on both sides of the common face match each other.
 * \end{itemize}
 *
 * For a use of this mechanism, take a look of the code in the error
 * estimation hierarchy, since there often the jump of a finite element
 * function's gradient across cell boundaries is computed.
 *
 *
 * \subsection{Other implementational subjects}
 *
 * It does not seem useful to ask for the off-points of the trial functions
 * (name #support_points# in the #FEValuesBase# class) for subfaces. These are
 * therefore not supported for this class and should throw an error if
 * accessed. Specifying #update_support_points# for the #UpdateFlags# in the
 * constructor is disallowed.
 *
 * The values of the trial functions on the subfaces are stored as an array
 * of matrices, each matrix representing the values of the trial functions at
 * the quadrature points at one subface. The ordering is as follows: the values
 * of the trial functions at face #face#, subface #subface# are stored in
 * #shape_values[face*(1<<(dim-1))+subface]#. The same order applies for the
 * quadrature points on the unit cell, which are stored in the
 * #unit_quadrature_points# array. Note that #1<<(dim-1)# is the number of
 * subfaces per face.
 *
 * One subtle problem is that if a face is at the boundary, then computation
 * of subfaces may be a bit tricky, since we do not know whether the user
 * intends to better approximate the boundary by the subfaces or only wants
 * to have the subfaces be one part of the mother face. However, it is hardly
 * conceivable what someone wants when using this class for faces at the
 * boundary, in the end this class was invented to facilitate integration
 * along faces with cells of different refinement levels on both sides,
 * integration along the boundary of the domain is better done through
 * the #FEFaceValues# class. For this reason, calling #reinit# with a
 * boundary face will result in an error.
 * 
 * @author Wolfgang Bangerth, 1998
 */
template <int dim>
class FESubfaceValues : public FEFaceValuesBase<dim>
{
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
    FESubfaceValues (const FiniteElement<dim> &fe,
		     const Quadrature<dim-1>  &face_quadrature,
		     const UpdateFlags         update_flags);

				     /**
				      * Reinitialize the gradients, Jacobi
				      * determinants, etc for the face with
				      * number #face_no# of #cell#
				      * and the given finite element.
				      */
    void reinit (const typename DoFHandler<dim>::cell_iterator &cell,
		 const unsigned int                    face_no,
		 const unsigned int                    subface_no);

				     /**
				      * Exception
				      */
    DeclException0 (ExcReinitCalledWithBoundaryFace);
};





/*------------------------ Inline functions: FEValuesBase ------------------------*/




template <int dim>
inline
const FullMatrix<double> & FEValuesBase<dim>::get_shape_values () const
{
  Assert (update_flags & update_values, ExcAccessToUninitializedField());
  Assert (selected_dataset<shape_values.size(),
	  ExcIndexRange (selected_dataset, 0, shape_values.size()));
  return shape_values[selected_dataset];
};



template <int dim>
inline
const vector<vector<Tensor<1,dim> > > &
FEValuesBase<dim>::get_shape_grads () const
{
  Assert (update_flags & update_gradients, ExcAccessToUninitializedField());
  return shape_gradients;
};



template <int dim>
inline
const vector<vector<Tensor<2,dim> > > &
FEValuesBase<dim>::get_shape_2nd_derivatives () const
{
  Assert (update_flags & update_second_derivatives, ExcAccessToUninitializedField());
  return shape_2nd_derivatives;
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
FEValuesBase<dim>::get_support_points () const {
  Assert (update_flags & update_support_points, ExcAccessToUninitializedField());
  return support_points;
};



template <int dim>
inline
const vector<double> &
FEValuesBase<dim>::get_JxW_values () const {
  Assert (update_flags & update_JxW_values, ExcAccessToUninitializedField());
  return JxW_values;
};



template <int dim>
inline
const FiniteElement<dim> & 
FEValuesBase<dim>::get_fe () const {
  return *fe;
};


/*------------------------ Inline functions: FEFaceValuesBase --------------------*/


template <int dim>
inline
const vector<Point<dim> > &
FEFaceValuesBase<dim>::get_normal_vectors () const {
  Assert (update_flags & update_normal_vectors, ExcAccessToUninitializedField());
  return normal_vectors;
};




/*----------------------------   fe_values.h     ---------------------------*/
/* end of #ifndef __fe_values_H */
#endif
/*----------------------------   fe_values.h     ---------------------------*/
