//---------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 1998, 1999, 2000, 2001, 2002, 2003, 2004, 2005 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------
#ifndef __deal2__fe_values_h
#define __deal2__fe_values_h


#include <base/config.h>
#include <base/exceptions.h>
#include <base/subscriptor.h>
#include <base/point.h>
#include <base/tensor.h>
#include <base/vector_slice.h>
#include <base/quadrature.h>
#include <base/table.h>
#include <lac/vector.h>
#include <lac/block_vector.h>
#include <lac/petsc_vector.h>
#include <lac/petsc_block_vector.h>
#include <grid/tria.h>
#include <grid/tria_iterator.h>
#include <dofs/dof_handler.h>
#include <dofs/dof_accessor.h>
#include <fe/fe.h>
#include <fe/fe_update_flags.h>
#include <fe/mapping.h>
#include <multigrid/mg_dof_handler.h>
#include <multigrid/mg_dof_accessor.h>

#include <algorithm>
#include <memory>

template <int dim> class Quadrature;

//TODO: Add access to mapping values to FEValuesBase
//TODO: Several FEValuesBase of a system should share Mapping

/*!@addtogroup feaccess */
/*@{*/

/**
 * Contains all data vectors for FEValues.
 * This class has been extracted from FEValuesBase to be handed
 * over to the fill functions of Mapping and
 * FiniteElement.
 *
 * @note All data fields are public, but this is not
 * critical, because access to this object is private in FEValues.
 *
 * @author Guido Kanschat, 2000
 */
template <int dim>
class FEValuesData
{
  public:
				     /**
				      * Initialize all vectors to
				      * correct size.
				      */
    void initialize (const unsigned int        n_quadrature_points,
		     const FiniteElement<dim> &fe,
		     const UpdateFlags         flags);

				     /**
				      * Storage type for shape
				      * values. Each row in the matrix
				      * denotes the values of a single
				      * shape function at the
				      * different points, columns are
				      * for a single point with the
				      * different shape functions.
				      *
				      * If a shape function has more
				      * than one non-zero component
				      * (in deal.II diction: it is
				      * non-primitive), then we
				      * allocate one row per non-zero
				      * component, and shift
				      * subsequent rows backward.
				      * Lookup of the correct row for
				      * a shape function is thus
				      * simple in case the entire
				      * finite element is primitive
				      * (i.e. all shape functions are
				      * primitive), since then the
				      * shape function number equals
				      * the row number. Otherwise, use
				      * the
				      * #shape_function_to_row_table
				      * array to get at the first row
				      * that belongs to this
				      * particular shape function, and
				      * navigate among all the rows
				      * for this shape function using
				      * the
				      * FiniteElement::get_nonzero_components()
				      * function which tells us which
				      * components are non-zero and
				      * thus have a row in the array
				      * presently under discussion.
				      */
    typedef Table<2,double> ShapeVector;

				     /**
				      * Storage type for
				      * gradients. The layout of data
				      * is the same as for the
				      * #ShapeVector data type.
				      */
    typedef std::vector<std::vector<Tensor<1,dim> > > GradientVector;

				     /**
				      * Likewise for second order
				      * derivatives.
				      */
    typedef std::vector<std::vector<Tensor<2,dim> > > GradGradVector;

				     /**
				      * Store the values of the shape
				      * functions at the quadrature
				      * points. See the description of
				      * the data type for the layout
				      * of the data in this field.
				      */
    ShapeVector shape_values;

    				     /**
				      * Store the gradients of the
				      * shape functions at the
				      * quadrature points. See the
				      * description of the data type
				      * for the layout of the data in
				      * this field.
				      */
    GradientVector shape_gradients;

				     /**
				      * Store the 2nd derivatives of
				      * the shape functions at the
				      * quadrature points.  See the
				      * description of the data type
				      * for the layout of the data in
				      * this field.
				      */
    GradGradVector shape_2nd_derivatives;

				     /**
				      * Store an array of weights
				      * times the Jacobi determinant
				      * at the quadrature points. This
				      * function is reset each time
				      * reinit() is called. The
				      * Jacobi determinant is actually
				      * the reciprocal value of the
				      * Jacobi matrices stored in this
				      * class, see the general
				      * documentation of this class
				      * for more information.
				      *
				      * However, if this object refers
				      * to an FEFaceValues of
				      * FESubfaceValues object, then
				      * the JxW_values correspond to
				      * the Jacobian of the
				      * transformation of the face,
				      * not the cell, i.e. the
				      * dimensionality is that of a
				      * surface measure, not of a
				      * volume measure. In this case,
				      * it is computed from the
				      * boundary forms, rather than
				      * the Jacobian matrix.
				      */
    std::vector<double>       JxW_values;

				     /**
				      * Array of quadrature points. This array
				      * is set up upon calling reinit() and
				      * contains the quadrature points on the
				      * real element, rather than on the
				      * reference element.
				      */
    std::vector<Point<dim> >  quadrature_points;

				     /**
				      * List of outward normal vectors at the
				      * quadrature points. This field is filled
				      * in by the finite element class.
				      */
    std::vector<Point<dim> >  normal_vectors;

                                     /**
				      * List of boundary forms at the
				      * quadrature points. This field is filled
				      * in by the finite element class.
				      */
    std::vector<Tensor<1,dim> >  boundary_forms;

				     /**
				      * Indicate the first row which a
				      * given shape function occupies
				      * in the #shape_values,
				      * #shape_gradients and
				      * #shape_2nd_derivatives
				      * arrays. If all shape functions
				      * are primitive, then this is
				      * the identity mapping. If, on
				      * the other hand some shape
				      * functions have more than one
				      * non-zero vector components,
				      * then they may occupy more than
				      * one row, and this array
				      * indicates which is the first
				      * one.
				      *
				      * The questions which particular
				      * vector component occupies
				      * which row for a given shape
				      * function is answered as
				      * follows: we allocate one row
				      * for each non-zero component as
				      * indicated by the
				      * FiniteElement::get_nonzero_components()
				      * function, and the rows are in
				      * ascending order exactly those
				      * non-zero components.
				      */
    std::vector<unsigned int> shape_function_to_row_table;
    
                                     /**
				      * Original update flags handed
				      * to the constructor of
				      * FEValues.
				      */
    UpdateFlags          update_flags;
};


/**
 * FEValues, FEFaceValues and FESubfaceValues objects are programming
 * interfaces to finite element and mapping classes on the one hand
 * side, to cells and quadrature rules on the other side. The reason
 * for their existence is possible optimization. Depending on the type
 * of finite element and mapping, some values can be computed once on
 * the unit cell. Others must be computed on each cell, but maybe
 * computation of several values at the same time offers ways for
 * optimization. Since this interlay may be complex and depends on the
 * actual finite element, it cannot be left to the applications
 * programmer.
 *
 * FEValues, FEFaceValues and FESubfaceValues provide only data
 * handling: computations are left to objects of type Mapping and
 * FiniteElement. These provide functions <tt>get_*_data</tt> and
 * <tt>fill_*_values</tt> which are called by the constructor and
 * <tt>reinit</tt> functions of <tt>FEValues*</tt>, respectively.
 *
 * <h3>General usage</h3>
 *
 * Usually, an object of <tt>FEValues*</tt> is used in integration loops
 * over all cells of a triangulation. To take full advantage of the
 * optimization features, it should be constructed before the
 * loop. Then, it must be re-initialized for each grid cell. This is
 * like a magnifying glass being used to look at one item after the
 * other. A typical piece of code looks like this:
 *
 * @code
 * FEValues values (mapping, finite_element, quadrature, flags);
 * for (cell = dof_handler.begin_active();
 *      cell != dof_handler.end();
 *      ++cell)
 *   {
 *     values.reinit(cell);
 *     ...
 *   }
 * @endcode
 *
 *
 *  <h3>Member functions</h3>
 *
 *  The functions of this class fall into different cathegories:
 *  <ul>
 *  <li> shape_value(), shape_grad(), etc: return one of the values 
 *    of this object at a time. These functions are inlined, so this
 *    is the suggested access to all finite element values. There
 *    should be no loss in performance with an optimizing compiler. If
 *    the finite element is vector valued, then these functions return
 *    the only non-zero component of the requested shape
 *    function. However, some finite elements have shape functions
 *    that have more than one non-zero component (we call them
 *    non-"primitive"), and in this case this set of functions will
 *    throw an exception since they cannot generate a useful
 *    result. Rather, use the next set of functions.
 *
 *  <li> shape_value_component(), shape_grad_component(), etc:
 *    This is the same set of functions as above, except that for vector
 *    valued finite elements they return only one vector component. This
 *    is useful for elements of which shape functions have more than one
 *    non-zero component, since then the above functions cannot be used,
 *    and you have to walk over all (or only the non-zero) components of
 *    the shape function using this set of functions.
 *   
 *  <li> get_function_values(), get_function_grads(), etc.: Compute a
 *    finite element function or its derivative in quadrature points.
 *
 *  <li> reinit: initialize the FEValues object for a certain cell.
 *    This function is not in the present class but only in the derived
 *    classes and has a variable call syntax. 
 *    See the docs for the derived classes for more information.
 * </ul>
 *
 *
 * <h3>UpdateFlags</h3>
 *
 * The UpdateFlags object handed to the constructor is used to
 * determine, which of the data fields to compute. This way, it is
 * possible to avoid expensive computations of useless derivatives.
 * In the beginning, these flags are processed through the functions
 * Mapping::update_once(), Mapping::update_each(),
 * FiniteElement::update_once() FiniteElement::update_each(). All the
 * results are bit-wise or'd and determine the fields actually
 * computed. This enables Mapping and FiniteElement to schedule
 * auxiliary data fields for updating. Still, it is recommended to
 * give <b>all</b> needed update flags to FEValues.
 *
 * @author Wolfgang Bangerth, 1998, 2003, Guido Kanschat, 2001
 */
template <int dim>
class FEValuesBase : protected FEValuesData<dim>
{
  public:
				     /**
				      * Number of quadrature points.
				      */
    const unsigned int n_quadrature_points;

				     /**
				      * Number of shape functions per
				      * cell. If we use this base
				      * class to evaluate a finite
				      * element on faces of cells,
				      * this is still the number of
				      * degrees of freedom per cell,
				      * not per face.
				      */
    const unsigned int dofs_per_cell;

    
				     /**
				      * Constructor. Set up the array
				      * sizes with <tt>n_q_points</tt>
				      * quadrature points, <tt>dofs_per_cell</tt>
				      * trial functions per cell and
				      * with the given pattern to
				      * update the fields when the
				      * <tt>reinit</tt> function of the
				      * derived classes is called. The
				      * fields themselves are not set
				      * up, this must happen in the
				      * constructor of the derived
				      * class.
				      */
    FEValuesBase (const unsigned int n_q_points,
		  const unsigned int dofs_per_cell,
		  const UpdateFlags         update_flags,
		  const Mapping<dim>       &mapping,
		  const FiniteElement<dim> &fe);


				     /**
				      * Destructor.
				      */
    ~FEValuesBase ();
				     /// @name ShapeAccess Access to shape function values
				     //@{
    
				     /**
				      * Value of a shape function at a
				      * quadrature point on the cell,
				      * face or subface selected the
				      * last time the <tt>reinit</tt>
				      * function of the derived class
				      * was called.
				      *
				      * If the shape function is
				      * vector-valued, then this
				      * returns the only non-zero
				      * component. If the shape
				      * function has more than one
				      * non-zero component (i.e. it is
				      * not primitive), then throw an
				      * exception of type
				      * ExcShapeFunctionNotPrimitive. In
				      * that case, use the
				      * shape_value_component()
				      * function.
				      *
				      * @arg function_no Number
				      * of the shape function to be
				      * computed
				      * @arg point_no Number of
				      * the quadrature point at which
				      * function is to be computed
				      */
    double shape_value (const unsigned int function_no,
			const unsigned int point_no) const;

				     /**
				      * Compute one vector component of
				      * the value of a shape function
				      * at a quadrature point. If the
				      * finite element is scalar, then
				      * only component zero is allowed
				      * and the return value equals
				      * that of the shape_value()
				      * function. If the finite
				      * element is vector valued but
				      * all shape functions are
				      * primitive (i.e. they are
				      * non-zero in only one
				      * component), then the value
				      * returned by shape_value()
				      * equals that of this function
				      * for exactly one
				      * component. This function is
				      * therefore only of greater
				      * interest if the shape function
				      * is not primitive, but then it
				      * is necessary since the other
				      * function cannot be used.
				      *
				      * @arg function_no Number
				      * of the shape function to be
				      * computed
				      * @arg point_no Number of
				      * the quadrature point at which
				      * function is to be computed
				      * @arg component vector component to be computed
				      */
    double shape_value_component (const unsigned int function_no,
				  const unsigned int point_no,
				  const unsigned int component) const;

    				     /**
				      * Compute the gradient of the
				      * <tt>i</tt>th shape function at the
				      * <tt>j</tt>th quadrature point with
				      * respect to real cell
				      * coordinates.  If you want to
				      * get the derivative in one of
				      * the coordinate directions, use
				      * the appropriate function of
				      * the Tensor class to
				      * extract one component. Since
				      * only a reference to the
				      * gradient's value is returned,
				      * there should be no major
				      * performance drawback.
				      *
				      * If the shape function is
				      * vector-valued, then this
				      * returns the only non-zero
				      * component. If the shape
				      * function has more than one
				      * non-zero component (i.e. it is
				      * not primitive), then throw an
				      * exception of type
				      * ExcShapeFunctionNotPrimitive. In
				      * that case, use the
				      * shape_grad_component()
				      * function.
				      */
    const Tensor<1,dim> &
    shape_grad (const unsigned int function,
		const unsigned int quadrature_point) const;

				     /**
				      * Return one vector component of
				      * the gradient of a shape function
				      * at a quadrature point. If the
				      * finite element is scalar, then
				      * only component zero is allowed
				      * and the return value equals
				      * that of the shape_grad()
				      * function. If the finite
				      * element is vector valued but
				      * all shape functions are
				      * primitive (i.e. they are
				      * non-zero in only one
				      * component), then the value
				      * returned by shape_grad()
				      * equals that of this function
				      * for exactly one
				      * component. This function is
				      * therefore only of greater
				      * interest if the shape function
				      * is not primitive, but then it
				      * is necessary since the other
				      * function cannot be used.
				      */
    Tensor<1,dim>
    shape_grad_component (const unsigned int function_no,
			  const unsigned int point_no,
			  const unsigned int component) const;

    				     /**
				      * Second derivatives of
				      * the <tt>function_no</tt>th shape function at
				      * the <tt>point_no</tt>th quadrature point
				      * with respect to real cell
				      * coordinates. If you want to
				      * get the derivatives in one of
				      * the coordinate directions, use
				      * the appropriate function of
				      * the Tensor class to
				      * extract one component. Since
				      * only a reference to the
				      * derivative values is returned,
				      * there should be no major
				      * performance drawback.
				      *
				      * If the shape function is
				      * vector-valued, then this
				      * returns the only non-zero
				      * component. If the shape
				      * function has more than one
				      * non-zero component (i.e. it is
				      * not primitive), then throw an
				      * exception of type
				      * ExcShapeFunctionNotPrimitive. In
				      * that case, use the
				      * shape_grad_grad_component()
				      * function.
				      */
    const Tensor<2,dim> &
    shape_2nd_derivative (const unsigned int function_no,
			  const unsigned int point_no) const;


				     /**
				      * Return one vector component of
				      * the gradient of a shape
				      * function at a quadrature
				      * point. If the finite element
				      * is scalar, then only component
				      * zero is allowed and the return
				      * value equals that of the
				      * shape_2nd_derivative()
				      * function. If the finite
				      * element is vector valued but
				      * all shape functions are
				      * primitive (i.e. they are
				      * non-zero in only one
				      * component), then the value
				      * returned by
				      * shape_2nd_derivative()
				      * equals that of this function
				      * for exactly one
				      * component. This function is
				      * therefore only of greater
				      * interest if the shape function
				      * is not primitive, but then it
				      * is necessary since the other
				      * function cannot be used.
				      */
    Tensor<2,dim>
    shape_2nd_derivative_component (const unsigned int function_no,
				    const unsigned int point_no,
				    const unsigned int component) const;
    

				     //@}
				     /// @name FunctionAccess Access to values of global finite element functions
				     //@{
    
				     /**
				      * Returns the values of the
				      * finite element function
				      * characterized by
				      * <tt>fe_function</tt> restricted to
				      * the cell, face or subface
				      * selected the last time the
				      * <tt>reinit</tt> function of the
				      * derived class was called, at
				      * the quadrature points.
				      *
				      * If the present cell is not an
				      * active one the interpolated
				      * function values are returned.
				      *
				      * To get values of
				      * multi-component elements,
				      * there is another
				      * get_function_values() below,
				      * returning a vector of vectors
				      * of results.
				      *
				      * This function may only be used if the
				      * finite element in use is a scalar one,
				      * i.e. has only one vector component. If
				      * it is a vector-valued one, then use
				      * the other get_function_values()
				      * function.
				      * 
				      * The function assumes that the
				      * <tt>values</tt> object already has the
				      * correct size. 
				      *
				      * The actual data type of the input
				      * vector may be either a Vector&lt;T&gt;,
				      * BlockVector&lt;T&gt;, or one of the
				      * PETSc vector wrapper classes. It
				      * represents a global vector of
				      * DoF values associated with the
				      * DofHandler object with
				      * which this FEValues
				      * object was last initialized.
				      */
    template <class InputVector, typename number>
    void get_function_values (const InputVector& fe_function,
			      std::vector<number>& values) const;

				     /**
				      * Access to vector valued finite
				      * element functions.
				      *
				      * This function does the same as
				      * the other
				      * get_function_values(), but
				      * applied to multi-component
				      * elements.
				      *
				      * The actual data type of the input
				      * vector may be either a Vector&lt;T&gt;,
				      * BlockVector&lt;T&gt;, or one of the
				      * PETSc vector wrapper classes. It
				      * represents a global vector of
				      * DoF values associated with the
				      * DofHandler object with
				      * which this FEValues
				      * object was last initialized.
				      */
    template <class InputVector, typename number>
    void get_function_values (const InputVector       &fe_function,
			      std::vector<Vector<number> > &values) const;

				     /**
				      * Generate function values from
				      * an arbitrary vector.
				      * 
				      * This function offers the
				      * possibility to extract
				      * function values in quadrature
				      * points from vectors not
				      * corresponding to a whole
				      * discretization.
				      *
				      * You may want to use this
				      * function, if you want to
				      * access just a single block
				      * from a BlockVector, if you
				      * have a multi-level vector or
				      * if you already have a local
				      * representation of your finite
				      * element data.
				      */
    template <class InputVector, typename number>
    void get_function_values (const InputVector& fe_function,
			      const VectorSlice<const std::vector<unsigned int> >& indices,
			      std::vector<number>& values) const;

				     /**
				      * Generate vector function
				      * values from an arbitrary
				      * vector.
				      *
				      * This function offers the
				      * possibility to extract
				      * function values in quadrature
				      * points from vectors not
				      * corresponding to a whole
				      * discretization.
				      *
				      * The length of the vector
				      * <tt>indices</tt> may even be a
				      * multiple of the number of dofs
				      * per cell. Then, the vectors in
				      * <tt>value</tt> should allow
				      * for the same multiple of the
				      * components of the finite
				      * element.
				      *
				      * You may want to use this
				      * function, if you want to
				      * access just a single block
				      * from a BlockVector, if you
				      * have a multi-level vector or
				      * if you already have a local
				      * representation of your finite
				      * element data.
				      *
				      * Since this function allows for
				      * fairly general combinations of
				      * argument sizes, be aware that
				      * the checks on the arguments
				      * may not detect errors.
				      */
    template <class InputVector, typename number>
    void get_function_values (const InputVector& fe_function,
			      const VectorSlice<const std::vector<unsigned int> >& indices,
			      std::vector<Vector<number> >& values) const;


				     /**
				      * Generate vector function
				      * values from an arbitrary
				      * vector.
				      *
				      * This function offers the
				      * possibility to extract
				      * function values in quadrature
				      * points from vectors not
				      * corresponding to a whole
				      * discretization.
				      *
				      * The length of the vector
				      * <tt>indices</tt> may even be a
				      * multiple of the number of dofs
				      * per cell. Then, the vectors in
				      * <tt>value</tt> should allow
				      * for the same multiple of the
				      * components of the finite
				      * element.
				      *
				      * Depending on the last
				      * argument, the outer vector of
				      * <tt>values</tt> has either the
				      * length of the quadrature rule
				      * (<tt>quadrature_points_fastest
				      * == false</tt>) or the length
				      * of components to be filled
				      * <tt>quadrature_points_fastest
				      * == true</tt>. If <tt>p</tt> is
				      * the xurrent quadrature point
				      * number and <tt>i</tt> is the
				      * vector component of the
				      * solution desired, the access
				      * to <tt>values</tt> is
				      * <tt>values[p][i]</tt> if
				      * <tt>quadrature_points_fastest
				      * == false</tt>, and
				      * <tt>values[i][p]</tt>
				      * otherwise.
				      *
				      * You may want to use this
				      * function, if you want to
				      * access just a single block
				      * from a BlockVector, if you
				      * have a multi-level vector or
				      * if you already have a local
				      * representation of your finite
				      * element data.
				      *
				      * Since this function allows for
				      * fairly general combinations of
				      * argument sizes, be aware that
				      * the checks on the arguments
				      * may not detect errors.
				      */
    template <class InputVector, typename number>
    void get_function_values (const InputVector& fe_function,
			      const VectorSlice<const std::vector<unsigned int> >& indices,
			      std::vector<std::vector<number> >& values,
			      bool quadrature_points_fastest) const;

				     /**
				      * Compute the gradients of the finite
				      * element function characterized
				      * by @p fe_function restricted to
				      * @p cell at the quadrature points.
				      *
				      * If the present cell is not an active
				      * one the interpolated function values
				      * are returned.
				      *
				      * This function may only be used if the
				      * finite element in use is a scalar one,
				      * i.e. has only one vector component. If
				      * it is a vector-valued one, then use
				      * the other get_function_grads()
				      * function.
				      * 
				      * The function assumes that the
				      * @p gradients object already has the
				      * right size.
				      *
				      * The actual data type of the input
				      * vector may be either a Vector&lt;T&gt;,
				      * BlockVector&lt;T&gt;, or one of the
				      * PETSc vector wrapper classes. It
				      * represents a global vector of
				      * DoF values associated with the
				      * DofHandler object with
				      * which this FEValues
				      * object was last initialized.
				      *
				      * The output are the gradients
				      * of the function represented by
				      * these DoF values, as computed
				      * in real space (as opposed to
				      * on the unit cell).
				      */
    template <class InputVector>
    void get_function_grads (const InputVector      &fe_function,
			     std::vector<Tensor<1,dim> > &gradients) const;

				     /**
				      * Compute the gradients of the finite
				      * element function characterized
				      * by @p fe_function restricted to
				      * @p cell at the quadrature points.
				      *
				      * If the present cell is not an active
				      * one the interpolated function values
				      * are returned.
				      *
				      * The function assumes that the
				      * @p gradients object already has the
				      * right size.
				      *
				      * This function does the same as
				      * the other get_function_values(),
				      * but applied to multi-component
				      * elements.
				      *
				      * The actual data type of the input
				      * vector may be either a Vector&lt;T&gt;,
				      * BlockVector&lt;T&gt;, or one of the
				      * PETSc vector wrapper classes. It
				      * represents a global vector of
				      * DoF values associated with the
				      * DofHandler object with
				      * which this FEValues
				      * object was last initialized.
				      *
				      * The output are the gradients
				      * of the function represented by
				      * these DoF values, as computed
				      * in real space (as opposed to
				      * on the unit cell).
				      */
    template <class InputVector>
    void get_function_grads (const InputVector               &fe_function,
			     std::vector<std::vector<Tensor<1,dim> > > &gradients) const;

				     /**
				      * Function gradient access with
				      * more flexibility. see
				      * get_function_values() with
				      * corresponding arguments.
				      */
    template <class InputVector>
    void get_function_grads (const InputVector& fe_function,
			     const VectorSlice<const std::vector<unsigned int> >& indices,
			     std::vector<Tensor<1,dim> >& gradients) const;

				     /**
				      * Function gradient access with
				      * more flexibility. see
				      * get_function_values() with
				      * corresponding arguments.
				      */
    template <class InputVector>
    void get_function_grads (const InputVector& fe_function,
			     const VectorSlice<const std::vector<unsigned int> >& indices,
			     std::vector<std::vector<Tensor<1,dim> > >& gradients,
			     bool quadrature_points_fastest = false) const;

				     /**
				      * Compute the tensor of second
				      * derivatives of the finite
				      * element function characterized
				      * by @p fe_function restricted
				      * to @p cell at the quadrature
				      * points.
				      *
				      * The function assumes that the
				      * @p second_derivatives object
				      * already has the correct size.
				      *
				      * This function may only be used if the
				      * finite element in use is a scalar one,
				      * i.e. has only one vector component. If
				      * it is a vector-valued one, then use
				      * the other
				      * get_function_2nd_derivatives()
				      * function.
				      * 
				      * The actual data type of the input
				      * vector may be either a Vector&lt;T&gt;,
				      * BlockVector&lt;T&gt;, or one of the
				      * PETSc vector wrapper classes..It
				      * represents a global vector of
				      * DoF values associated with the
				      * DofHandler object with
				      * which this FEValues
				      * object was last initialized.
				      *
				      * The output are the second derivatives
				      * of the function represented by
				      * these DoF values, as computed
				      * in real space (as opposed to
				      * on the unit cell).
				      */
    template <class InputVector>
    void
    get_function_2nd_derivatives (const InputVector& fe_function,
                                  std::vector<Tensor<2,dim> >& second_derivatives) const;

    
				     /**
				      * Compute the tensor of second
				      * derivatives of the finite
				      * element function characterized
				      * by @p fe_function restricted to
				      * @p cell at the quadrature points.
				      *
				      * The function assumes that the
				      * @p second_derivatives object already has
				      * the right size.
				      *
				      * This function does the same as
				      * the other one with the same
				      * name, but applies to
				      * vector-valued finite elements.
				      *
				      * The actual data type of the input
				      * vector may be either a Vector&lt;T&gt;,
				      * BlockVector&lt;T&gt;, or one of the
				      * PETSc vector wrapper classes. It
				      * represents a global vector of
				      * DoF values associated with the
				      * DofHandler object with
				      * which this FEValues
				      * object was last initialized.
				      *
				      * The output are the second derivatives
				      * of the function represented by
				      * these DoF values, as computed
				      * in real space (as opposed to
				      * on the unit cell).
				      */
    template <class InputVector>
    void
    get_function_2nd_derivatives (const InputVector      &fe_function,
                                  std::vector<std::vector<Tensor<2,dim> > > &second_derivatives,
				  bool quadrature_points_fastest = false) const;
				     //@}
    
				     /**
				      * Position of the <tt>i</tt>th
				      * quadrature point in real space.
				      */
    const Point<dim> & quadrature_point (const unsigned int i) const;

				     /**
				      * Return a pointer to the vector of
				      * quadrature points.
				      */
    const std::vector<Point<dim> > & get_quadrature_points () const;

				     /**
				      * Mapped quadrature weight. If
				      * this object refers to a volume
				      * evaluation (i.e. the derived
				      * class is of type FEValues),
				      * then this is the Jacobi
				      * determinant times the weight
				      * of the *<tt>i</tt>th unit
				      * quadrature point.
				      *
				      * For surface evaluations
				      * (i.e. classes FEFaceValues or
				      * FESubfaceValues), is the
				      * mapped surface element times
				      * the weight of the quadrature
				      * point.
				      */
    double JxW (const unsigned int quadrature_point) const;

				     /**
				      * Pointer to the array holding
				      * the values returned by JxW().
				      */
    const std::vector<double> & get_JxW_values () const;
    
				     /**
				      * Constant reference to the
				      * selected mapping object.
				      */
    const Mapping<dim> & get_mapping () const;

				     /**
				      * Constant reference to the
				      * selected finite element
				      * object.
				      */
    const FiniteElement<dim> & get_fe () const;

				     /**
				      * Return a triangulation
				      * iterator to the current cell.
				      */
    const typename Triangulation<dim>::cell_iterator get_cell () const;

				     /**
				      * Determine an estimate for the
				      * memory consumption (in bytes)
				      * of this object.
				      */
    unsigned int memory_consumption () const;

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
    DeclException0 (ExcInvalidUpdateFlag);
				     /**
				      * Exception
				      */
    DeclException0 (ExcFEDontMatch);
				     /**
				      * Exception
				      */
    DeclException1 (ExcShapeFunctionNotPrimitive,
		    int,
		    << "The shape function with index " << arg1
		    << " is not primitive, i.e. it is vector-valued and "
		    << "has more than one non-zero vector component. This "
		    << "function cannot be called for these shape functions. "
		    << "Maybe you want to use the same function with the "
		    << "_component suffix?");
				     /**
				      * Exception
				      */
    DeclException0 (ExcFENotPrimitive);
    
  protected:
                                     /**
                                      * Objects of the FEValues
                                      * class need to store a pointer
                                      * (i.e. an iterator) to the
                                      * present cell in order to be
                                      * able to extract the values of
                                      * the degrees of freedom on this
                                      * cell in the
                                      * get_function_values() and
                                      * assorted functions. On the
                                      * other hand, this class should
                                      * also work for different
                                      * iterators, as long as they
                                      * have the same interface to
                                      * extract the DoF values (i.e.,
                                      * for example, they need to have
                                      * a @p get_interpolated_dof_values
                                      * function).
                                      *
                                      * This calls for a common base
                                      * class of iterator classes, and
                                      * making the functions we need
                                      * here @p virtual. On the other
                                      * hand, this is the only place
                                      * in the library where we need
                                      * this, and introducing a base
                                      * class of iterators and making
                                      * a function virtual penalizes
                                      * <em>all</em> users of the
                                      * iterators, which are basically
                                      * intended as very fast accessor
                                      * functions. So we do not want
                                      * to do this. Rather, what we do
                                      * here is making the functions
                                      * we need virtual only for use
                                      * with <em>this class</em>. The idea
                                      * is the following: have a
                                      * common base class which
                                      * declares some pure virtual
                                      * functions, and for each
                                      * possible iterator type, we
                                      * have a derived class which
                                      * stores the iterator to the
                                      * cell and implements these
                                      * functions. Since the iterator
                                      * classes have the same
                                      * interface, we can make the
                                      * derived classes a template,
                                      * templatized on the iterator
                                      * type.
                                      *
                                      * This way, the use of virtual
                                      * functions is restricted to
                                      * only this class, and other
                                      * users of iterators do not have
                                      * to bear the negative effects.
                                      *
                                      * @author Wolfgang Bangerth, 2003
                                      */
    class CellIteratorBase 
    {
      public:
                                         /**
                                          * Destructor. Made virtual
                                          * since we store only
                                          * pointers to the base
                                          * class.
                                          */
        virtual ~CellIteratorBase ();
        
                                         /**
                                          * Conversion operator to an
                                          * iterator for
                                          * triangulations. This
                                          * conversion is implicit for
                                          * the original iterators,
                                          * since they are derived
                                          * classes. However, since
                                          * here we have kind of a
                                          * parallel class hierarchy,
                                          * we have to have a
                                          * conversion operator.
                                          */
        virtual
        operator const typename Triangulation<dim>::cell_iterator () const = 0;

                                         /**
                                          * Return the number of
                                          * degrees of freedom the DoF
                                          * handler object has to
                                          * which the iterator belongs
                                          * to.
                                          */
        virtual
        unsigned int
        n_dofs_for_dof_handler () const = 0;

                                         /**
                                          * Call
                                          * @p get_interpolated_dof_values
                                          * of the iterator with the
                                          * given arguments.
                                          */
        virtual
        void
        get_interpolated_dof_values (const Vector<double> &in,
                                     Vector<double>       &out) const = 0;

                                         /**
                                          * Call
                                          * @p get_interpolated_dof_values
                                          * of the iterator with the
                                          * given arguments.
                                          */
        virtual
        void
        get_interpolated_dof_values (const Vector<float> &in,
                                     Vector<float>       &out) const = 0;

                                         /**
                                          * Call
                                          * @p get_interpolated_dof_values
                                          * of the iterator with the
                                          * given arguments.
                                          */
        virtual
        void
        get_interpolated_dof_values (const BlockVector<double> &in,
                                     Vector<double>       &out) const = 0;

                                         /**
                                          * Call
                                          * @p get_interpolated_dof_values
                                          * of the iterator with the
                                          * given arguments.
                                          */
        virtual
        void
        get_interpolated_dof_values (const BlockVector<float> &in,
                                     Vector<float>       &out) const = 0;

#ifdef DEAL_II_USE_PETSC
                                         /**
                                          * Call
                                          * @p get_interpolated_dof_values
                                          * of the iterator with the
                                          * given arguments.
                                          */
        virtual
        void
        get_interpolated_dof_values (const PETScWrappers::Vector &in,
                                     Vector<PetscScalar>         &out) const = 0;

                                         /**
                                          * Call
                                          * @p get_interpolated_dof_values
                                          * of the iterator with the
                                          * given arguments.
                                          */
        virtual
        void
        get_interpolated_dof_values (const PETScWrappers::BlockVector &in,
                                     Vector<PetscScalar>              &out) const = 0;
#endif
    };


                                     /**
                                      * Implementation of derived
                                      * classes of the
                                      * CellIteratorBase
                                      * interface. See there for a
                                      * description of the use of
                                      * these classes.
                                      * 
                                      * @author Wolfgang Bangerth, 2003
                                      */
    template <typename CI>
    class CellIterator : public CellIteratorBase
    {
      public:
                                         /**
                                          * Constructor. Take an
                                          * iterator and store it in
                                          * this class.
                                          */
        CellIterator (const CI &cell);
        
                                         /**
                                          * Conversion operator to an
                                          * iterator for
                                          * triangulations. This
                                          * conversion is implicit for
                                          * the original iterators,
                                          * since they are derived
                                          * classes. However, since
                                          * here we have kind of a
                                          * parallel class hierarchy,
                                          * we have to have a
                                          * conversion operator.
                                          */
        virtual
        operator const typename Triangulation<dim>::cell_iterator () const;
        
                                         /**
                                          * Return the number of
                                          * degrees of freedom the DoF
                                          * handler object has to
                                          * which the iterator belongs
                                          * to.
                                          */
        virtual
        unsigned int
        n_dofs_for_dof_handler () const;

                                         /**
                                          * Call
                                          * @p get_interpolated_dof_values
                                          * of the iterator with the
                                          * given arguments.
                                          */
        virtual
        void
        get_interpolated_dof_values (const Vector<double> &in,
                                     Vector<double>       &out) const;

                                         /**
                                          * Call
                                          * @p get_interpolated_dof_values
                                          * of the iterator with the
                                          * given arguments.
                                          */
        virtual
        void
        get_interpolated_dof_values (const Vector<float> &in,
                                     Vector<float>       &out) const;

                                         /**
                                          * Call
                                          * @p get_interpolated_dof_values
                                          * of the iterator with the
                                          * given arguments.
                                          */
        virtual
        void
        get_interpolated_dof_values (const BlockVector<double> &in,
                                     Vector<double>       &out) const;

                                         /**
                                          * Call
                                          * @p get_interpolated_dof_values
                                          * of the iterator with the
                                          * given arguments.
                                          */
        virtual
        void
        get_interpolated_dof_values (const BlockVector<float> &in,
                                     Vector<float>       &out) const;

#ifdef DEAL_II_USE_PETSC
                                         /**
                                          * Call
                                          * @p get_interpolated_dof_values
                                          * of the iterator with the
                                          * given arguments.
                                          */
        virtual
        void
        get_interpolated_dof_values (const PETScWrappers::Vector &in,
                                     Vector<PetscScalar>         &out) const;

                                         /**
                                          * Call
                                          * @p get_interpolated_dof_values
                                          * of the iterator with the
                                          * given arguments.
                                          */
        virtual
        void
        get_interpolated_dof_values (const PETScWrappers::BlockVector &in,
                                     Vector<PetscScalar>              &out) const;
#endif
      private:
                                         /**
                                          * Copy of the iterator which
                                          * we use in this object.
                                          */
        const CI cell;
    };


                                     /**
                                      * Implementation of a derived
                                      * class of the
                                      * CellIteratorBase
                                      * interface. See there for a
                                      * description of the use of
                                      * these classes.
                                      *
                                      * This class is basically a
                                      * specialization of the general
                                      * template for iterators into
                                      * Triangulation objects (but
                                      * since C++ does not allow
                                      * something like this for nested
                                      * classes, it runs under a
                                      * separate name). Since these do
                                      * not implement the interface
                                      * that we would like to call,
                                      * the functions of this class
                                      * cannot be implemented
                                      * meaningfully. However, most
                                      * functions of the FEValues
                                      * class do not make any use of
                                      * degrees of freedom at all, so
                                      * it should be possible to call
                                      * FEValues::reinit() with a tria
                                      * iterator only; this class
                                      * makes this possible, but
                                      * whenever one of the functions
                                      * of FEValues tries to call
                                      * any of the functions of this
                                      * class, an exception will be
                                      * raised reminding the user that
                                      * if she wants to use these
                                      * features, then the
                                      * FEValues object has to be
                                      * reinitialized with a cell
                                      * iterator that allows to
                                      * extract degree of freedom
                                      * information.
                                      * 
                                      * @author Wolfgang Bangerth, 2003
                                      */
    class TriaCellIterator : public CellIteratorBase
    {
      public:
                                         /**
                                          * Constructor. Take an
                                          * iterator and store it in
                                          * this class.
                                          */
        TriaCellIterator (const typename Triangulation<dim>::cell_iterator &cell);
        
                                         /**
                                          * Conversion operator to an
                                          * iterator for
                                          * triangulations. This
                                          * conversion is implicit for
                                          * the original iterators,
                                          * since they are derived
                                          * classes. However, since
                                          * here we have kind of a
                                          * parallel class hierarchy,
                                          * we have to have a
                                          * conversion operator. Here,
                                          * the conversion is trivial,
                                          * from and to the same time.
                                          */
        virtual
        operator const typename Triangulation<dim>::cell_iterator () const;

                                         /**
                                          * Implement the respective
                                          * function of the base
                                          * class. Since this is not
                                          * possible, we just raise an
                                          * error.
                                          */
        virtual
        unsigned int
        n_dofs_for_dof_handler () const;

                                         /**
                                          * Implement the respective
                                          * function of the base
                                          * class. Since this is not
                                          * possible, we just raise an
                                          * error.
                                          */
        virtual
        void
        get_interpolated_dof_values (const Vector<double> &in,
                                     Vector<double>       &out) const;

                                         /**
                                          * Implement the respective
                                          * function of the base
                                          * class. Since this is not
                                          * possible, we just raise an
                                          * error.
                                          */
        virtual
        void
        get_interpolated_dof_values (const Vector<float> &in,
                                     Vector<float>       &out) const;

                                         /**
                                          * Implement the respective
                                          * function of the base
                                          * class. Since this is not
                                          * possible, we just raise an
                                          * error.
                                          */
        virtual
        void
        get_interpolated_dof_values (const BlockVector<double> &in,
                                     Vector<double>       &out) const;

                                         /**
                                          * Implement the respective
                                          * function of the base
                                          * class. Since this is not
                                          * possible, we just raise an
                                          * error.
                                          */
        virtual
        void
        get_interpolated_dof_values (const BlockVector<float> &in,
                                     Vector<float>       &out) const;

#ifdef DEAL_II_USE_PETSC
                                         /**
                                          * Call
                                          * @p get_interpolated_dof_values
                                          * of the iterator with the
                                          * given arguments.
                                          */
        virtual
        void
        get_interpolated_dof_values (const PETScWrappers::Vector &in,
                                     Vector<PetscScalar>         &out) const;

                                         /**
                                          * Call
                                          * @p get_interpolated_dof_values
                                          * of the iterator with the
                                          * given arguments.
                                          */
        virtual
        void
        get_interpolated_dof_values (const PETScWrappers::BlockVector &in,
                                     Vector<PetscScalar>              &out) const;
#endif
      private:
                                         /**
                                          * Copy of the iterator which
                                          * we use in this object.
                                          */
        const typename Triangulation<dim>::cell_iterator cell;

                                         /**
                                          * String to be displayed
                                          * whenever one of the
                                          * functions of this class is
                                          * called. Make it a static
                                          * member variable, since we
                                          * show the same message for
                                          * all member functions.
                                          */
        static const char * const message_string;
    };
    
				     /**
				      * Store the cell selected last time
				      * the reinit() function was called
				      * to make access
				      * to the <tt>get_function_*</tt> functions
				      * safer.
				      */
    std::auto_ptr<const CellIteratorBase> present_cell;
    
				     /**
				      * Storage for the mapping object.
				      */
    const SmartPointer<const Mapping<dim> > mapping;
    
				     /**
				      * Store the finite element for later use.
				      */
    const SmartPointer<const FiniteElement<dim> > fe;

    
				     /**
				      * Internal data of mapping.
				      */
    SmartPointer<typename Mapping<dim>::InternalDataBase> mapping_data;

				     /**
				      * Internal data of finite element.
				      */
    SmartPointer<typename Mapping<dim>::InternalDataBase> fe_data;

				     /**
				      * Initialize some update
				      * flags. Called from the
				      * @p initialize functions of
				      * derived classes, which are in
				      * turn called from their
				      * constructors.
				      *
				      * Basically, this function finds
				      * out using the finite element
				      * and mapping object already
				      * stored which flags need to be
				      * set to compute everything the
				      * user wants, as expressed
				      * through the flags passed as
				      * argument.
				      */
    UpdateFlags compute_update_flags (const UpdateFlags update_flags) const;

				     /**
				      * Returns reference to default
				      * MappingQ1 object. Needed
				      * by constructors of derived
				      * classes that uses
				      * MappingQ1 implicitly.
				      */
    static const Mapping<dim> &get_default_mapping();

  private:
                                     /**
                                      * Copy constructor. Since
                                      * objects of this class are not
                                      * copyable, we make it private,
                                      * and also do not implement it.
                                      */
    FEValuesBase (const FEValuesBase &);
    
                                     /**
                                      * Copy operator. Since
                                      * objects of this class are not
                                      * copyable, we make it private,
                                      * and also do not implement it.
                                      */
    FEValuesBase & operator= (const FEValuesBase &);    
};



/**
 * Finite element evaluated in quadrature points of a cell.
 *
 * This function implements the initialization routines for
 * FEValuesBase, if values in quadrature points of a cell are
 * needed. For further documentation see this class.
 * 
 * @author Wolfgang Bangerth, 1998, Guido Kanschat, 2001
 */
template <int dim>
class FEValues : public FEValuesBase<dim>
{
  public:
				     /**
				      * Constructor. Gets cell
				      * independent data from mapping
				      * and finite element objects,
				      * matching the quadrature rule
				      * and update flags.
				      */
    FEValues (const Mapping<dim>       &mapping,
	      const FiniteElement<dim> &fe,
	      const Quadrature<dim>    &quadrature,
	      const UpdateFlags         update_flags);

                                     /**
				      * Constructor. Uses MappingQ1
				      * implicitly.
				      */
    FEValues (const FiniteElement<dim> &fe,
	      const Quadrature<dim>    &quadrature,
	      const UpdateFlags         update_flags);
    
				     /**
				      * Reinitialize the gradients,
				      * Jacobi determinants, etc for
				      * the given cell of type
				      * "iterator into a DoFHandler
				      * object", and the finite
				      * element associated with this
				      * object. It is assumed that the
				      * finite element used by the
				      * given cell is also the one
				      * used by this FEValues
				      * object.
				      */
    void reinit (const typename DoFHandler<dim>::cell_iterator &cell);

				     /**
				      * Reinitialize the gradients,
				      * Jacobi determinants, etc for
				      * the given cell of type
				      * "iterator into a MGDoFHandler
				      * object", and the finite
				      * element associated with this
				      * object. It is assumed that the
				      * finite element used by the
				      * given cell is also the one
				      * used by this FEValues
				      * object.
				      */
    void reinit (const typename MGDoFHandler<dim>::cell_iterator &cell);

				     /**
				      * Reinitialize the gradients,
				      * Jacobi determinants, etc for
				      * the given cell of type
				      * "iterator into a Triangulation
				      * object", and the given finite
				      * element. Since iterators into
				      * triangulation alone only
				      * convey information about the
				      * geometry of a cell, but not
				      * about degrees of freedom
				      * possibly associated with this
				      * cell, you will not be able to
				      * call some functions of this
				      * class if they need information
				      * about degrees of
				      * freedom. These functions are,
				      * above all, the
				      * <tt>get_function_value/grads/2nd_derivatives</tt>
				      * functions. If you want to call
				      * these functions, you have to
				      * call the @p reinit variants
				      * that take iterators into
				      * DoFHandler or other DoF
				      * handler type objects.
				      */
    void reinit (const typename Triangulation<dim>::cell_iterator &cell);

				     /**
				      * Return a reference to the copy
				      * of the quadrature formula
				      * stored by this object.
				      */
    const Quadrature<dim> & get_quadrature () const;
    
				     /**
				      * Determine an estimate for the
				      * memory consumption (in bytes)
				      * of this object.
				      */
    unsigned int memory_consumption () const;

                                     /**
                                      * Return a reference to this
                                      * very object.
                                      *
                                      * Though it seems that it is not very
                                      * useful, this function is there to
                                      * provide capability to the hpFEValues
                                      * class, in which case it provides the
                                      * FEValues object for the present cell
                                      * (remember that for hp finite elements,
                                      * the actual FE object used may change
                                      * from cell to cell, so we also need
                                      * different FEValues objects for
                                      * different cells; once you reinitialize
                                      * the hpFEValues object for a specific
                                      * cell, it retrieves the FEValues object
                                      * for the FE on that cell and returns it
                                      * through a function of the same name as
                                      * this one; this function here therefore
                                      * only provides the same interface so
                                      * that one can templatize on
                                      * FEValues/hpFEValues).
                                      */
    const FEValues<dim> & get_present_fe_values () const;
    
  private:
				     /**
				      * Store a copy of the quadrature
				      * formula here.
				      */
    const Quadrature<dim> quadrature;

				     /**
				      * Do work common to the two
				      * constructors.
				      */
    void initialize (const UpdateFlags update_flags);

                                     /**
                                      * The reinit() functions do
                                      * only that part of the work
                                      * that requires knowledge of the
                                      * type of iterator. After
                                      * setting present_cell(),
                                      * they pass on to this function,
                                      * which does the real work, and
                                      * which is independent of the
                                      * actual type of the cell
                                      * iterator.
                                      */
    void do_reinit ();
};


/**
 * Extend the interface of FEValuesBase to values that only make sense
 * when evaluating something on the surface of a cell. All the data
 * that is available in the interior of cells is also available here.
 *
 * On surfaces of mesh cells, normal vectors and boundary forms are
 * additional values that can be computed. This class provides the
 * interface to access those. Implementations are in derived classes
 * FEFaceValues and FESubfaceValues.
 *
 * The boundary form is the cross product of the images of the unit
 * tangential vectors. Therefore, it is the unit normal vector
 * multiplied with the surface element. Since it may be cheaper to
 * compute the boundary form immediately, use this value to integrate
 * <tt>n.ds</tt>.
 *
 * See FEValuesBase
 *
 *  @author Wolfgang Bangerth, 1998, Guido Kanschat, 2000, 2001
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
				      * @p n_faces_or_subfaces is the number
				      * of faces or subfaces that this object
				      * is to store. The actual number depends
				      * on the derived class, for
				      * FEFaceValues it is <tt>2*dim</tt>, while for
				      * the FESubfaceValues class it is
				      * <tt>2*dim*(1<<(dim-1))</tt>, i.e. the number
				      * of faces times the number of subfaces
				      * per face.
				      */
    FEFaceValuesBase (const unsigned int n_q_points,
		      const unsigned int dofs_per_cell,
		      const UpdateFlags         update_flags,
		      const Mapping<dim>       &mapping,
		      const FiniteElement<dim> &fe,
		      const Quadrature<dim-1>& quadrature);

    				     /**
				      * Return the outward normal vector to
				      * the cell at the <tt>i</tt>th quadrature
				      * point. The length of the vector
				      * is normalized to one.
				      */
    const Point<dim> & normal_vector (const unsigned int i) const;
    
    				     /**
				      * Boundary form of the
				      * transformation of the cell at
				      * the <tt>i</tt>th quadrature point.
				      *
				      * The boundary form is the cross
				      * product of the images of the
				      * unit tangential
				      * vectors. Therefore, it is the
				      * unit normal vector multiplied
				      * with the surface
				      * element. Since it may be
				      * cheaper to compute the
				      * boundary form immediately, use
				      * this value to integrate
				      * <tt>n.ds</tt>.
				      */
    const Tensor<1,dim> & boundary_form (const unsigned int i) const;
    
				     /**
				      * Return the list of outward normal
				      * vectors to the cell at the
				      * quadrature points.
				      */
    const std::vector<Point<dim> > & get_normal_vectors () const;

				     /**
				      * Return the list of outward
				      * normal vectors times the
				      * Jacobian of the surface
				      * mapping.
				      */
    const std::vector<Tensor<1,dim> > & get_boundary_forms () const;

				     /**
				      * Return a reference to the copy
				      * of the quadrature formula
				      * stored by this object.
				      */
    const Quadrature<dim-1> & get_quadrature () const;
    
				     /**
				      * Determine an estimate for the
				      * memory consumption (in bytes)
				      * of this object.
				      */
    unsigned int memory_consumption () const;
    
  protected:
				     /**
				      * Store a copy of the quadrature
				      * formula here.
				      */
    const Quadrature<dim-1> quadrature;
};



/**
 * Finite element evaluated in quadrature points on a face.
 *
 * This class adds the functionality of FEFaceValuesBase to
 * FEValues; see there for more documentation.
 *
 * Since finite element functions and their derivatives may be
 * discontinuous at cell boundaries, there is no restriction of this
 * function to a mesh face. But, there are limits of these values
 * approaching the face from either of the neighboring cells.
 *
 * @author Wolfgang Bangerth, 1998, Guido Kanschat, 2000, 2001
 */
template <int dim>
class FEFaceValues : public FEFaceValuesBase<dim>
{
  public:
				     /**
				      * Constructor. Gets cell
				      * independent data from mapping
				      * and finite element objects,
				      * matching the quadrature rule
				      * and update flags.
				      */
    FEFaceValues (const Mapping<dim>       &mapping,
		  const FiniteElement<dim> &fe,
		  const Quadrature<dim-1>  &quadrature,
		  const UpdateFlags         update_flags);

                                     /**
				      * Constructor. Uses MappingQ1
				      * implicitly.
				      */
    FEFaceValues (const FiniteElement<dim> &fe,
		  const Quadrature<dim-1>  &quadrature,
		  const UpdateFlags         update_flags);

				     /**
				      * Reinitialize the gradients, Jacobi
				      * determinants, etc for the face with
				      * number @p face_no of @p cell
				      * and the given finite element.
				      */
    void reinit (const typename DoFHandler<dim>::cell_iterator &cell,
		 const unsigned int                            face_no);
    
				     /**
				      * Reinitialize the gradients,
				      * Jacobi determinants, etc for
				      * the given cell of type
				      * "iterator into a MGDoFHandler
				      * object", and the finite
				      * element associated with this
				      * object. It is assumed that the
				      * finite element used by the
				      * given cell is also the one
				      * used by this FEValues
				      * object.
				      */
    void reinit (const typename MGDoFHandler<dim>::cell_iterator &cell,
		 const unsigned int                               face_no);

				     /**
				      * Reinitialize the gradients,
				      * Jacobi determinants, etc for
				      * the given face on given cell
				      * of type "iterator into a
				      * Triangulation object", and the
				      * given finite element. Since
				      * iterators into triangulation
				      * alone only convey information
				      * about the geometry of a cell,
				      * but not about degrees of
				      * freedom possibly associated
				      * with this cell, you will not
				      * be able to call some functions
				      * of this class if they need
				      * information about degrees of
				      * freedom. These functions are,
				      * above all, the
				      * <tt>get_function_value/grads/2nd_derivatives</tt>
				      * functions. If you want to call
				      * these functions, you have to
				      * call the @p reinit variants
				      * that take iterators into
				      * DoFHandler or other DoF
				      * handler type objects.
				      */
    void reinit (const typename Triangulation<dim>::cell_iterator &cell,
		 const unsigned int                    face_no);
    
                                     /**
                                      * Return a reference to this
                                      * very object.
                                      *
                                      * Though it seems that it is not very
                                      * useful, this function is there to
                                      * provide capability to the hpFEValues
                                      * class, in which case it provides the
                                      * FEValues object for the present cell
                                      * (remember that for hp finite elements,
                                      * the actual FE object used may change
                                      * from cell to cell, so we also need
                                      * different FEValues objects for
                                      * different cells; once you reinitialize
                                      * the hpFEValues object for a specific
                                      * cell, it retrieves the FEValues object
                                      * for the FE on that cell and returns it
                                      * through a function of the same name as
                                      * this one; this function here therefore
                                      * only provides the same interface so
                                      * that one can templatize on
                                      * FEValues/hpFEValues).
                                      */
    const FEFaceValues<dim> & get_present_fe_values () const;
  private:

				     /**
				      * Do work common to the two
				      * constructors.
				      */
    void initialize (const UpdateFlags update_flags);    

                                     /**
                                      * The reinit() functions do
                                      * only that part of the work
                                      * that requires knowledge of the
                                      * type of iterator. After
                                      * setting present_cell(),
                                      * they pass on to this function,
                                      * which does the real work, and
                                      * which is independent of the
                                      * actual type of the cell
                                      * iterator.
                                      */
    void do_reinit (const unsigned int face_no);
};


/**
 * Finite element evaluated in quadrature points on a face.
 *
 * This class adds the functionality of FEFaceValuesBase to
 * FEValues; see there for more documentation.
 *
 * This class is used for faces lying on a refinement edge. In this
 * case, the neighboring cell is refined. To be able to compute
 * differences between interior and exterior function values, the
 * refinement of the neighboring cell must be simulated on this
 * cell. This is achieved by applying a quadrature rule that simulates
 * the refinement. The resulting data fields are split up to reflect
 * the refinement structure of the neighbor: a subface number
 * corresponds to the number of the child of the neighboring face.
 *
 * @author Wolfgang Bangerth, 1998, Guido Kanschat, 2000, 2001
 */
template <int dim>
class FESubfaceValues : public FEFaceValuesBase<dim>
{
  public:
				     /**
				      * Constructor. Gets cell
				      * independent data from mapping
				      * and finite element objects,
				      * matching the quadrature rule
				      * and update flags.
				      */
    FESubfaceValues (const Mapping<dim>       &mapping,
		     const FiniteElement<dim> &fe,
		     const Quadrature<dim-1>  &face_quadrature,
		     const UpdateFlags         update_flags);

                                     /**
				      * Constructor. Uses MappingQ1
				      * implicitly.
				      */
    FESubfaceValues (const FiniteElement<dim> &fe,
		     const Quadrature<dim-1>  &face_quadrature,
		     const UpdateFlags         update_flags);

				     /**
				      * Reinitialize the gradients,
				      * Jacobi determinants, etc for
				      * the given cell of type
				      * "iterator into a DoFHandler
				      * object", and the finite
				      * element associated with this
				      * object. It is assumed that the
				      * finite element used by the
				      * given cell is also the one
				      * used by this
				      * FESubfaceValues object.
				      */
    void reinit (const typename DoFHandler<dim>::cell_iterator &cell,
		 const unsigned int                    face_no,
		 const unsigned int                    subface_no);

    				     /**
				      * Reinitialize the gradients,
				      * Jacobi determinants, etc for
				      * the given cell of type
				      * "iterator into a MGDoFHandler
				      * object", and the finite
				      * element associated with this
				      * object. It is assumed that the
				      * finite element used by the
				      * given cell is also the one
				      * used by this FEValues
				      * object.
				      */
    void reinit (const typename MGDoFHandler<dim>::cell_iterator &cell,
		 const unsigned int                    face_no,
		 const unsigned int                    subface_no);

				     /**
				      * Reinitialize the gradients,
				      * Jacobi determinants, etc for
				      * the given subface on given
				      * cell of type "iterator into a
				      * Triangulation object", and the
				      * given finite element. Since
				      * iterators into triangulation
				      * alone only convey information
				      * about the geometry of a cell,
				      * but not about degrees of
				      * freedom possibly associated
				      * with this cell, you will not
				      * be able to call some functions
				      * of this class if they need
				      * information about degrees of
				      * freedom. These functions are,
				      * above all, the
				      * <tt>get_function_value/grads/2nd_derivatives</tt>
				      * functions. If you want to call
				      * these functions, you have to
				      * call the @p reinit variants
				      * that take iterators into
				      * DoFHandler or other DoF
				      * handler type objects.
				      */
    void reinit (const typename Triangulation<dim>::cell_iterator &cell,
		 const unsigned int                    face_no,
		 const unsigned int                    subface_no);
    
                                     /**
                                      * Return a reference to this
                                      * very object.
                                      *
                                      * Though it seems that it is not very
                                      * useful, this function is there to
                                      * provide capability to the hpFEValues
                                      * class, in which case it provides the
                                      * FEValues object for the present cell
                                      * (remember that for hp finite elements,
                                      * the actual FE object used may change
                                      * from cell to cell, so we also need
                                      * different FEValues objects for
                                      * different cells; once you reinitialize
                                      * the hpFEValues object for a specific
                                      * cell, it retrieves the FEValues object
                                      * for the FE on that cell and returns it
                                      * through a function of the same name as
                                      * this one; this function here therefore
                                      * only provides the same interface so
                                      * that one can templatize on
                                      * FEValues/hpFEValues).
                                      */
    const FESubfaceValues<dim> & get_present_fe_values () const;

                                     /**
				      * Exception
				      */
    DeclException0 (ExcReinitCalledWithBoundaryFace);

				     /**
				      * Exception
				      */
    DeclException0 (ExcFaceHasNoSubfaces);

  private:

				     /**
				      * Do work common to the two
				      * constructors.
				      */
    void initialize (const UpdateFlags update_flags);

                                     /**
                                      * The reinit() functions do
                                      * only that part of the work
                                      * that requires knowledge of the
                                      * type of iterator. After
                                      * setting present_cell(),
                                      * they pass on to this function,
                                      * which does the real work, and
                                      * which is independent of the
                                      * actual type of the cell
                                      * iterator.
                                      */
    void do_reinit (const unsigned int face_no,
                    const unsigned int subface_no);
};

/*@}*/

///@if NoDoc
/*------------------------ Inline functions: FEValuesBase ------------------------*/


template <int dim>
inline
double
FEValuesBase<dim>::shape_value (const unsigned int i,
				const unsigned int j) const
{
  Assert (this->update_flags & update_values,
	  ExcAccessToUninitializedField());
  Assert (fe->is_primitive (i),
	  ExcShapeFunctionNotPrimitive(i));

				   // if the entire FE is primitive,
				   // then we can take a short-cut:
  if (fe->is_primitive())
    return this->shape_values(i,j);
  else
				     // otherwise, use the mapping
				     // between shape function numbers
				     // and rows. note that by the
				     // assertions above, we know that
				     // this particular shape function
				     // is primitive, so there is no
				     // question to which vector
				     // component the call of this
				     // function refers
    return this->shape_values(this->shape_function_to_row_table[i], j);
}



template <int dim>
inline
double
FEValuesBase<dim>::shape_value_component (const unsigned int i,
					  const unsigned int j,
					  const unsigned int component) const
{
  Assert (this->update_flags & update_values,
	  ExcAccessToUninitializedField());
  Assert (component < fe->n_components(),
	  ExcIndexRange(component, 0, fe->n_components()));
			
				   // if this particular shape
				   // function is primitive, then we
				   // can take a short-cut by checking
				   // whether the requested component
				   // is the only non-zero one (note
				   // that calling
				   // system_to_component_table only
				   // works if the shape function is
				   // primitive):
  if (fe->is_primitive(i))
    {
      if (component == fe->system_to_component_index(i).first)
	return this->shape_values(this->shape_function_to_row_table[i],j);
      else
	return 0;
    }
  else
    {
				       // no, this shape function is
				       // not primitive. then we have
				       // to loop over its components
				       // to find the corresponding
				       // row in the arrays of this
				       // object. before that check
				       // whether the shape function
				       // is non-zero at all within
				       // this component:
      if (fe->get_nonzero_components(i)[component] == false)
	return 0.;

				       // count how many non-zero
				       // component the shape function
				       // has before the one we are
				       // looking for, and add this to
				       // the offset of the first
				       // non-zero component of this
				       // shape function in the arrays
				       // we index presently:
      const unsigned int
	row = (this->shape_function_to_row_table[i]
	       +
	       std::count (fe->get_nonzero_components(i).begin(),
			   fe->get_nonzero_components(i).begin()+component,
			   true));
      return this->shape_values(row, j);
    };
}



template <int dim>
inline
const Tensor<1,dim> &
FEValuesBase<dim>::shape_grad (const unsigned int i,
			       const unsigned int j) const
{
  Assert (this->update_flags & update_gradients,
	  ExcAccessToUninitializedField());
  Assert (fe->is_primitive (i),
	  ExcShapeFunctionNotPrimitive(i));
  Assert (i<this->shape_gradients.size(),
	  ExcIndexRange (i, 0, this->shape_gradients.size()));
  Assert (j<this->shape_gradients[0].size(),
	  ExcIndexRange (j, 0, this->shape_gradients[0].size()));

				   // if the entire FE is primitive,
				   // then we can take a short-cut:
  if (fe->is_primitive())
    return this->shape_gradients[i][j];
  else
				     // otherwise, use the mapping
				     // between shape function numbers
				     // and rows. note that by the
				     // assertions above, we know that
				     // this particular shape function
				     // is primitive, so there is no
				     // question to which vector
				     // component the call of this
				     // function refers
    return this->shape_gradients[this->shape_function_to_row_table[i]][j];
}



template <int dim>
inline
Tensor<1,dim>
FEValuesBase<dim>::shape_grad_component (const unsigned int i,
					 const unsigned int j,
					 const unsigned int component) const
{
  Assert (this->update_flags & update_gradients,
	  ExcAccessToUninitializedField());
  Assert (component < fe->n_components(),
	  ExcIndexRange(component, 0, fe->n_components()));
			
				   // if this particular shape
				   // function is primitive, then we
				   // can take a short-cut by checking
				   // whether the requested component
				   // is the only non-zero one (note
				   // that calling
				   // system_to_component_table only
				   // works if the shape function is
				   // primitive):
  if (fe->is_primitive(i))
    {
      if (component == fe->system_to_component_index(i).first)
	return this->shape_gradients[this->shape_function_to_row_table[i]][j];
      else
	return Tensor<1,dim>();
    }
  else
    {
				       // no, this shape function is
				       // not primitive. then we have
				       // to loop over its components
				       // to find the corresponding
				       // row in the arrays of this
				       // object. before that check
				       // whether the shape function
				       // is non-zero at all within
				       // this component:
      if (fe->get_nonzero_components(i)[component] == false)
	return Tensor<1,dim>();

				       // count how many non-zero
				       // component the shape function
				       // has before the one we are
				       // looking for, and add this to
				       // the offset of the first
				       // non-zero component of this
				       // shape function in the arrays
				       // we index presently:
      const unsigned int
	row = (this->shape_function_to_row_table[i]
	       +
	       std::count (fe->get_nonzero_components(i).begin(),
			   fe->get_nonzero_components(i).begin()+component,
			   true));
      return this->shape_gradients[row][j];
    };
}



template <int dim>
inline
const Tensor<2,dim> &
FEValuesBase<dim>::shape_2nd_derivative (const unsigned int i,
					 const unsigned int j) const
{
  Assert (this->update_flags & update_second_derivatives,
	  ExcAccessToUninitializedField());
  Assert (fe->is_primitive (i),
	  ExcShapeFunctionNotPrimitive(i));
  Assert (i<this->shape_2nd_derivatives.size(),
	  ExcIndexRange (i, 0, this->shape_2nd_derivatives.size()));
  Assert (j<this->shape_2nd_derivatives[0].size(),
	  ExcIndexRange (j, 0, this->shape_2nd_derivatives[0].size()));

				   // if the entire FE is primitive,
				   // then we can take a short-cut:
  if (fe->is_primitive())
    return this->shape_2nd_derivatives[i][j];
  else
				     // otherwise, use the mapping
				     // between shape function numbers
				     // and rows. note that by the
				     // assertions above, we know that
				     // this particular shape function
				     // is primitive, so there is no
				     // question to which vector
				     // component the call of this
				     // function refers
    return this->shape_2nd_derivatives[this->shape_function_to_row_table[i]][j];
}



template <int dim>
inline
Tensor<2,dim>
FEValuesBase<dim>::shape_2nd_derivative_component (const unsigned int i,
						   const unsigned int j,
						   const unsigned int component) const
{
  Assert (this->update_flags & update_second_derivatives,
	  ExcAccessToUninitializedField());
  Assert (component < fe->n_components(),
	  ExcIndexRange(component, 0, fe->n_components()));
			
				   // if this particular shape
				   // function is primitive, then we
				   // can take a short-cut by checking
				   // whether the requested component
				   // is the only non-zero one (note
				   // that calling
				   // system_to_component_table only
				   // works if the shape function is
				   // primitive):
  if (fe->is_primitive(i))
    {
      if (component == fe->system_to_component_index(i).first)
	return this->shape_2nd_derivatives[this->shape_function_to_row_table[i]][j];
      else
	return Tensor<2,dim>();
    }
  else
    {
				       // no, this shape function is
				       // not primitive. then we have
				       // to loop over its components
				       // to find the corresponding
				       // row in the arrays of this
				       // object. before that check
				       // whether the shape function
				       // is non-zero at all within
				       // this component:
      if (fe->get_nonzero_components(i)[component] == false)
	return Tensor<2,dim>();

				       // count how many non-zero
				       // component the shape function
				       // has before the one we are
				       // looking for, and add this to
				       // the offset of the first
				       // non-zero component of this
				       // shape function in the arrays
				       // we index presently:
      const unsigned int
	row = (this->shape_function_to_row_table[i]
	       +
	       std::count (fe->get_nonzero_components(i).begin(),
			   fe->get_nonzero_components(i).begin()+component,
			   true));
      return this->shape_2nd_derivatives[row][j];
    };
}



template <int dim>
inline
const FiniteElement<dim> & 
FEValuesBase<dim>::get_fe () const
{
  return *fe;
}


template <int dim>
inline
const Mapping<dim> & 
FEValuesBase<dim>::get_mapping () const
{
  return *mapping;
}



template <int dim>
inline
const typename Triangulation<dim>::cell_iterator
FEValuesBase<dim>::get_cell () const
{
  return *present_cell;
}


template <int dim>
inline
const std::vector<Point<dim> > &
FEValuesBase<dim>::get_quadrature_points () const
{
  Assert (this->update_flags & update_q_points, ExcAccessToUninitializedField());
  return this->quadrature_points;
}



template <int dim>
inline
const std::vector<double> &
FEValuesBase<dim>::get_JxW_values () const
{
  Assert (this->update_flags & update_JxW_values, ExcAccessToUninitializedField());
  return this->JxW_values;
}



template <int dim>
inline
const Point<dim> &
FEValuesBase<dim>::quadrature_point (const unsigned int i) const
{
  Assert (this->update_flags & update_q_points, ExcAccessToUninitializedField());
  Assert (i<this->quadrature_points.size(), ExcIndexRange(i, 0, this->quadrature_points.size()));
  
  return this->quadrature_points[i];
}




template <int dim>
inline
double
FEValuesBase<dim>::JxW (const unsigned int i) const
{
  Assert (this->update_flags & update_JxW_values, ExcAccessToUninitializedField());
  Assert (i<this->JxW_values.size(), ExcIndexRange(i, 0, this->JxW_values.size()));
  
  return this->JxW_values[i];
}


/*------------------------ Inline functions: FEValues ----------------------------*/


template <int dim>
inline
const Quadrature<dim> &
FEValues<dim>::get_quadrature () const 
{
  return quadrature;
}



template <int dim>
inline
const FEValues<dim> &
FEValues<dim>::get_present_fe_values () const 
{
  return *this;
}


/*------------------------ Inline functions: FEFaceValuesBase --------------------*/


template <int dim>
inline
const Point<dim> &
FEFaceValuesBase<dim>::normal_vector (const unsigned int i) const
{
  Assert (i<this->normal_vectors.size(),
	  ExcIndexRange(i, 0, this->normal_vectors.size()));
  Assert (this->update_flags & update_normal_vectors,
	  typename FEValuesBase<dim>::ExcAccessToUninitializedField());
  
  return this->normal_vectors[i];
}



/*------------------------ Inline functions: FE*FaceValues --------------------*/

template <int dim>
inline
const Quadrature<dim-1> &
FEFaceValuesBase<dim>::get_quadrature () const 
{
  return quadrature;
}



template <int dim>
inline
const FEFaceValues<dim> &
FEFaceValues<dim>::get_present_fe_values () const 
{
  return *this;
}



template <int dim>
inline
const FESubfaceValues<dim> &
FESubfaceValues<dim>::get_present_fe_values () const 
{
  return *this;
}



template <int dim>
inline
const Tensor<1,dim> &
FEFaceValuesBase<dim>::boundary_form (const unsigned int i) const
{
  Assert (i<this->boundary_forms.size(),
	  ExcIndexRange(i, 0, this->boundary_forms.size()));
  Assert (this->update_flags & update_boundary_forms,
	  typename FEValuesBase<dim>::ExcAccessToUninitializedField());
  
  return this->boundary_forms[i];
}

///@endif


#endif
