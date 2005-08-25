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
#ifndef __deal2__fe_h
#define __deal2__fe_h

#include <base/config.h>
#include <fe/fe_base.h>
#include <dofs/dof_handler.h>

template <int dim> class FEValuesData;
template <int dim> class FEValuesBase;
template <int dim> class FEValues;
template <int dim> class FEFaceValues;
template <int dim> class FESubfaceValues;
template <int dim> class FESystem;
template <int dim> class FECollection;

/*!@addtogroup febase */
/*@{*/

/**
 * Common interface of all finite elements. Here, the functions to
 * fill the data fields of FEValues are declared.
 *
 * The interface of this class is very restrictive. The reason is that
 * finite element values should be accessed only by use of
 * FEValues objects. These, together with @p FiniteElement are
 * responsible to provide an optimized implementation.
 *
 * This even holds for evaluating finite elements at their support
 * points (provided the element is based on Lagrangian interpolation):
 * first, it is necessary to construct a quadrature rule from the
 * support points. This is then fed into an object of class
 * FEValues. Even for evaluation on the unit cell, you will need
 * a triangulation containing that single cell.
 *
 * Basically, this class just declares the shape function and their
 * derivatives on the unit cell $[0,1]^d$, and the means to transform
 * them onto a given cell in physical space if provided by the
 * FEValues class with a Mapping object.
 *
 * @author Wolfgang Bangerth, Guido Kanschat, Ralf Hartmann, 1998, 2000, 2001
 */
template <int dim>
class FiniteElement : public Subscriptor,
		      public FiniteElementData<dim>
{
  public:
				   /**
				    * Base class for internal data.
				    * Adds data for second derivatives to
				    * Mapping::InternalDataBase()
				    *
				    * For information about the
				    * general purpose of this class,
				    * see the documentation of the
				    * base class.
				    *
				    * @author Guido Kanschat, 2001
				    */
  class InternalDataBase : public Mapping<dim>::InternalDataBase
    {
      public:      
					 /**
					  * Destructor. Needed to
					  * avoid memory leaks with
					  * difference quotients.
					  */
	virtual ~InternalDataBase ();

					 /**
					  * Initialize some pointers
					  * used in the computation of
					  * second derivatives by
					  * finite differencing of
					  * gradients.
					  */
	void initialize_2nd (const FiniteElement<dim> *element,
			     const Mapping<dim>       &mapping,
			     const Quadrature<dim>    &quadrature);
	
					 /**
					  * Storage for @p FEValues
					  * objects needed to
					  * approximate second
					  * derivatives.
					  *
					  * The ordering is <tt>p+hx</tt>,
					  * <tt>p+hy</tt>, <tt>p+hz</tt>,
					  * @p p-hx, @p p-hy,
					  * @p p-hz, where unused
					  * entries in lower dimensions
					  * are missing.
					  */
	std::vector<FEValues<dim>*> differences;
    };

  public:
                                     /**
                                      * Copy constructor. This one is declared
                                      * as a public constructor to avoid
                                      * certain compiler errors when a copy
                                      * constructor is required even if it is
                                      * not executed (for example when binding
                                      * a temporary object to a constant
                                      * reference). However, if you try to
                                      * actually call it, it will throw an
                                      * exception, since copying finite
                                      * element objects is not really
                                      * supported. If you want to copy such an
                                      * object, use the @p clone function.
                                      */
    FiniteElement (const FiniteElement &);
    
				     /**
				      * Constructor
				      */
    FiniteElement (const FiniteElementData<dim> &fe_data,
		   const std::vector<bool> &restriction_is_additive_flags,
		   const std::vector<std::vector<bool> > &nonzero_components);

				     /**
				      * Virtual destructor. Makes sure
				      * that pointers to this class
				      * are deleted properly.
				      */
    virtual ~FiniteElement ();
    
				     /**
				      * Return a string that uniquely
				      * identifies a finite
				      * element. The general
				      * convention is that this is the
				      * class name, followed by the
				      * space dimension in angle
				      * brackets, and the polynomial
				      * degree and whatever else is
				      * necessary in parentheses. For
				      * example, <tt>FE_Q<2>(3)</tt> is the
				      * value returned for a cubic
				      * element in 2d.
				      *
				      * Systems of elements have their
				      * own naming convention, see the
				      * FESystem class.
				      */
    virtual std::string get_name () const = 0;

				     /**
				      * @name Shape function access
				      * @{
				      */
    
				     /**
				      * Return the value of the
				      * @p ith shape function at the
				      * point @p p. @p p is a point
				      * on the reference element. If
				      * the finite element is
				      * vector-valued, then return the
				      * value of the only non-zero
				      * component of the vector value
				      * of this shape function. If the
				      * shape function has more than
				      * one non-zero component (which
				      * we refer to with the term
				      * non-primitive), then derived
				      * classes implementing this
				      * function should throw an
				      * exception of type
				      * @p ExcShapeFunctionNotPrimitive. In
				      * that case, use the
				      * shape_value_component()
				      * function.
				      *
				      * An
				      * @p ExcUnitShapeValuesDoNotExist
				      * is thrown if the shape values
				      * of the @p FiniteElement under
				      * consideration depends on the
				      * shape of the cell in real
				      * space.
				      */
    virtual double shape_value (const unsigned int  i,
			        const Point<dim>   &p) const;

				     /**
				      * Just like for @p shape_value,
				      * but this function will be
				      * called when the shape function
				      * has more than one non-zero
				      * vector component. In that
				      * case, this function should
				      * return the value of the
				      * @p component-th vector
				      * component of the @p ith shape
				      * function at point @p p.
				      */
    virtual double shape_value_component (const unsigned int i,
					  const Point<dim>   &p,
					  const unsigned int component) const;
    
				     /**
				      * Return the gradient of the
				      * @p ith shape function at the
				      * point @p p. @p p is a point
				      * on the reference element, and
				      * likewise the gradient is the
				      * gradient on the unit cell with
				      * respect to unit cell
				      * coordinates. If
				      * the finite element is
				      * vector-valued, then return the
				      * value of the only non-zero
				      * component of the vector value
				      * of this shape function. If the
				      * shape function has more than
				      * one non-zero component (which
				      * we refer to with the term
				      * non-primitive), then derived
				      * classes implementing this
				      * function should throw an
				      * exception of type
				      * @p ExcShapeFunctionNotPrimitive. In
				      * that case, use the
				      * shape_grad_component()
				      * function.
				      *
				      * An
				      * @p ExcUnitShapeValuesDoNotExist
				      * is thrown if the shape values
				      * of the @p FiniteElement under
				      * consideration depends on the
				      * shape of the cell in real
				      * space.
				      */
    virtual Tensor<1,dim> shape_grad (const unsigned int  i,
				      const Point<dim>   &p) const;

				     /**
				      * Just like for @p shape_grad,
				      * but this function will be
				      * called when the shape function
				      * has more than one non-zero
				      * vector component. In that
				      * case, this function should
				      * return the gradient of the
				      * @p component-th vector
				      * component of the @p ith shape
				      * function at point @p p.
				      */
    virtual Tensor<1,dim> shape_grad_component (const unsigned int i,
						const Point<dim>   &p,
						const unsigned int component) const;

				     /**
				      * Return the tensor of second
				      * derivatives of the @p ith
				      * shape function at point @p p
				      * on the unit cell. The
				      * derivatives are derivatives on
				      * the unit cell with respect to
				      * unit cell coordinates. If
				      * the finite element is
				      * vector-valued, then return the
				      * value of the only non-zero
				      * component of the vector value
				      * of this shape function. If the
				      * shape function has more than
				      * one non-zero component (which
				      * we refer to with the term
				      * non-primitive), then derived
				      * classes implementing this
				      * function should throw an
				      * exception of type
				      * @p ExcShapeFunctionNotPrimitive. In
				      * that case, use the
				      * shape_grad_grad_component()
				      * function.
				      *
				      * An
				      * @p ExcUnitShapeValuesDoNotExist
				      * is thrown if the shape values
				      * of the @p FiniteElement under
				      * consideration depends on the
				      * shape of the cell in real
				      * space.
				      */
    virtual Tensor<2,dim> shape_grad_grad (const unsigned int  i,
					   const Point<dim>   &p) const;

				     /**
				      * Just like for @p shape_grad_grad,
				      * but this function will be
				      * called when the shape function
				      * has more than one non-zero
				      * vector component. In that
				      * case, this function should
				      * return the gradient of the
				      * @p component-th vector
				      * component of the @p ith shape
				      * function at point @p p.
				      */
    virtual Tensor<2,dim> shape_grad_grad_component (const unsigned int i,
						     const Point<dim>   &p,
						     const unsigned int component) const;
				     /**
				      * Check for non-zero values on a face.
				      *
				      * This function returns
				      * @p true, if the shape
				      * function @p shape_index has
				      * non-zero values on the face
				      * @p face_index.
				      */
    virtual bool has_support_on_face (const unsigned int shape_index,
				      const unsigned int face_index) const = 0;
    
				     //@}
				     /**
				      * @name Transfer and constraint matrices
				      * @{
				      */
    
				     /**
				      * Projection from a fine grid
				      * space onto a coarse grid
				      * space. If this projection
				      * operator is associated with a
				      * matrix @p P, then the
				      * restriction of this matrix
				      * @p P_i to a single child cell
				      * is returned here.
				      *
				      * The matrix @p P is the
				      * concatenation or the sum of
				      * the cell matrices @p P_i,
				      * depending on the
				      * @p restriction_is_additive_flags
				      * given to the constructor. This
				      * distinguishes interpolation
				      * (concatenation) and projection
				      * with respect to scalar
				      * products (summation).
				      *
				      * Row and column indices are
				      * related to coarse grid and
				      * fine grid spaces,
				      * respectively, consistent with
				      * the definition of the
				      * associated operator.
				      *
				      * If projection matrices are not
				      * implemented in the derived
				      * finite element class, this
				      * function aborts with
				      * @p ExcProjectionVoid.
				      */
    const FullMatrix<double> &
    get_restriction_matrix (const unsigned int child) const;

				     /**
				      * Embedding matrix between grids.
				      * 
				      * The identity operator from a
				      * coarse grid space into a fine
				      * grid space is associated with
				      * a matrix @p P. The
				      * restriction of this matrix @p P_i to
				      * a single child cell is
				      * returned here.
				      *
				      * The matrix @p P is the
				      * concatenation, not the sum of
				      * the cell matrices
				      * @p P_i. That is, if the same
				      * non-zero entry <tt>j,k</tt> exists
				      * in in two different child
				      * matrices @p P_i, the value
				      * should be the same in both
				      * matrices and it is copied into
				      * the matrix @p P only once.
				      *
				      * Row and column indices are
				      * related to fine grid and
				      * coarse grid spaces,
				      * respectively, consistent with
				      * the definition of the
				      * associated operator.
				      *
				      * These matrices are used by
				      * routines assembling the
				      * prolongation matrix for
				      * multi-level methods.  Upon
				      * assembling the transfer matrix
				      * between cells using this
				      * matrix array, zero elements in
				      * the prolongation matrix are
				      * discarded and will not fill up
				      * the transfer matrix.
				      *
				      * If projection matrices are not
				      * implemented in the derived
				      * finite element class, this
				      * function aborts with
				      * @p ExcEmbeddingVoid. You can
				      * check whether this is the case
				      * by calling the
				      * prolongation_is_implemented().
				      */
    const FullMatrix<double> &
    get_prolongation_matrix (const unsigned int child) const;

                                     /**
                                      * Return whether this element implements
                                      * its prolongation matrices. The return
                                      * value also indicates whether a call to
                                      * the @p get_prolongation_matrix
                                      * function will generate an error or
                                      * not.
                                      *
                                      * This function is mostly here in order
                                      * to allow us to write more efficient
                                      * test programs which we run on all
                                      * kinds of weird elements, and for which
                                      * we simply need to exclude certain
                                      * tests in case something is not
                                      * implemented. It will in general
                                      * probably not be a great help in
                                      * applications, since there is not much
                                      * one can do if one needs these features
                                      * and they are not implemented. This
                                      * function could be used to check
                                      * whether a call to
                                      * <tt>get_prolongation_matrix()</tt> will
                                      * succeed; however, one then still needs
                                      * to cope with the lack of information
                                      * this just expresses.
                                      */
    bool prolongation_is_implemented () const;

                                     /**
                                      * Return whether this element implements
                                      * its restriction matrices. The return
                                      * value also indicates whether a call to
                                      * the @p get_restriction_matrix
                                      * function will generate an error or
                                      * not.
                                      *
                                      * This function is mostly here in order
                                      * to allow us to write more efficient
                                      * test programs which we run on all
                                      * kinds of weird elements, and for which
                                      * we simply need to exclude certain
                                      * tests in case something is not
                                      * implemented. It will in general
                                      * probably not be a great help in
                                      * applications, since there is not much
                                      * one can do if one needs these features
                                      * and they are not implemented. This
                                      * function could be used to check
                                      * whether a call to
                                      * <tt>get_restriction_matrix()</tt> will
                                      * succeed; however, one then still needs
                                      * to cope with the lack of information
                                      * this just expresses.
                                      */
    bool restriction_is_implemented () const;

				     /**
				      * Access the
				      * @p restriction_is_additive_flag
				      * field. See there for more
				      * information on its contents.
				      *
				      * The index must be between zero
				      * and the number of shape
				      * functions of this element.
				      */
    bool restriction_is_additive (const unsigned int index) const;

    				     /**
				      * Return a readonly reference to
				      * the matrix which describes the
				      * constraints at the interface
				      * between a refined and an
				      * unrefined cell.
				      * 
				      * The matrix is obviously empty
				      * in only one space dimension,
				      * since there are no constraints
				      * then.
				      *
				      * Note that some finite elements
				      * do not (yet) implement hanging
				      * node constraints. If this is
				      * the case, then this function
				      * will generate an exception,
				      * since no useful return value
				      * can be generated. If you
				      * should have a way to live with
				      * this, then you might want to
				      * use the
				      * @p constraints_are_implemented
				      * function to check up front
				      * whethehr this function will
				      * succeed or generate the
				      * exception.
				      */
    const FullMatrix<double> & constraints () const;

                                     /**
                                      * Return whether this element
                                      * implements its hanging node
                                      * constraints. The return value
                                      * also indicates whether a call
                                      * to the @p constraint function
                                      * will generate an error or not.
                                      *
                                      * This function is mostly here
                                      * in order to allow us to write
                                      * more efficient test programs
                                      * which we run on all kinds of
                                      * weird elements, and for which
                                      * we simply need to exclude
                                      * certain tests in case hanging
                                      * node constraints are not
                                      * implemented. It will in
                                      * general probably not be a
                                      * great help in applications,
                                      * since there is not much one
                                      * can do if one needs hanging
                                      * node constraints and they are
                                      * not implemented. This function
                                      * could be used to check whether
                                      * a call to <tt>constraints()</tt>
                                      * will succeed; however, one
                                      * then still needs to cope with
                                      * the lack of information this
                                      * just expresses.
                                      */
    bool constraints_are_implemented () const;

				     /**
				      * Return the matrix
				      * interpolating from the given
				      * finite element to the present
				      * one. The size of the matrix is
				      * then @p dofs_per_cell times
				      * <tt>source.dofs_per_cell</tt>.
				      *
				      * Derived elements will have to
				      * implement this function. They
				      * may only provide interpolation
				      * matrices for certain source
				      * finite elements, for example
				      * those from the same family. If
				      * they don't implement
				      * interpolation from a given
				      * element, then they must throw
				      * an exception of type
				      * FiniteElement<dim>::ExcInterpolationNotImplemented.
				      */
    virtual void
    get_interpolation_matrix (const FiniteElement<dim> &source,
			      FullMatrix<double>           &matrix) const;
				     //@}
    
				     /**
				      * Comparison operator. We also
				      * check for equality of the
				      * constraint matrix, which is
				      * quite an expensive operation.
				      * Do therefore use this function
				      * with care, if possible only
				      * for debugging purposes.
				      *
				      * Since this function is not
				      * that important, we avoid an
				      * implementational question
				      * about comparing arrays and do
				      * not compare the matrix arrays
				      * @p restriction and
				      * @p prolongation.
				      */
    bool operator == (const FiniteElement<dim> &) const;

				     /**
				      * @name Index computations
				      * @{
				      */
				     /**
				      * Compute vector component and
				      * index of this shape function
				      * within the shape functions
				      * corresponding to this
				      * component from the index of a
				      * shape function within this
				      * finite element.
				      *
				      * If the element is scalar, then
				      * the component is always zero,
				      * and the index within this
				      * component is equal to the
				      * overall index.
				      *
				      * If the shape function
				      * referenced has more than one
				      * non-zero component, then it
				      * cannot be associated with one
				      * vector component, and an
				      * exception of type
				      * @p ExcShapeFunctionNotPrimitive
				      * will be raised.
				      *
				      * Note that if the element is
				      * composed of other (base)
				      * elements, and a base element
				      * has more than one component
				      * but all its shape functions
				      * are primitive (i.e. are
				      * non-zero in only one
				      * component), then this mapping
				      * contains valid
				      * information. However, the
				      * index of a shape function of
				      * this element within one
				      * component (i.e. the second
				      * number of the respective entry
				      * of this array) does not
				      * indicate the index of the
				      * respective shape function
				      * within the base element (since
				      * that has more than one
				      * vector-component). For this
				      * information, refer to the
				      * @p system_to_base_table field
				      * and the
				      * @p system_to_base_index
				      * function.
				      */
    std::pair<unsigned int, unsigned int>
    system_to_component_index (const unsigned int index) const;

				     /**
				      * Compute the shape function for
				      * the given vector component and
				      * index.
				      *
				      * If the element is scalar, then
				      * the component must be zero,
				      * and the index within this
				      * component is equal to the
				      * overall index.
				      *
				      * This is the opposite operation
				      * from the @p system_to_component_index
				      * function.
				      */
   unsigned int component_to_system_index(const unsigned int component,
                                          const unsigned int index) const;
  
				     /**
				      * Same as above, but do it for
				      * shape functions and their
				      * indices on a face.
				      */
    std::pair<unsigned int, unsigned int>
    face_system_to_component_index (const unsigned int index) const;

                                     /**
                                      * Return for shape function
                                      * @p index the base element it
                                      * belongs to, the number of the
                                      * copy of this base element
                                      * (which is between zero and the
                                      * multiplicity of this element),
                                      * and the index of this shape
                                      * function within this base
                                      * element.
                                      *
                                      * If the element is not composed of
				      * others, then base and instance
				      * are always zero, and the index
				      * is equal to the number of the
				      * shape function. If the element
				      * is composed of single
				      * instances of other elements
				      * (i.e. all with multiplicity
				      * one) all of which are scalar,
				      * then base values and dof
				      * indices within this element
				      * are equal to the
				      * @p system_to_component_table. It
				      * differs only in case the
				      * element is composed of other
				      * elements and at least one of
				      * them is vector-valued itself.
				      *
				      * This function returns valid
				      * values also in the case of
				      * vector-valued
				      * (i.e. non-primitive) shape
				      * functions, in contrast to the
				      * @p system_to_component_index
				      * function.
                                      */
    std::pair<std::pair<unsigned int, unsigned int>, unsigned int>
    system_to_base_index (const unsigned int index) const;

                                     /**
                                      * Same as
                                      * @p system_to_base_index, but
                                      * for degrees of freedom located
                                      * on a face.
                                      */
    std::pair<std::pair<unsigned int, unsigned int>, unsigned int>
    face_system_to_base_index (const unsigned int index) const;
    
				     /**
				      * Return in which of the vector
				      * components of this finite
				      * element the @p ithe shape
				      * function is non-zero. The
				      * length of the returned array
				      * is equal to the number of
				      * vector components of this
				      * element.
				      *
				      * For most finite element
				      * spaces, the result of this
				      * function will be a vector with
				      * exactly one element being
				      * @p true, since for most
				      * spaces the individual vector
				      * components are independent. In
				      * that case, the component with
				      * the single zero is also the
				      * first element of what
				      * <tt>system_to_component_index(i)</tt>
				      * returns.
				      *
				      * Only for those
				      * spaces that couple the
				      * components, for example to
				      * make a shape function
				      * divergence free, will there be
				      * more than one @p true entry.
				      */
    const std::vector<bool> &
    get_nonzero_components (const unsigned int i) const;

				     /**
				      * Return in how many vector
				      * components the @p ith shape
				      * function is non-zero. This
				      * value equals the number of
				      * entries equal to @p true in
				      * the result of the
				      * @p get_nonzero_components
				      * function.
				      *
				      * For most finite element
				      * spaces, the result will be
				      * equal to one. It is not equal
				      * to one only for those ansatz
				      * spaces for which vector-valued
				      * shape functions couple the
				      * individual components, for
				      * example in order to make them
				      * divergence-free.
				      */
    unsigned int
    n_nonzero_components (const unsigned int i) const;

				     /**
				      * Return whether the @p ith
				      * shape function is primitive in
				      * the sense that the shape
				      * function is non-zero in only
				      * one vector
				      * component. Non-primitive shape
				      * functions would then, for
				      * example, be those of
				      * divergence free ansatz spaces,
				      * in which the individual vector
				      * components are coupled.
				      *
				      * The result of the function is
				      * @p true if and only if the
				      * result of
				      * <tt>n_nonzero_components(i)</tt> is
				      * equal to one.
				      */
    bool
    is_primitive (const unsigned int i) const;

				     /**
				      * Return whether the entire
				      * finite element is primitive,
				      * in the sense that all its
				      * shape functions are
				      * primitive. If the finite
				      * element is scalar, then this
				      * is always the case.
				      *
				      * Since this is an extremely
				      * common operation, the result
				      * is cached in the
				      * @p cached_primitivity
				      * variable which is computed in
				      * the constructor.
				      */
    bool
    is_primitive () const;
    
				     /**
				      * Number of base elements in a
				      * mixed discretization.
				      *
				      * Note that even for vector
				      * valued finite elements, the
				      * number of components needs not
				      * coincide with the number of
				      * base elements, since they may
				      * be reused. For example, if you
				      * create a FESystem with
				      * three identical finite element
				      * classes by using the
				      * constructor that takes one
				      * finite element and a
				      * multiplicity, then the number
				      * of base elements is still one,
				      * although the number of
				      * components of the finite
				      * element is equal to the
				      * multiplicity.
				      */
    virtual unsigned int n_base_elements () const = 0;

				     /**
				      * Access to base element
				      * objects. If the element is
				      * scalar, then
				      * <tt>base_element(0)</tt> is
				      * @p this.
				      */
    virtual
    const FiniteElement<dim> &
    base_element (const unsigned int index) const = 0;

                                     /**
                                      * This index denotes how often
                                      * the base element @p index is
                                      * used in a composed element. If
                                      * the element is scalar, then
                                      * the result is always equal to
                                      * one. See the documentation for
                                      * the @p n_base_elements
                                      * function for more details.
                                      */
    virtual
    unsigned int
    element_multiplicity (const unsigned int index) const = 0;
    
 				     /**
				      * Given a vector component,
				      * return an index which base
				      * element implements this
				      * component, and which vector
				      * component in this base element
				      * this is. This information is
				      * only of interest for
				      * vector-valued finite elements
				      * which are composed of several
				      * sub-elements. In that case,
				      * one may want to obtain
				      * information about the element
				      * implementing a certain vector
				      * component, which can be done
				      * using this function and the
				      * FESystem::@p base_element
				      * function.
				      *
				      * If this is a scalar finite
				      * element, then the return value
				      * is always equal to a pair of
				      * zeros.
				      */
    std::pair<unsigned int,unsigned int>
    component_to_base (const unsigned int component) const;
				     //@}
    
				     /**
				      * @name Support points and interpolation
				      * @{
				      */
    
				     /**
				      * Return the support points of
				      * the trial functions on the
				      * unit cell, if the derived
				      * finite element defines some.
				      * Finite elements that allow
				      * some kind of interpolation
				      * operation usually have support
				      * points. On the other hand,
				      * elements that define their
				      * degrees of freedom by, for
				      * example, moments on faces, or
				      * as derivatives, don't have
				      * support points. In that case,
				      * the returned field is empty.
				      *
				      * If the finite element defines
				      * support points, then their
				      * number equals the number of
				      * degrees of freedom of the
				      * element.  The order of points
				      * in the array matches that
				      * returned by the
				      * <tt>cell->get_dof_indices</tt>
				      * function.
				      *
				      * See the class documentation
				      * for details on support points.
				      */
    const std::vector<Point<dim> > &
    get_unit_support_points () const;    

				     /**
				      * Return whether a finite
				      * element has defined support
				      * points. If the result is true,
				      * then a call to the
				      * @p get_unit_support_points
				      * yields a non-empty array.
				      *
				      * The result may be false if an
				      * element is not defined by
				      * interpolating shape functions,
				      * for example by P-elements on
				      * quadrilaterals. It will
				      * usually only be true if the
				      * element constructs its shape
				      * functions by the requirement
				      * that they be one at a certain
				      * point and zero at all the
				      * points associated with the
				      * other shape functions.
				      *
				      * In composed elements (i.e. for
				      * the FESystem class, the
				      * result will be true if all all
				      * the base elements have defined
				      * support points.
				      */
    bool has_support_points () const;

                                     /**
                                      * Return the position of the
                                      * support point of the
                                      * @p indexth shape function. If
                                      * it does not exist, raise an
                                      * exception.
                                      *
                                      * The default implementation
                                      * simply returns the respective
                                      * element from the array you get
                                      * from
                                      * get_unit_support_points(),
                                      * but derived elements may
                                      * overload this function. In
                                      * particular, note that the
                                      * FESystem class overloads
                                      * it so that it can return the
                                      * support points of individual
                                      * base elements, of not all the
                                      * base elements define support
                                      * points. In this way, you can
                                      * still ask for certain support
                                      * points, even if
                                      * @p get_unit_support_points
                                      * only returns an empty array.
                                      */
    virtual
    Point<dim>
    unit_support_point (const unsigned int index) const;
    
				     /**
				      * Return the support points of
				      * the trial functions on the
				      * unit face, if the derived
				      * finite element defines some.
				      * Finite elements that allow
				      * some kind of interpolation
				      * operation usually have support
				      * points. On the other hand,
				      * elements that define their
				      * degrees of freedom by, for
				      * example, moments on faces, or
				      * as derivatives, don't have
				      * support points. In that case,
				      * the returned field is empty
				      *
				      * Note that elements that have
				      * support points need not
				      * necessarily have some on the
				      * faces, even if the
				      * interpolation points are
				      * located physically on a
				      * face. For example, the
				      * discontinuous elements have
				      * interpolation points on the
				      * vertices, and for higher
				      * degree elements also on the
				      * faces, but they are not
				      * defined to be on faces since
				      * in that case degrees of
				      * freedom from both sides of a
				      * face (or from all adjacent
				      * elements to a vertex) would be
				      * identified with each other,
				      * which is not what we would
				      * like to have). Logically,
				      * these degrees of freedom are
				      * therefore defined to belong to
				      * the cell, rather than the face
				      * or vertex. In that case, the
				      * returned element would
				      * therefore have length zero.
				      *
				      * If the finite element defines
				      * support points, then their
				      * number equals the number of
				      * degrees of freedom on the face
				      * (@p dofs_per_face). The order
				      * of points in the array matches
				      * that returned by the
				      * <tt>cell->get_dof_indices</tt>
				      * function.
				      *
				      * See the class documentation
				      * for details on support points.
				      */
    const std::vector<Point<dim-1> > &
    get_unit_face_support_points () const;    

				     /**
				      * Return whether a finite
				      * element has defined support
				      * points on faces. If the result
				      * is true, then a call to the
				      * @p get_unit_support_points
				      * yields a non-empty array.
				      *
				      * For more information, see the
				      * documentation for the
				      * has_support_points()
				      * function.
				      */
    bool has_face_support_points () const;

                                     /**
                                      * The function corresponding to
                                      * the unit_support_point()
                                      * function, but for faces. See
                                      * there for more information.
                                      */
    virtual
    Point<dim-1>
    unit_face_support_point (const unsigned int index) const;
    
				     /**
				      * Return a support point vector
				      * for generalized interpolation.
				      */
    const std::vector<Point<dim> > &
    get_generalized_support_points () const;    

				     /**
				      *
				      */
    bool has_generalized_support_points () const;

				     /**
				      *
				      */
    const std::vector<Point<dim-1> > &
    get_generalized_face_support_points () const;

				     /**
				      * Return whether a finite
				      * element has defined support
				      * points on faces. If the result
				      * is true, then a call to the
				      * @p get_unit_support_points
				      * yields a non-empty array.
				      *
				      * For more information, see the
				      * documentation for the
				      * has_support_points()
				      * function.
				      */
    bool has_generalized_face_support_points () const;

				     /**
				      * Interpolate a set of scalar
				      * values, computed in the
				      * generalized support points.
				      *
				      * @note This function is
				      * implemented in
				      * FiniteElement for the case
				      * that the element has support
				      * points. In this case, the
				      * resulting coefficients are
				      * just the values in the suport
				      * points. All other elements
				      * must reimplement it.
				      */
    virtual void interpolate(std::vector<double>&       local_dofs,
			     const std::vector<double>& values) const;
      
				     /**
				      * Interpolate a set of vector
				      * values, computed in the
				      * generalized support points.
				      *
				      * Since a finite element often
				      * only interpolates part of a
				      * vector, <tt>offset</tt> is
				      * used to determine the first
				      * component of the vector to be
				      * interpolated. Maybe consider
				      * changing your data structures
				      * to use the next function.
				      */
    virtual void interpolate(std::vector<double>&                local_dofs,
			     const std::vector<Vector<double> >& values,
			     unsigned int offset = 0) const;
      
				     /**
				      * Interpolate a set of vector
				      * values, computed in the
				      * generalized support points.
				      */
    virtual void interpolate(
      std::vector<double>& local_dofs,
      const VectorSlice<const std::vector<std::vector<double> > >& values) const;
      
				     //@}
    
				     /**
				      * Determine an estimate for the
				      * memory consumption (in bytes)
				      * of this object.
				      *
				      * This function is made virtual,
				      * since finite element objects
				      * are usually accessed through
				      * pointers to their base class,
				      * rather than the class itself.
				      */
    virtual unsigned int memory_consumption () const;
				     /**
				      * Exception
				      *
				      * @ingroup Exceptions
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
				      *
				      * @ingroup Exceptions
				      */
    DeclException0 (ExcFENotPrimitive);
				     /**
				      * Exception
				      *
				      * @ingroup Exceptions
				      */
    DeclException0 (ExcUnitShapeValuesDoNotExist);

				     /**
				      * Attempt to access support
				      * points of a finite element
				      * which is not Lagrangian.
				      *
				      * @ingroup Exceptions
				      */
    DeclException0 (ExcFEHasNoSupportPoints);

				     /**
				      * Attempt to access embedding
				      * matrices of a finite element
				      * which did not implement these
				      * matrices.
				      *
				      * @ingroup Exceptions
				      */
    DeclException0 (ExcEmbeddingVoid);
    
				     /**
				      * Attempt to access restriction
				      * matrices of a finite element
				      * which did not implement these
				      * matrices.
				      *
				      * Exception
				      * @ingroup Exceptions
				      */
    DeclException0 (ExcProjectionVoid);
    
				     /**
				      * Attempt to access constraint
				      * matrices of a finite element
				      * which did not implement these
				      * matrices.
				      *
				      * Exception
				      * @ingroup Exceptions
				      */
    DeclException0 (ExcConstraintsVoid);
    
				     /**
				      * Exception
				      * @ingroup Exceptions
				      */
    DeclException2 (ExcWrongInterfaceMatrixSize,
		    int, int,
		    << "The interface matrix has a size of " << arg1
		    << "x" << arg2
		    << ", which is not reasonable in the present dimension.");
				     /**
				      * Exception
				      * @ingroup Exceptions
				      */
    DeclException2 (ExcComponentIndexInvalid,
		    int, int,
		    << "The component-index pair (" << arg1 << ", " << arg2
		    << ") is invalid, i.e. non-existent");
                                     /**
                                      * Exception
				      * @ingroup Exceptions
                                      */
    DeclException0 (ExcInterpolationNotImplemented);
    
  protected:  
 				     /**
				      * Array of projection matrices. See
				      * get_restriction_matrix() above.
				      *
				      * Matrices in this array are
				      * automatically initialized to
				      * correct size. If the derived
				      * finite element class does not
				      * implement these matrices, they
				      * should be resized to zero
				      * size.
				      */
    FullMatrix<double> restriction[GeometryInfo<dim>::children_per_cell];

    				     /**
				      * Array of embedding matrices. See
				      * <tt>get_prolongation_matrix()</tt> above.
				      *
				      * Matrices in this array are
				      * automatically initialized to
				      * correct size. If the derived
				      * finite element class does not
				      * implement these matrices, they
				      * should be resized to zero
				      * size.
				      */
    FullMatrix<double> prolongation[GeometryInfo<dim>::children_per_cell];

   				     /**
				      * Specify the constraints which
				      * the dofs on the two sides of a
				      * cell interface underly if the
				      * line connects two cells of
				      * which one is refined once.
				      *
				      * For further details see the
				      * general description of the
				      * derived class.
				      *
				      * This field is obviously
				      * useless in one space dimension
				      * and has there a zero size.
				      */
    FullMatrix<double> interface_constraints;

                                     /**
                                      * Return the size of interface
                                      * constraint matrices. Since
                                      * this is needed in every
                                      * derived finite element class
                                      * when initializing their size,
                                      * it is placed into this
                                      * function, to avoid having to
                                      * recompute the
                                      * dimension-dependent size of
                                      * these matrices each time.
                                      *
                                      * Note that some elements do not
                                      * implement the interface
                                      * constraints for certain
                                      * polynomial degrees. In this
                                      * case, this function still
                                      * returns the size these
                                      * matrices should have when
                                      * implemented, but the actual
                                      * matrices are empty.
                                      */
    TableIndices<2>
    interface_constraints_size () const;
    
				     /**
				      * Store what
				      * @p system_to_component_index
				      * will return.
				      */
    std::vector< std::pair<unsigned int, unsigned int> > system_to_component_table;

                                     /**
  				      * Map between linear dofs and
 				      * component dofs on face. This
 				      * is filled with default values
 				      * in the constructor, but
 				      * derived classes will have to
 				      * overwrite the information if
 				      * necessary.
 				      *
 				      * By component, we mean the
 				      * vector component, not the base
 				      * element. The information thus
 				      * makes only sense if a shape
 				      * function is non-zero in only
 				      * one component.
				      */
    std::vector< std::pair<unsigned int, unsigned int> > face_system_to_component_table;

				     /**
				      * For each shape function, store
				      * to which base element and
				      * which instance of this base
				      * element (in case its
				      * multiplicity is greater than
				      * one) it belongs, and its index
				      * within this base element. If
				      * the element is not composed of
				      * others, then base and instance
				      * are always zero, and the index
				      * is equal to the number of the
				      * shape function. If the element
				      * is composed of single
				      * instances of other elements
				      * (i.e. all with multiplicity
				      * one) all of which are scalar,
				      * then base values and dof
				      * indices within this element
				      * are equal to the
				      * @p system_to_component_table. It
				      * differs only in case the
				      * element is composed of other
				      * elements and at least one of
				      * them is vector-valued itself.
				      *
				      * This array has valid values
				      * also in the case of
				      * vector-valued
				      * (i.e. non-primitive) shape
				      * functions, in contrast to the
				      * @p system_to_component_table.
				      */
    std::vector<std::pair<std::pair<unsigned int,unsigned int>,unsigned int> >
    system_to_base_table;

				     /**
				      * Likewise for the indices on
				      * faces.
				      */
    std::vector<std::pair<std::pair<unsigned int,unsigned int>,unsigned int> >
    face_system_to_base_table;
    
				     /**
				      * The base element establishing
				      * a component.
				      *
				      * This table converts a
				      * component number to a pair
				      * consisting of the
				      * @p base_element number, and
				      * the component within this base
				      * element. While component
				      * information contains
				      * multiplicity of base elements,
				      * the result allows access to
				      * shape functions of the base
				      * element.
				      *
				      * This variable is set to the
				      * correct size by the
				      * constructor of this class, but
				      * needs to be initialized by
				      * derived classes, unless its
				      * size is one and the only entry
				      * is a zero, which is the case
				      * for scalar elements. In that
				      * case, the initialization by
				      * the base class is sufficient.
				      */
    std::vector<std::pair<unsigned int, unsigned int> > component_to_base_table;
    
				     /**
				      * Projection matrices are
				      * concatenated or summed up.
				      *
				      * This flags decides on how the
				      * projection matrices of the
				      * children of the same father
				      * are put together to one
				      * operator. The possible modes
				      * are concatenation and
				      * summation.
				      *
				      * If the projection is defined
				      * by an interpolation operator,
				      * the child matrices are
				      * concatenated, i.e. values
				      * belonging to the same node
				      * functional are identified and
				      * enter the interpolated value
				      * only once. In this case, the
				      * flag must be @p false.
				      *
				      * For projections with respect
				      * to scalar products, the child
				      * matrices must be summed up to
				      * build the complete matrix. The
				      * flag should be @p true.
				      *
				      * For examples of use of these
				      * flags, see the places in the
				      * library where it is queried.
				      * 
				      * There is one flag per shape
				      * function, indicating whether
				      * it belongs to the class of
				      * shape functions that are
				      * additive in the restriction or
				      * not.
				      *
				      * Note that in previous versions
				      * of the library, there was one
				      * flag per vector component of
				      * the element. This is based on
				      * the fact that all the shape
				      * functions that belong to the
				      * same vector component must
				      * necessarily behave in the same
				      * way, to make things
				      * reasonable. However, the
				      * problem is that it is
				      * sometimes impossible to query
				      * this flag in the vector-valued
				      * case: this used to be done
				      * with the
				      * @p system_to_component_index
				      * function that returns which
				      * vector component a shape
				      * function is associated
				      * with. The point is that since
				      * we now support shape functions
				      * that are associated with more
				      * than one vector component (for
				      * example the shape functions of
				      * Raviart-Thomas, or Nedelec
				      * elements), that function can
				      * no more be used, so it can be
				      * difficult to find out which
				      * for vector component we would
				      * like to query the
				      * restriction-is-additive flags.
				      */
    const std::vector<bool> restriction_is_additive_flags;

				     /**
				      * List of support points on the
				      * unit cell, in case the finite
				      * element has any. The
				      * constructor leaves this field
				      * empty, derived classes may
				      * write in some contents.
				      *
				      * Finite elements that allow
				      * some kind of interpolation
				      * operation usually have support
				      * points. On the other hand,
				      * elements that define their
				      * degrees of freedom by, for
				      * example, moments on faces, or
				      * as derivatives, don't have
				      * support points. In that case,
				      * this field remains empty.
				      */
    std::vector<Point<dim> > unit_support_points;

				     /**
				      * Same for the faces. See the
				      * description of the
				      * @p get_unit_face_support_points
				      * function for a discussion of
				      * what contributes a face
				      * support point.
				      */
    std::vector<Point<dim-1> > unit_face_support_points;
    
				     /**
				      * Support points used for
				      * interpolation functions of
				      * non-Lagrangian elements.
				      */
    std::vector<Point<dim> > generalized_support_points;
    
				     /**
				      * Face support points used for
				      * interpolation functions of
				      * non-Lagrangian elements.
				      */    
    std::vector<Point<dim-1> > generalized_face_support_points;
    
				     /**
				      * For each shape function, give
				      * a vector of bools (with size
				      * equal to the number of vector
				      * components which this finite
				      * element has) indicating in
				      * which component each of these
				      * shape functions is non-zero.
				      *
				      * For primitive elements, there
				      * is only one non-zero
				      * component.
				      */
    const std::vector<std::vector<bool> > nonzero_components;

				     /**
				      * This array holds how many
				      * values in the respective entry
				      * of the @p nonzero_components
				      * element are non-zero. The
				      * array is thus a short-cut to
				      * allow faster access to this
				      * information than if we had to
				      * count the non-zero entries
				      * upon each request for this
				      * information. The field is
				      * initialized in the constructor
				      * of this class.
				      */
    const std::vector<unsigned int> n_nonzero_components_table;

				     /**
				      * Store whether all shape
				      * functions are primitive. Since
				      * finding this out is a very
				      * common operation, we cache the
				      * result, i.e. compute the value
				      * in the constructor for simpler
				      * access.
				      */
    const bool cached_primitivity;

                                     /**
				      * Compute second derivatives by
				      * finite differences of
				      * gradients.
				      */
    void compute_2nd (const Mapping<dim>                      &mapping,
		      const typename Triangulation<dim>::cell_iterator    &cell,
		      const unsigned int                       offset,
		      typename Mapping<dim>::InternalDataBase &mapping_internal,
		      InternalDataBase                        &fe_internal,
		      FEValuesData<dim>                       &data) const;

				     /**
				      * Given the pattern of nonzero
				      * components for each shape
				      * function, compute for each
				      * entry how many components are
				      * non-zero for each shape
				      * function. This function is
				      * used in the constructor of
				      * this class.
				      */
    static
    std::vector<unsigned int>
    compute_n_nonzero_components (const std::vector<std::vector<bool> > &nonzero_components);
    
				     /**
				      * Allow the FESystem class to access the
				      * restriction and prolongation matrices
				      * directly. Hence, FESystem has the
				      * possibility to see if these matrices
				      * are initialized without accessing
				      * these matrices through the
				      * @p get_restriction_matrix and
				      * @p get_prolongation_matrix
				      * functions. This is important as these
				      * functions include assertions that
				      * throw if the matrices are not already
				      * initialized.
				      */
    template <int dim_> friend class FESystem;

                                     /**
                                      * Make the inner class a
                                      * friend. This is not strictly
                                      * necessary, but the Intel
                                      * compiler seems to want this.
                                      */
    friend class InternalDataBase;
    

				     /**
				      * Exception
				      *
				      * @ingroup Exceptions
				      */
    DeclException0 (ExcBoundaryFaceUsed);
				     /**
				      * Exception
				      *
				      * @ingroup Exceptions
				      */
    DeclException0 (ExcJacobiDeterminantHasWrongSign);

  protected:

				     /**
				      * Determine the values a finite
				      * element should compute on
				      * initialization of data for
				      * @p FEValues.
				      *
				      * Given a set of flags
				      * indicating what quantities are
				      * requested from a @p FEValues
				      * object, @p update_once and
				      * @p update_each compute which
				      * values must really be
				      * computed. Then, the
				      * <tt>fill_*_values</tt> functions
				      * are called with the result of
				      * these.
				      *
				      * Furthermore, values must be
				      * computed either on the unit
				      * cell or on the physical
				      * cell. For instance, the
				      * function values of @p FE_Q do
				      * only depend on the quadrature
				      * points on the unit
				      * cell. Therefore, this flags
				      * will be returned by
				      * @p update_once. The gradients
				      * require computation of the
				      * covariant transformation
				      * matrix. Therefore,
				      * @p update_covariant_transformation
				      * and @p update_gradients will
				      * be returned by
				      * @p update_each.
				      *
				      * For an example see the same
				      * function in the derived class
				      * @p FE_Q.
				      */
    virtual UpdateFlags update_once (const UpdateFlags flags) const = 0;
  
				     /**
				      * Complementary function for
				      * @p update_once.
				      *
				      * While @p update_once returns
				      * the values to be computed on
				      * the unit cell for yielding the
				      * required data, this function
				      * determines the values that
				      * must be recomputed on each
				      * cell.
				      *
				      * Refer to @p update_once for
				      * more details.
				      */
    virtual UpdateFlags update_each (const UpdateFlags flags) const = 0;
  
				     /**
				      * @p clone function instead of
				      * a copy constructor.
				      *
				      * This function is needed by the
				      * constructors of FESystem as well
				      * as by the FECollection class.
				      */
    virtual FiniteElement<dim> *clone() const = 0;
    
  private:
				     /**
				      * Second derivatives of shapes
				      * functions are not computed
				      * analytically, but by finite
				      * differences of the
				      * gradients. This static
				      * variable denotes the step
				      * length to be used for
				      * that. It's value is set to
				      * 1e-6.
				      */
    static const double fd_step_length;

				     /**
				      * Prepare internal data
				      * structures and fill in values
				      * independent of the
				      * cell. Returns a pointer to an
				      * object of which the caller of
				      * this function then has to
				      * assume ownership (which
				      * includes destruction when it
				      * is no more needed).
				      */
    virtual typename Mapping<dim>::InternalDataBase*
    get_data (const UpdateFlags      flags,
	      const Mapping<dim>    &mapping,
	      const Quadrature<dim> &quadrature) const = 0;

				     /**
				      * Prepare internal data
				      * structure for transformation
				      * of faces and fill in values
				      * independent of the
				      * cell. Returns a pointer to an
				      * object of which the caller of
				      * this function then has to
				      * assume ownership (which
				      * includes destruction when it
				      * is no more needed).
				      */
    virtual typename Mapping<dim>::InternalDataBase*
    get_face_data (const UpdateFlags        flags,
		   const Mapping<dim>      &mapping,
		   const Quadrature<dim-1> &quadrature) const;

				     /**
				      * Prepare internal data
				      * structure for transformation
				      * of children of faces and fill
				      * in values independent of the
				      * cell. Returns a pointer to an
				      * object of which the caller of
				      * this function then has to
				      * assume ownership (which
				      * includes destruction when it
				      * is no more needed).
				      */
    virtual typename Mapping<dim>::InternalDataBase*
    get_subface_data (const UpdateFlags        flags,
		      const Mapping<dim>      &mapping,
		      const Quadrature<dim-1> &quadrature) const;

				     /**
				      * Fill the fields of
				      * FEValues. This function
				      * performs all the operations
				      * needed to compute the data of an
				      * @p FEValues object.
				      *
				      * The same function in
				      * @p mapping must have been
				      * called for the same cell first!
				      */				      
    virtual void
    fill_fe_values (const Mapping<dim>                   &mapping,
		    const typename Triangulation<dim>::cell_iterator &cell,
		    const Quadrature<dim>                &quadrature,
		    typename Mapping<dim>::InternalDataBase       &mapping_internal,
		    typename Mapping<dim>::InternalDataBase       &fe_internal,
		    FEValuesData<dim>                    &data) const = 0;
    
				     /**
				      * Fill the fields of
				      * FEFaceValues. This function
				      * performs all the operations
				      * needed to compute the data of an
				      * @p FEFaceValues object.
				      *
				      * The same function in
				      * @p mapping must have been
				      * called for the same cell first!
				      */				      
    virtual void
    fill_fe_face_values (const Mapping<dim>                   &mapping,
			 const typename Triangulation<dim>::cell_iterator &cell,
			 const unsigned int                    face_no,
			 const Quadrature<dim-1>              &quadrature,
			 typename Mapping<dim>::InternalDataBase       &mapping_internal,
			 typename Mapping<dim>::InternalDataBase       &fe_internal,
			 FEValuesData<dim>                    &data) const = 0;
    
				     /**
				      * Fill the fields of
				      * FESubfaceValues. This function
				      * performs all the operations
				      * needed to compute the data of an
				      * @p FESubfaceValues object.
				      *
				      * The same function in
				      * @p mapping must have been
				      * called for the same cell first!
				      */				      
    virtual void
    fill_fe_subface_values (const Mapping<dim>                   &mapping,
			    const typename Triangulation<dim>::cell_iterator &cell,
			    const unsigned int                    face_no,
			    const unsigned int                    sub_no,
			    const Quadrature<dim-1>              &quadrature,
			    typename Mapping<dim>::InternalDataBase &mapping_internal,
			    typename Mapping<dim>::InternalDataBase &fe_internal,
			    FEValuesData<dim>                    &data) const = 0;

				     /**
				      * Declare some other classes as
				      * friends of this class.
				      */
    friend class FEValuesBase<dim>;
    friend class FEValues<dim>;
    friend class FEFaceValues<dim>;
    friend class FESubfaceValues<dim>;
    friend class FESystem<dim>;
    friend class FECollection<dim>;
};

/*@}*/
//----------------------------------------------------------------------//

template <int dim>  
inline
std::pair<unsigned int,unsigned int>
FiniteElement<dim>::system_to_component_index (const unsigned int index) const
{
  Assert (index < system_to_component_table.size(),
	 ExcIndexRange(index, 0, system_to_component_table.size()));
  Assert (is_primitive (index),
	  typename FiniteElement<dim>::ExcShapeFunctionNotPrimitive(index));
  return system_to_component_table[index];
}

template <int dim>  
inline
unsigned int
FiniteElement<dim>::component_to_system_index (const unsigned int component,
                                                   const unsigned int index) const
{
   std::vector< std::pair<unsigned int, unsigned int> >::const_iterator
      it = std::find(system_to_component_table.begin(), system_to_component_table.end(),
                     std::pair<unsigned int, unsigned int>(component, index));

   Assert(it != system_to_component_table.end(), ExcComponentIndexInvalid(component, index));
   return std::distance(system_to_component_table.begin(), it);
}



template <int dim>  
inline
std::pair<unsigned int,unsigned int>
FiniteElement<dim>::face_system_to_component_index (const unsigned int index) const
{
  Assert(index < face_system_to_component_table.size(),
	 ExcIndexRange(index, 0, face_system_to_component_table.size()));

                                   // in debug mode, check whether the
                                   // function is primitive, since
                                   // otherwise the result may have no
                                   // meaning
                                   //
                                   // since the primitivity tables are
                                   // all geared towards cell dof
                                   // indices, rather than face dof
                                   // indices, we have to work a
                                   // little bit...
                                   //
                                   // in 1d, the face index is equal
                                   // to the cell index
  Assert (((dim == 1) && is_primitive(index))
          ||
                                           // in 2d, construct it like
                                           // this:
          ((dim == 2) &&
           is_primitive (index < (GeometryInfo<2>::vertices_per_face *
                                  this->dofs_per_vertex)
                         ?
                         index
                         :
                         GeometryInfo<2>::vertices_per_cell *
                         this->dofs_per_vertex +
                         (index -
                          GeometryInfo<2>::vertices_per_face *
                          this->dofs_per_vertex)))
          ||
                                           // likewise in 3d, but more
                                           // complicated
          ((dim == 3) &&
           is_primitive (index < (GeometryInfo<3>::vertices_per_face *
                                  this->dofs_per_vertex)
                         ?
                         index
                         :
                         (index < (GeometryInfo<3>::vertices_per_face *
                                   this->dofs_per_vertex
                                   +
                                   GeometryInfo<3>::lines_per_face *
                                   this->dofs_per_line)
                          ?
                          GeometryInfo<3>::vertices_per_cell *
                          this->dofs_per_vertex +
                          (index -
                           GeometryInfo<3>::vertices_per_face *
                           this->dofs_per_vertex)
                          :
                          GeometryInfo<3>::vertices_per_cell *
                          this->dofs_per_vertex +
                          GeometryInfo<3>::lines_per_cell *
                          this->dofs_per_line +
                          (index -
                           GeometryInfo<3>::vertices_per_face *
                           this->dofs_per_vertex
                           -
                           GeometryInfo<3>::lines_per_face *
                           this->dofs_per_line)))),
          typename FiniteElement<dim>::ExcShapeFunctionNotPrimitive(index));

  return face_system_to_component_table[index];
}



template <int dim>  
inline
std::pair<std::pair<unsigned int,unsigned int>,unsigned int>
FiniteElement<dim>::system_to_base_index (const unsigned int index) const
{
  Assert (index < system_to_base_table.size(),
	 ExcIndexRange(index, 0, system_to_base_table.size()));
  return system_to_base_table[index];
}




template <int dim>  
inline
std::pair<std::pair<unsigned int,unsigned int>,unsigned int>
FiniteElement<dim>::face_system_to_base_index (const unsigned int index) const
{
  Assert(index < face_system_to_base_table.size(),
	 ExcIndexRange(index, 0, face_system_to_base_table.size()));
  return face_system_to_base_table[index];
}



template <int dim>  
inline
std::pair<unsigned int,unsigned int>
FiniteElement<dim>::component_to_base (const unsigned int index) const
{
  Assert(index < component_to_base_table.size(),
	 ExcIndexRange(index, 0, component_to_base_table.size()));

  return component_to_base_table[index];
}


template <int dim>
inline
bool
FiniteElement<dim>::restriction_is_additive (const unsigned int index) const
{
  Assert(index < this->dofs_per_cell,
	 ExcIndexRange(index, 0, this->dofs_per_cell));
  return restriction_is_additive_flags[index];
}


template <int dim>
inline
const std::vector<bool> &
FiniteElement<dim>::get_nonzero_components (const unsigned int i) const
{
  Assert (i < this->dofs_per_cell, ExcIndexRange (i, 0, this->dofs_per_cell));
  return nonzero_components[i];
}



template <int dim>
inline
unsigned int
FiniteElement<dim>::n_nonzero_components (const unsigned int i) const
{
  Assert (i < this->dofs_per_cell, ExcIndexRange (i, 0, this->dofs_per_cell));
  return n_nonzero_components_table[i];
}



template <int dim>
inline
bool
FiniteElement<dim>::is_primitive (const unsigned int i) const
{
  Assert (i < this->dofs_per_cell, ExcIndexRange (i, 0, this->dofs_per_cell));

				   // return primitivity of a shape
				   // function by checking whether it
				   // has more than one non-zero
				   // component or not. we could cache
				   // this value in an array of bools,
				   // but accessing a bit-vector (as
				   // std::vector<bool> is) is
				   // probably more expensive than
				   // just comparing against 1
  return (n_nonzero_components_table[i] == 1);
}


template <int dim>
inline
bool
FiniteElement<dim>::is_primitive () const
{
  return cached_primitivity;
}


#endif
