//---------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2004 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------
#ifndef __deal2__fe_poly_h
#define __deal2__fe_poly_h


#include <fe/fe.h>


/**
 * This class gives a unified framework for the implementation of
 * FiniteElement classes based on a polynomial spaces like the
 * TensorProductPolynomials or a PolynomialSpace classes.
 *
 * Every class that implements following functions can be used as
 * template parameter POLY.
 *
 * @code
 * double compute_value (const unsigned int i,
 *                       const Point<dim> &p) const;
 *
 * Tensor<1,dim> compute_grad (const unsigned int i,
 *                             const Point<dim> &p) const;
 *
 * Tensor<2,dim> compute_grad_grad (const unsigned int i,
 *                                  const Point<dim> &p) const;
 * @endcode
 * Example classes are TensorProductPolynomials, PolynomialSpace or
 * PolynomialsP.
 *
 * This class is not a fully implemented FiniteElement class. Instead
 * there are several pure virtual functions declared in the
 * FiniteElement and FiniteElementBase classes which cannot
 * implemented by this class but are left for implementation in
 * derived classes.
 *
 * Furthermore, this class assumes that shape functions of the
 * FiniteElement under consideration do <em>not</em> depend on the
 * actual shape of the cells in real space, i.e. update_once()
 * includes update_values. For FiniteElements whose shape functions
 * depend on the cells in real space, the update_once() and
 * update_each() functions must be overloaded.
 *
 * Todos:
 * - checke dim of POLY
 * - templatisiere nur auf POLY, dim ergibt sich durch POLY::dim
 *
 * @author Ralf Hartmann 2004
 **/
template <class POLY, int dim/*=POLY::dim*/>
class FE_Poly : public FiniteElement<dim>
{
  public:
				     /**
				      * Constructor.
				      */
    FE_Poly (unsigned int degree,
	     const POLY& poly_space,
	     const FiniteElementData<dim> &fe_data,
	     const std::vector<bool> &restriction_is_additive_flags,
	     const std::vector<std::vector<bool> > &nonzero_components);

				     /**
				      * Return the polynomial degree
				      * of this finite element,
				      * i.e. the value passed to the
				      * constructor.
				      */
    unsigned int get_degree () const;

				     /**
				      * Return the value of the
				      * <tt>i</tt>th shape function at
				      * the point <tt>p</tt>. See the
				      * FiniteElementBase base class
				      * for more information about the
				      * semantics of this function.
				      */
    virtual double shape_value (const unsigned int i,
			        const Point<dim> &p) const;
    
				     /**
				      * Return the value of the
				      * <tt>component</tt>th vector
				      * component of the <tt>i</tt>th
				      * shape function at the point
				      * <tt>p</tt>. See the
				      * FiniteElementBase base class
				      * for more information about the
				      * semantics of this function.
				      *
				      * Since this element is scalar,
				      * the returned value is the same
				      * as if the function without the
				      * <tt>_component</tt> suffix
				      * were called, provided that the
				      * specified component is zero.
				      */
    virtual double shape_value_component (const unsigned int i,
					  const Point<dim> &p,
					  const unsigned int component) const;

				     /**
				      * Return the gradient of the
				      * <tt>i</tt>th shape function at
				      * the point <tt>p</tt>. See the
				      * FiniteElementBase base class
				      * for more information about the
				      * semantics of this function.
				      */
    virtual Tensor<1,dim> shape_grad (const unsigned int  i,
				      const Point<dim>   &p) const;

				     /**
				      * Return the gradient of the
				      * <tt>component</tt>th vector
				      * component of the <tt>i</tt>th
				      * shape function at the point
				      * <tt>p</tt>. See the
				      * FiniteElementBase base class
				      * for more information about the
				      * semantics of this function.
				      *
				      * Since this element is scalar,
				      * the returned value is the same
				      * as if the function without the
				      * <tt>_component</tt> suffix
				      * were called, provided that the
				      * specified component is zero.
				      */
    virtual Tensor<1,dim> shape_grad_component (const unsigned int i,
						const Point<dim> &p,
						const unsigned int component) const;

				     /**
				      * Return the tensor of second
				      * derivatives of the
				      * <tt>i</tt>th shape function at
				      * point <tt>p</tt> on the unit
				      * cell. See the
				      * FiniteElementBase base class
				      * for more information about the
				      * semantics of this function.
				      */
    virtual Tensor<2,dim> shape_grad_grad (const unsigned int  i,
					   const Point<dim> &p) const;

				     /**
				      * Return the second derivative
				      * of the <tt>component</tt>th
				      * vector component of the
				      * <tt>i</tt>th shape function at
				      * the point <tt>p</tt>. See the
				      * FiniteElementBase base class
				      * for more information about the
				      * semantics of this function.
				      *
				      * Since this element is scalar,
				      * the returned value is the same
				      * as if the function without the
				      * <tt>_component</tt> suffix
				      * were called, provided that the
				      * specified component is zero.
				      */
    virtual Tensor<2,dim> shape_grad_grad_component (const unsigned int i,
						     const Point<dim> &p,
						     const unsigned int component) const;

                                     /**
				      * Number of base elements in a
				      * mixed discretization. Since
				      * this is a scalar element,
				      * return one.
				      */
    virtual unsigned int n_base_elements () const;
    
				     /**
				      * Access to base element
				      * objects. Since this element is
				      * scalar,
				      * <tt>base_element(0)</tt> is
				      * <tt>this</tt>, and all other
				      * indices throw an error.
				      */
    virtual const FiniteElement<dim> &
    base_element (const unsigned int index) const;

                                     /**
                                      * Multiplicity of base element
                                      * <tt>index</tt>. Since this is
                                      * a scalar element,
                                      * <tt>element_multiplicity(0)</tt>
                                      * returns one, and all other
                                      * indices will throw an error.
                                      */
    virtual unsigned int element_multiplicity (const unsigned int index) const;

    
  protected:
      
				     /**
				      * Prepare internal data
				      * structures and fill in values
				      * independent of the cell.
				      */
    virtual
    typename Mapping<dim>::InternalDataBase *
    get_data (const UpdateFlags,
	      const Mapping<dim>& mapping,
	      const Quadrature<dim>& quadrature) const ;

				     /**
				      * Implementation of the same
				      * function in FiniteElement.
				      */
    virtual void
    fill_fe_values (const Mapping<dim> &mapping,
		    const typename Triangulation<dim>::cell_iterator &cell,
		    const Quadrature<dim>                &quadrature,
		    typename Mapping<dim>::InternalDataBase      &mapping_internal,
		    typename Mapping<dim>::InternalDataBase      &fe_internal,
		    FEValuesData<dim>& data) const;
    
				     /**
				      * Implementation of the same
				      * function in FiniteElement.
				      */
    virtual void
    fill_fe_face_values (const Mapping<dim> &mapping,
			 const typename Triangulation<dim>::cell_iterator &cell,
			 const unsigned int                    face_no,
			 const Quadrature<dim-1>                &quadrature,
			 typename Mapping<dim>::InternalDataBase      &mapping_internal,
			 typename Mapping<dim>::InternalDataBase      &fe_internal,
			 FEValuesData<dim>& data) const ;
    
				     /**
				      * Implementation of the same
				      * function in FiniteElement.
				      */
    virtual void
    fill_fe_subface_values (const Mapping<dim> &mapping,
			    const typename Triangulation<dim>::cell_iterator &cell,
			    const unsigned int                    face_no,
			    const unsigned int                    sub_no,
			    const Quadrature<dim-1>                &quadrature,
			    typename Mapping<dim>::InternalDataBase      &mapping_internal,
			    typename Mapping<dim>::InternalDataBase      &fe_internal,
			    FEValuesData<dim>& data) const ;

    
				     /**
				      * Determine the values that need
				      * to be computed on the unit
				      * cell to be able to compute all
				      * values required by
				      * <tt>flags</tt>.
				      *
				      * For the purpuse of this
				      * function, refer to the
				      * documentation in
				      * FiniteElement.
				      *
				      * This class assumes that shape
				      * functions of this
				      * FiniteElement do <em>not</em>
				      * depend on the actual shape of
				      * the cells in real
				      * space. Therefore, the effect
				      * in this element is as follows:
				      * if <tt>update_values</tt> is
				      * set in <tt>flags</tt>, copy it
				      * to the result. All other flags
				      * of the result are cleared,
				      * since everything else must be
				      * computed for each cell.
				      */
    virtual UpdateFlags update_once (const UpdateFlags flags) const;
  
				     /**
				      * Determine the values that need
				      * to be computed on every cell
				      * to be able to compute all
				      * values required by
				      * <tt>flags</tt>.
				      *
				      * For the purpuse of this
				      * function, refer to the
				      * documentation in
				      * FiniteElement.
				      *
				      * This class assumes that shape
				      * functions of this
				      * FiniteElement do <em>not</em>
				      * depend on the actual shape of
				      * the cells in real
				      * space.
				      *
				      * The effect in this element is
				      * as follows:
				      * <ul>

				      * <li> if
				      * <tt>update_gradients</tt> is
				      * set, the result will contain
				      * <tt>update_gradients</tt> and
				      * <tt>update_covariant_transformation</tt>.
				      * The latter is required to
				      * transform the gradient on the
				      * unit cell to the real
				      * cell. Remark, that the action
				      * required by
				      * <tt>update_covariant_transformation</tt>
				      * is actually performed by the
				      * Mapping object used in
				      * conjunction with this finite
				      * element.  <li> if
				      * <tt>update_second_derivatives</tt>
				      * is set, the result will
				      * contain
				      * <tt>update_second_derivatives</tt>
				      * and
				      * <tt>update_covariant_transformation</tt>.
				      * The rationale is the same as
				      * above and no higher
				      * derivatives of the
				      * transformation are required,
				      * since we use difference
				      * quotients for the actual
				      * computation.
				      *
				      * </ul>
				      */
    virtual UpdateFlags update_each (const UpdateFlags flags) const;


				     /**
				      * Fields of cell-independent data.
				      *
				      * For information about the
				      * general purpose of this class,
				      * see the documentation of the
				      * base class.
				      */
    class InternalData : public FiniteElementBase<dim>::InternalDataBase
    {
      public:
					 /**
					  * Array with shape function
					  * values in quadrature
					  * points. There is one
					  * row for each shape
					  * function, containing
					  * values for each quadrature
					  * point.
					  *
					  * In this array, we store
					  * the values of the shape
					  * function in the quadrature
					  * points on the unit
					  * cell. Since these values
					  * do not change under
					  * transformation to the real
					  * cell, we only need to copy
					  * them over when visiting a
					  * concrete cell.
					  */
	std::vector<std::vector<double> > shape_values;

					 /**
					  * Array with shape function
					  * gradients in quadrature
					  * points. There is one
					  * row for each shape
					  * function, containing
					  * values for each quadrature
					  * point.
					  *
					  * We store the gradients in
					  * the quadrature points on
					  * the unit cell. We then
					  * only have to apply the
					  * transformation (which is a
					  * matrix-vector
					  * multiplication) when
					  * visiting an actual cell.
					  */      
	std::vector<std::vector<Tensor<1,dim> > > shape_gradients;
    };
    
				     /**
				      * Degree of the polynomials.
				      */  
    const unsigned int degree;

                                     /**
                                      * The polynomial space. Its type
                                      * is given by the template
                                      * parameter POLY.
                                      */    
    POLY poly_space;
};



#endif
