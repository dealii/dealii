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
 * This class gives a unified framework for implementing a
 * FiniteElement class based on a polynomial space like a
 * @p TensorProductSpace or a @p PolynomialSpace.
 *
 * Every class that implements following functions can be used as
 * template parameter @p POLY. Example classes are @p
 * TensorProductSpace and @p PolynomialSpace.
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
 *
 * This class is not a fully implemented FiniteElement class. Instead
 * there are several pure virtual functions declared in the
 * FiniteElement and FiniteElementBase classes which cannot
 * implemented by this class but are left for implementation in
 * derived classes.
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
    FE_Poly (const POLY& poly_space,
	     const FiniteElementData<dim> &fe_data,
	     const std::vector<bool> &restriction_is_additive_flags,
	     const std::vector<std::vector<bool> > &nonzero_components);

				     /**
				      * Return the value of the
				      * @p{i}th shape function at the
				      * point @p{p}. See the
				      * @ref{FiniteElementBase} base
				      * class for more information
				      * about the semantics of this
				      * function.
				      */
    virtual double shape_value (const unsigned int i,
			        const Point<dim> &p) const;
    
				     /**
				      * Return the value of the
				      * @p{component}th vector
				      * component of the @p{i}th shape
				      * function at the point
				      * @p{p}. See the
				      * @ref{FiniteElementBase} base
				      * class for more information
				      * about the semantics of this
				      * function.
				      *
				      * Since this element is scalar,
				      * the returned value is the same
				      * as if the function without the
				      * @p{_component} suffix were
				      * called, provided that the
				      * specified component is zero.
				      */
    virtual double shape_value_component (const unsigned int i,
					  const Point<dim> &p,
					  const unsigned int component) const;

				     /**
				      * Return the gradient of the
				      * @p{i}th shape function at the
				      * point @p{p}. See the
				      * @ref{FiniteElementBase} base
				      * class for more information
				      * about the semantics of this
				      * function.
				      */
    virtual Tensor<1,dim> shape_grad (const unsigned int  i,
				      const Point<dim>   &p) const;

				     /**
				      * Return the gradient of the
				      * @p{component}th vector
				      * component of the @p{i}th shape
				      * function at the point
				      * @p{p}. See the
				      * @ref{FiniteElementBase} base
				      * class for more information
				      * about the semantics of this
				      * function.
				      *
				      * Since this element is scalar,
				      * the returned value is the same
				      * as if the function without the
				      * @p{_component} suffix were
				      * called, provided that the
				      * specified component is zero.
				      */
    virtual Tensor<1,dim> shape_grad_component (const unsigned int i,
						const Point<dim> &p,
						const unsigned int component) const;

				     /**
				      * Return the tensor of second
				      * derivatives of the @p{i}th
				      * shape function at point @p{p}
				      * on the unit cell. See the
				      * @ref{FiniteElementBase} base
				      * class for more information
				      * about the semantics of this
				      * function.
				      */
    virtual Tensor<2,dim> shape_grad_grad (const unsigned int  i,
					   const Point<dim> &p) const;

				     /**
				      * Return the second derivative
				      * of the @p{component}th vector
				      * component of the @p{i}th shape
				      * function at the point
				      * @p{p}. See the
				      * @ref{FiniteElementBase} base
				      * class for more information
				      * about the semantics of this
				      * function.
				      *
				      * Since this element is scalar,
				      * the returned value is the same
				      * as if the function without the
				      * @p{_component} suffix were
				      * called, provided that the
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
				      * scalar, @p{base_element(0)} is
				      * @p{this}, and all other
				      * indices throw an error.
				      */
    virtual const FiniteElement<dim> &
    base_element (const unsigned int index) const;

                                     /**
                                      * Multiplicity of base element
                                      * @p{index}. Since this is a
                                      * scalar element,
                                      * @p{element_multiplicity(0)}
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
				      * function in
				      * @ref{FiniteElement}.
				      */
    virtual void
    fill_fe_values (const Mapping<dim> &mapping,
		    const typename DoFHandler<dim>::cell_iterator &cell,
		    const Quadrature<dim>                &quadrature,
		    typename Mapping<dim>::InternalDataBase      &mapping_internal,
		    typename Mapping<dim>::InternalDataBase      &fe_internal,
		    FEValuesData<dim>& data) const;
    
				     /**
				      * Implementation of the same
				      * function in
				      * @ref{FiniteElement}.
				      */
    virtual void
    fill_fe_face_values (const Mapping<dim> &mapping,
			 const typename DoFHandler<dim>::cell_iterator &cell,
			 const unsigned int                    face_no,
			 const Quadrature<dim-1>                &quadrature,
			 typename Mapping<dim>::InternalDataBase      &mapping_internal,
			 typename Mapping<dim>::InternalDataBase      &fe_internal,
			 FEValuesData<dim>& data) const ;
    
				     /**
				      * Implementation of the same
				      * function in
				      * @ref{FiniteElement}.
				      */
    virtual void
    fill_fe_subface_values (const Mapping<dim> &mapping,
			    const typename DoFHandler<dim>::cell_iterator &cell,
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
				      * values required by @p{flags}.
				      *
				      * For the purpuse of this
				      * function, refer to the
				      * documentation in
				      * @p{FiniteElement}.
				      *
				      * The effect in this element is
				      * as follows: if
				      * @p{update_values} is set in
				      * @p{flags}, copy it to the
				      * result. All other flags of the
				      * result are cleared, since
				      * everything else must be
				      * computed for each cell.
				      */
    virtual UpdateFlags update_once (const UpdateFlags flags) const;
  
				     /**
				      * Determine the values that need
				      * to be computed on every
				      * cell to be able to compute all
				      * values required by @p{flags}.
				      *
				      * For the purpuse of this
				      * function, refer to the
				      * documentation in
				      * @p{FiniteElement}.
				      *
				      * The effect in this element is
				      * as follows:
				      * @begin{itemize}
				      * @item if @p{update_gradients}
				      * is set, the result will
				      * contain @p{update_gradients}
				      * and
				      * @p{update_covariant_transformation}.
				      * The latter is required to
				      * transform the gradient on the
				      * unit cell to the real
				      * cell. Remark, that the action
				      * required by
				      * @p{update_covariant_transformation}
				      * is actually performed by the
				      * @p{Mapping} object used in
				      * conjunction with this finite
				      * element.
				      * @item if
				      * @p{update_second_derivatives}
				      * is set, the result will
				      * contain
				      * @p{update_second_derivatives}
				      * and
				      * @p{update_covariant_transformation}.
				      * The rationale is the same as
				      * above and no higher
				      * derivatives of the
				      * transformation are required,
				      * since we use difference
				      * quotients for the actual
				      * computation.
				      * @end{itemize}
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
	Table<2,double> shape_values;

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
	Table<2,Tensor<1,dim> > shape_gradients;
    };

                                     /**
                                      * The polynomial space.
                                      */    
    POLY poly_space;
};



#endif
