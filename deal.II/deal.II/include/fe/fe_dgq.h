//---------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 1998, 1999, 2000, 2001, 2002 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------
#ifndef __deal2__fe_dgq_h
#define __deal2__fe_dgq_h

#include <base/config.h>
#include <base/polynomial.h>
#include <base/tensor_product_polynomials.h>
#include <fe/fe.h>

template <int dim> class TensorProductPolynomials;
template <int dim> class MappingQ;


/**
 * Discontinuous tensor product elements based on equidistant support points.
 *
 * This is a discontinuous finite element using interpolating tensor
 * product polynomials. The shape functions are Lagrangian
 * interpolants of an equidistant grid of points on the unit cell. The
 * points are numbered in lexicographical order, @p{x} running fastest.
 *
 * @author Guido Kanschat, Ralf Hartmann, 2001
 */
template <int dim>
class FE_DGQ : public FiniteElement<dim>
{
  public:
				     /**
				      * Constructor for tensor product
				      * polynomials of degree @p{k}.
				      */
    FE_DGQ (const unsigned int k);
    
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
				      * Return the polynomial degree
				      * of this finite element,
				      * i.e. the value passed to the
				      * constructor.
				      */
    unsigned int get_degree () const;

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
    virtual const FiniteElement<dim> & base_element (const unsigned int index) const;
    
				     /**
				      * Check for non-zero values on a face.
				      *
				      * This function returns
				      * @p{true}, if the shape
				      * function @p{shape_index} has
				      * non-zero values on the face
				      * @p{face_index}.
				      *
				      * Implementation of the
				      * interface in
				      * @ref{FiniteElement}
				      */
    virtual bool has_support_on_face (const unsigned int shape_index,
				      const unsigned int face_index) const;

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
				      * Declare a nested class which
				      * will hold static definitions of
				      * various matrices such as
				      * constraint and embedding
				      * matrices. The definition of
				      * the various static fields are
				      * in the files @p{fe_dgq_[123]d.cc}
				      * in the source directory.
				      */
    struct Matrices
    {
					 /**
					  * Pointers to the embedding
					  * matrices, one for each
					  * polynomial degree starting
					  * from constant elements
					  */
	static const double * const embedding[];

					 /**
					  * Number of elements (first
					  * index) the above field
					  * has. Equals the highest
					  * polynomial degree plus one
					  * for which the embedding
					  * matrices have been
					  * computed.
					  */
	static const unsigned int n_embedding_matrices;

					 /**
					  * As @p{embedding} but for
					  * projection matrices.
					  */
	static const double * const projection_matrices[];

					 /**
					  * As
					  * @p{n_embedding_matrices}
					  * but for projection
					  * matrices.
					  */
	static const unsigned int n_projection_matrices;
    };

    
  protected:

				     /**
				      * @p{clone} function instead of
				      * a copy constructor.
				      *
				      * This function is needed by the
				      * constructors of @p{FESystem}.
				      */
    virtual FiniteElement<dim> *clone() const;
  
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

  private:
				     /**
				      * Only for internal use. Its
				      * full name is
				      * @p{get_dofs_per_object_vector}
				      * function and it creates the
				      * @p{dofs_per_object} vector that is
				      * needed within the constructor to
				      * be passed to the constructor of
				      * @p{FiniteElementData}.
				      */
    static std::vector<unsigned int> get_dpo_vector(unsigned int degree);
    
				     /**
				      * Given a set of flags indicating
				      * what quantities are requested
				      * from a @p{FEValues} object,
				      * return which of these can be
				      * precomputed once and for
				      * all. Often, the values of
				      * shape function at quadrature
				      * points can be precomputed, for
				      * example, in which case the
				      * return value of this function
				      * would be the logical and of
				      * the input @p{flags} and
				      * @p{update_values}.
				      *
				      * For the present kind of finite
				      * element, this is exactly the
				      * case.
				      */
    virtual UpdateFlags update_once (const UpdateFlags flags) const;
  
				     /**
				      * This is the opposite to the
				      * above function: given a set of
				      * flags indicating what we want
				      * to know, return which of these
				      * need to be computed each time
				      * we visit a new cell.
				      *
				      * If for the computation of one
				      * quantity something else is
				      * also required (for example, we
				      * often need the covariant
				      * transformation when gradients
				      * need to be computed), include
				      * this in the result as well.
				      */
    virtual UpdateFlags update_each (const UpdateFlags flags) const;
  
				     /**
				      * Compute renumbering for rotation
				      * of degrees of freedom.
				      *
				      * Rotates a tensor product
				      * numbering of degrees of
				      * freedom by 90 degrees. It is
				      * used to compute the transfer
				      * matrices of the children by
				      * using only the matrix for the
				      * first child.
				      *
				      * The direction parameter
				      * determines the type of
				      * rotation. It is one character
				      * of @p{xXyYzZ}. The character
				      * determines the axis of
				      * rotation, case determines the
				      * direction. Lower case is
				      * counter-clockwise seen in
				      * direction of the axis.
				      *
				      * Since rotation around the
				      * y-axis is not used, it is not
				      * implemented either.
				      */
    void rotate_indices (std::vector<unsigned int> &indices,
			 const char                 direction) const;
  
				     /**
				      * Degree of the polynomials.
				      */  
    const unsigned int degree;

				     /**
				      * Pointer to the tensor
				      * product polynomials.
				      */
    const TensorProductPolynomials<dim> polynomial_space;

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
					  * vector for each shape
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
					  * vector for each shape
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
	typename std::vector<std::vector<Tensor<1,dim> > > shape_gradients;
    };
    
				     /**
				      * Allow access from other dimensions.
				      */
    template <int dim1> friend class FE_DGQ;

				     /**
				      * Allows @p{MappingQ} class to
				      * access to build_renumbering
				      * function.
				      */
    template <int dim1> friend class MappingQ;
};


// declaration of explicit specializations of member variables, if the
// compiler allows us to do that (the standard says we must)
#ifndef DEAL_II_MEMBER_VAR_SPECIALIZATION_BUG
template <> 
const double * const FE_DGQ<1>::Matrices::embedding[];

template <>
const unsigned int FE_DGQ<1>::Matrices::n_embedding_matrices;

template <>
const double * const FE_DGQ<1>::Matrices::projection_matrices[];

template <>
const unsigned int FE_DGQ<1>::Matrices::n_projection_matrices;

template <> 
const double * const FE_DGQ<2>::Matrices::embedding[];

template <>
const unsigned int FE_DGQ<2>::Matrices::n_embedding_matrices;

template <>
const double * const FE_DGQ<2>::Matrices::projection_matrices[];

template <>
const unsigned int FE_DGQ<2>::Matrices::n_projection_matrices;

template <> 
const double * const FE_DGQ<3>::Matrices::embedding[];

template <>
const unsigned int FE_DGQ<3>::Matrices::n_embedding_matrices;

template <>
const double * const FE_DGQ<3>::Matrices::projection_matrices[];

template <>
const unsigned int FE_DGQ<3>::Matrices::n_projection_matrices;
#endif

#endif
