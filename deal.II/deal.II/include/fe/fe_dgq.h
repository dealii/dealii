//---------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 1998, 1999, 2000, 2001 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------
#ifndef __deal2__fe_dgq_h
#define __deal2__fe_dgq_h

#include <base/polynomial.h>
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
    FE_DGQ (unsigned int k);
				     /**
				      * Destructor.
				      */
    ~FE_DGQ ();
    
				     /**
				      * Return the value of the
				      * @p{i}th shape function at the
				      * point @p{p}.  @p{p} is a point
				      * on the reference element.
				      */
    virtual double shape_value (const unsigned int i,
			        const Point<dim> &p) const;
    
				     /**
				      * Return the gradient of the
				      * @p{i}th shape function at the
				      * point @p{p}. @p{p} is a point
				      * on the reference element, and
				      * likewise the gradient is the
				      * gradient on the unit cell with
				      * respect to unit cell
				      * coordinates.
				      */
    virtual Tensor<1,dim> shape_grad (const unsigned int  i,
				      const Point<dim>   &p) const;

				     /**
				      * Return the tensor of second
				      * derivatives of the @p{i}th
				      * shape function at point @p{p}
				      * on the unit cell. The
				      * derivatives are derivatives on
				      * the unit cell with respect to
				      * unit cell coordinates.
				      */
    virtual Tensor<2,dim> shape_grad_grad (const unsigned int  i,
					   const Point<dim> &p) const;
    
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
				      * Declare a nested class which
				      * will has static definitions of
				      * various matrices such as
				      * constraint and embedding
				      * matrices. The definition of
				      * the various static fields are
				      * in the files @p{fe_q_[123]d.cc}
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
				      * Compute flags for initial update only.
				      */
    virtual UpdateFlags update_once (const UpdateFlags flags) const;
  
				     /**
				      * Compute flags for update on each cell.
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
    TensorProductPolynomials<dim>* poly;

				     /**
				      * Fields of cell-independent data.
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
    friend class MappingQ<dim>;
};


// declaration of explicit specializations

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
