//---------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 1998, 1999, 2000, 2001, 2002, 2003, 2004 by the deal.II authors
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
#include <base/tensor_product_polynomials.h>
#include <fe/fe_poly.h>

template <int dim> class MappingQ;

/*!@addtogroup fe */
/*@{*/

/**
 * Discontinuous tensor product elements based on equidistant support
 * points.
 *
 * This is a discontinuous finite element based on tensor products of
 * Lagrangian polynomials. The shape functions are Lagrangian
 * interpolants of an equidistant grid of points on the unit cell. The
 * points are numbered in lexicographical order, with @p x running
 * fastest.
 *
 * @author Ralf Hartmann, 2001, 2004, Guido Kanschat 2001
 */
template <int dim>
class FE_DGQ : public FE_Poly<TensorProductPolynomials<dim>,dim>
{
  public:
				     /**
				      * Constructor for tensor product
				      * polynomials of degree @p p.
				      */
    FE_DGQ (const unsigned int p);
    
				     /**
				      * Return a string that uniquely
				      * identifies a finite
				      * element. This class returns
				      * @p FE_DGQ<dim>(degree) , with
				      * @p dim and @p degree
				      * replaced by appropriate
				      * values.
				      */
    virtual std::string get_name () const;
    
				     /**
				      * Return the matrix
				      * interpolating from the given
				      * finite element to the present
				      * one. The size of the matrix is
				      * then @p dofs_per_cell times
				      * <tt>source.dofs_per_cell</tt>.
				      *
				      * These matrices are only
				      * available if the source
				      * element is also a @p FE_DGQ
				      * element. Otherwise, an
				      * exception of type
				      * @ref{FiniteElementBase<dim>::ExcInterpolationNotImplemented}
				      * is thrown.
				      */
    virtual void
    get_interpolation_matrix (const FiniteElementBase<dim> &source,
			      FullMatrix<double>           &matrix) const;
    
				     /**
				      * Check for non-zero values on a face.
				      *
				      * This function returns
				      * @p true, if the shape
				      * function @p shape_index has
				      * non-zero values on the face
				      * @p face_index.
				      *
				      * Implementation of the
				      * interface in
				      * FiniteElement
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
				      * in the files <tt>fe_dgq_[123]d.cc</tt>
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
					  * As @p embedding but for
					  * projection matrices.
					  */
	static const double * const projection_matrices[];

					 /**
					  * As
					  * @p n_embedding_matrices
					  * but for projection
					  * matrices.
					  */
	static const unsigned int n_projection_matrices;
    };

    
  protected:

				     /**
				      * @p clone function instead of
				      * a copy constructor.
				      *
				      * This function is needed by the
				      * constructors of @p FESystem.
				      */
    virtual FiniteElement<dim> *clone() const;

  private:
				     /**
				      * Only for internal use. Its
				      * full name is
				      * @p get_dofs_per_object_vector
				      * function and it creates the
				      * @p dofs_per_object vector that is
				      * needed within the constructor to
				      * be passed to the constructor of
				      * @p FiniteElementData.
				      */
    static std::vector<unsigned int> get_dpo_vector(unsigned int degree);
  
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
				      * of @p xXyYzZ. The character
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
				      * Allow access from other dimensions.
				      */
    template <int dim1> friend class FE_DGQ;

				     /**
				      * Allows @p MappingQ class to
				      * access to build_renumbering
				      * function.
				      */
    template <int dim1> friend class MappingQ;
};

/*@}*/

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
