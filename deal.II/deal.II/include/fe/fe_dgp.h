//----------------------------------------------------------------
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
#ifndef __deal2__fe_dgp_h
#define __deal2__fe_dgp_h

#include <base/config.h>
#include <base/polynomial_space.h>
#include <fe/fe_poly.h>

template <int dim> class MappingQ;

/*!@addtogroup fe */
/*@{*/

/**
 * Discontinuous finite elements based on Legendre polynomials.
 *
 * This finite element implements complete polynomial spaces, that is,
 * dim-dimensional polynomials of degree p. The underlying
 * polynomials form a Legendre basis on the unit square. Thus, the
 * mass matrix is diagonal, if the grid cells are parallelograms.
 *
 * @author Guido Kanschat, 2001, 2002, Ralf Hartmann 2004
 */
template <int dim>
class FE_DGP : public FE_Poly<PolynomialSpace<dim>,dim>
{
  public:
				     /**
				      * Constructor for tensor product
				      * polynomials of degree @p p.
				      */
    FE_DGP (const unsigned int p);
    
				     /**
				      * Return a string that uniquely
				      * identifies a finite
				      * element. This class returns
				      * @p{FE_DGP<dim>(degree)}, with
				      * @p dim and @p degree
				      * replaced by appropriate
				      * values.
				      */
    virtual std::string get_name () const;
    
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
				      * in the files @p{fe_dgp_[123]d.cc}
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
	static const double * const embedding[][GeometryInfo<dim>::children_per_cell];

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
	static const double * const projection_matrices[][GeometryInfo<dim>::children_per_cell];

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
				      * Allow access from other dimensions.
				      */
    template <int dim1> friend class FE_DGP;

				     /**
				      * Allows @p MappingQ class
				      * access to build_renumbering
				      * function.
				      */
    friend class MappingQ<dim>;
};

/* @} */
/// @if NoDoc


// declaration of explicit specializations of member variables, if the
// compiler allows us to do that (the standard says we must)
#ifndef DEAL_II_MEMBER_VAR_SPECIALIZATION_BUG
template <> 
const double * const FE_DGP<1>::Matrices::embedding[][GeometryInfo<1>::children_per_cell];

template <>
const unsigned int FE_DGP<1>::Matrices::n_embedding_matrices;

template <>
const double * const FE_DGP<1>::Matrices::projection_matrices[][GeometryInfo<1>::children_per_cell];

template <>
const unsigned int FE_DGP<1>::Matrices::n_projection_matrices;

template <> 
const double * const FE_DGP<2>::Matrices::embedding[][GeometryInfo<2>::children_per_cell];

template <>
const unsigned int FE_DGP<2>::Matrices::n_embedding_matrices;

template <>
const double * const FE_DGP<2>::Matrices::projection_matrices[][GeometryInfo<2>::children_per_cell];

template <>
const unsigned int FE_DGP<2>::Matrices::n_projection_matrices;

template <> 
const double * const FE_DGP<3>::Matrices::embedding[][GeometryInfo<3>::children_per_cell];

template <>
const unsigned int FE_DGP<3>::Matrices::n_embedding_matrices;

template <>
const double * const FE_DGP<3>::Matrices::projection_matrices[][GeometryInfo<3>::children_per_cell];

template <>
const unsigned int FE_DGP<3>::Matrices::n_projection_matrices;
#endif

/// @endif
#endif
