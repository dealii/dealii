//---------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2002, 2003, 2004, 2005, 2006 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------
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
 * dim-dimensional polynomials of degree p. For example, in 2d the
 * element FE_DGP(1) would represent the span of the functions
 * $\{1,\hat x,\hat y\}$, which is in contrast to the element FE_DGQ(1)
 * that is formed by the span of $\{1,\hat x,\hat y,\hat x\hat y\}$. Since the
 * DGP space has only three unknowns for each quadrilateral, it is
 * immediately clear that this element can not be continuous.
 *
 * The basis functions for this element are chosen to form a Legendre
 * basis on the unit square. Thus, the mass matrix is diagonal, if the
 * grid cells are parallelograms. Note that this is in contrast to the
 * FE_DGPMonomial class that actually uses the monomial basis listed
 * above as basis functions.
 *
 * The shape functions are defined in the class PolynomialSpace. The
 * polynomials used inside PolynomialSpace are Polynomials::Legendre
 * up to degree <tt>p</tt> given in FE_DGP. For the ordering of the
 * basis functions, refer to PolynomialSpace, remebering that the
 * Legendre polynomials are ordered by ascending degree.
 *
 *
 * <h3>Transformation properties</h3>
 *
 * It is worth noting that under a (bi-, tri-)linear mapping, the
 * space described by this element does not contain $P(k)$, even if we
 * use a basis of polynomials of degree $k$. Consequently, for
 * example, on meshes with non-affine cells, a linear function can not
 * be exactly represented by elements of type FE_DGP(1) or
 * FE_DGPMonomial(1).
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
				      * <tt>FE_DGP<dim>(degree)</tt>, with
				      * @p dim and @p degree
				      * replaced by appropriate
				      * values.
				      */
    virtual std::string get_name () const;
    
                                     /**
                                      * Return whether this element
                                      * implements its hanging node
                                      * constraints in the new way,
				      * which has to be used to make
				      * elements "hp compatible".
                                      *
				      * For the FE_DGP class the result is
				      * always true (independent of the degree
				      * of the element), as it has no hanging
				      * nodes (being a discontinuous element).
                                      */
    virtual bool hp_constraints_are_implemented () const;

				     /**
				      * Return the matrix
				      * interpolating from a face of
				      * of one element to the face of
				      * the neighboring element. 
				      * The size of the matrix is
				      * then @p dofs_per_face times
				      * <tt>source.dofs_per_face</tt>.
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
    get_face_interpolation_matrix (const FiniteElement<dim> &source,
				   FullMatrix<double>       &matrix) const;    

				     /**
				      * Return the matrix
				      * interpolating from a face of
				      * of one element to the face of
				      * the neighboring element. 
				      * The size of the matrix is
				      * then @p dofs_per_face times
				      * <tt>source.dofs_per_face</tt>.
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
    get_subface_interpolation_matrix (const FiniteElement<dim> &source,
				      const unsigned int        subface,
				      FullMatrix<double>       &matrix) const;

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
				      * in the files <tt>fe_dgp_[123]d.cc</tt>
				      * in the source directory.
				      */
    struct Matrices
    {
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
#ifndef DOXYGEN


// declaration of explicit specializations of member variables, if the
// compiler allows us to do that (the standard says we must)
#ifndef DEAL_II_MEMBER_VAR_SPECIALIZATION_BUG
template <>
const double * const FE_DGP<1>::Matrices::projection_matrices[][GeometryInfo<1>::children_per_cell];

template <>
const unsigned int FE_DGP<1>::Matrices::n_projection_matrices;

template <>
const double * const FE_DGP<2>::Matrices::projection_matrices[][GeometryInfo<2>::children_per_cell];

template <>
const unsigned int FE_DGP<2>::Matrices::n_projection_matrices;

template <>
const double * const FE_DGP<3>::Matrices::projection_matrices[][GeometryInfo<3>::children_per_cell];

template <>
const unsigned int FE_DGP<3>::Matrices::n_projection_matrices;
#endif

#endif // DOXYGEN
#endif
