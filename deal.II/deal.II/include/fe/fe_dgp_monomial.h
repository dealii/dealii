//---------------------------------------------------------------------------
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
//---------------------------------------------------------------------------
#ifndef __deal2__fe_dgp_monomial_h
#define __deal2__fe_dgp_monomial_h

#include <base/config.h>
#include <base/polynomials_p.h>
#include <fe/fe_poly.h>

template <int dim> class MappingQ;


/*!@addtogroup fe */
/*@{*/

/**
 * Discontinuous finite elements based on monomials of degree up to
 * <tt>p</tt>.
 *
 * This finite element makes use of the PolynomialsP class which
 * implements <tt>dim</tt>-dimensional polynomials of degree
 * <tt>p</tt> based the Polynomials::Polynomial and the
 * PolynomialSpace classes.
 *
 * @author Ralf Hartmann, 2004
 */
template <int dim>
class FE_DGPMonomial : public FE_Poly<PolynomialsP<dim>,dim>
{
  public:
				     /**
				      * Constructor for the polynomial
				      * space of degree <tt>p</tt>.
				      */
    FE_DGPMonomial (const unsigned int p);
    
				     /**
				      * Return a string that uniquely
				      * identifies a finite
				      * element. This class returns
				      * <tt>FE_DGPMonomial<dim>(degree)</tt>,
				      * with <tt>dim</tt> and
				      * <tt>p</tt> replaced by
				      * appropriate values.
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
				      * element is also a @p FE_Q
				      * element. Otherwise, an
				      * exception of type
				      * FiniteElementBase<dim>::ExcInterpolationNotImplemented
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
				      * Initialize the embedding
				      * matrices. Called from the
				      * constructor.
				      */
    void initialize_embedding ();

				     /**
				      * Initialize the restriction
				      * matrices. Called from the
				      * constructor.
				      */
    void initialize_restriction ();

				     /**
				      * Allows @p MappingQ class
				      * access to build_renumbering
				      * function.
				      */
    friend class MappingQ<dim>;
};

/*@}*/


#endif
