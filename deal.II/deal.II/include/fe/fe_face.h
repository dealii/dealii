//---------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2009 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------
#ifndef __deal2__fe_face_h
#define __deal2__fe_face_h

#include <base/config.h>
#include <base/tensor_product_polynomials.h>
#include <fe/fe_poly_face.h>

DEAL_II_NAMESPACE_OPEN


/**
 * @warning This class has not been tested
 *
 * A finite element, which is a tensor product polynomial on each face
 * and undefined in the interior of the cells.
 *
 * This finite element is the trace space of FE_RaviartThomas on the
 * faces and serves in hybridized methods.
 *
 * @author Guido Kanschat, 2009
 */
template <int dim, int spacedim=dim>
class FE_FaceQ : public FE_PolyFace<TensorProductPolynomials<dim-1>, dim, spacedim>
{
  public:
				     /**
				      * Constructor for tensor product
				      * polynomials of degree
				      * <tt>p</tt>. The shape
				      * functions created using this
				      * constructor correspond to
				      * Legendre polynomials in each
				      * coordinate direction.
				      */
    FE_FaceQ(unsigned int p);
    
    virtual FiniteElement<dim,spacedim>* clone() const;

				     /**
				      * Return a string that uniquely
				      * identifies a finite
				      * element. This class returns
				      * <tt>FE_DGQ<dim>(degree)</tt> , with
				      * <tt>dim</tt> and <tt>degree</tt>
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

  private:
				     /**
				      * Return vector with dofs per
				      * vertex, line, quad, hex.
				      */
    static std::vector<unsigned int> get_dpo_vector (const unsigned int deg);
};

DEAL_II_NAMESPACE_CLOSE

#endif
