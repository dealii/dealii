// ---------------------------------------------------------------------
//
// Copyright (C) 2000 - 2014 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------

#ifndef __deal2__fe_trace_h
#define __deal2__fe_trace_h

#include <deal.II/base/config.h>
#include <deal.II/base/tensor_product_polynomials.h>
#include <deal.II/fe/fe_poly_face.h>

DEAL_II_NAMESPACE_OPEN

/**
 * A finite element, which is the trace of FE_Q elements, that is
 * a tensor product of polynomials on the faces,
 * undefined in the interior of the cells and continuous. The basis functions on
 * the faces are from Polynomials::LagrangeEquidistant
 *
 * This finite element is the trace space of FE_Q on the
 * faces.
 *
 * @note Since these are only finite elements on faces, only
 * FEFaceValues and FESubfaceValues will be able to extract reasonable
 * values from any face polynomial. In order to make the use of
 * FESystem simpler, FEValues objects will not fail using this finite
 * element space, but all shape function values extracted will equal
 * to zero.
 *
 * @todo Polynomials::LagrangeEquidistant should be and will be
 * replaced by Polynomials::LagrangeGaussLobatto as soon as such a
 * polynomial set exists.
 * @todo so far, hanging nodes are not implemented
 *
 */

template <int dim, int spacedim=dim>
class FE_TraceQ : public FE_PolyFace<TensorProductPolynomials<dim-1>, dim, spacedim>
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
  FE_TraceQ(unsigned int p);

  virtual FiniteElement<dim,spacedim> *clone() const;

  /**
   * Return a string that uniquely
   * identifies a finite
   * element. This class returns
   * <tt>FE_DGQ<dim>(degree)</tt>, with
   * <tt>dim</tt> and <tt>degree</tt>
   * replaced by appropriate
   * values.
   */
  virtual std::string get_name () const;

  /**
   * This function returns @p true, if the shape function @p shape_index has
   * non-zero function values somewhere on the face @p face_index.
   */
  virtual bool has_support_on_face (const unsigned int shape_index,
                                    const unsigned int face_index) const;

  /**
   * Returns a list of constant modes of the element. For this element, it
   * simply returns one row with all entries set to true.
   */
  virtual std::pair<Table<2,bool>, std::vector<unsigned int> >
  get_constant_modes () const;

private:
  /**
   * Return vector with dofs per
   * vertex, line, quad, hex.
   */
  static std::vector<unsigned int> get_dpo_vector (const unsigned int deg);
};

DEAL_II_NAMESPACE_CLOSE

#endif
