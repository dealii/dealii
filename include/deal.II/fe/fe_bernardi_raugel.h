// ---------------------------------------------------------------------
//
// Copyright (C) 2004 - 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE.md at
// the top level directory of deal.II.
//
// ---------------------------------------------------------------------

#ifndef dealii_fe_bernardi_raugel_h
#define dealii_fe_bernardi_raugel_h

#include <deal.II/base/config.h>

#include <deal.II/base/geometry_info.h>
#include <deal.II/base/polynomial.h>
#include <deal.II/base/polynomials_bernardi_raugel.h>
#include <deal.II/base/tensor_product_polynomials.h>

#include <deal.II/fe/fe.h>
#include <deal.II/fe/fe_poly_tensor.h>

#include <vector>

DEAL_II_NAMESPACE_OPEN

/**
 * The Bernardi-Raugel element.
 *
 * <h3>Degrees of freedom</h3>
 * The BR1 element has <i>dim</i> degrees of freedom on each vertex and 1 on
 * each face. The shape functions are ordered by the $(Q_1)^d$ shape functions
 * supported on each vertex, increasing according to vertex ordering on the
 * element in GeometryInfo, then the bubble functions follow in the ordering
 * given in PolynomialsBernardiRaugel.
 *
 * This element only has 1 degree (degree $p=1$) because it yields an LBB stable
 * pair BR1-P0 for Stokes problems which is lower degree than the Taylor-Hood
 * element. The pair is sometimes referred to as an enriched P1-P0 element or a
 * reduced P2-P0 element.
 *
 * This element does not support hanging nodes or multigrid in the current
 * implementation.
 *
 * Some numerical experiments have shown that this element may converge with
 * first-order accuracy when using the BR1-Q0 pair for the mixed Laplace
 * equation in step-20.
 *
 */
template <int dim>
class FE_BernardiRaugel
  : public FE_PolyTensor<PolynomialsBernardiRaugel<dim>, dim>
{
public:
  /**
   * Constructor for the Bernardi-Raugel element of degree @p p. The only
   * supported degree is 1.
   *
   * @arg p: The degree of the element $p=1$ for $BR_1$.
   */
  FE_BernardiRaugel(const unsigned int p = 1);

  /**
   * Return a string that uniquely identifies a finite element. This class
   * returns <tt>FE_BR<dim>(degree)</tt>, with @p dim and @p degree replaced
   * by appropriate values.
   */
  virtual std::string
  get_name() const override;

  virtual std::unique_ptr<FiniteElement<dim, dim>>
  clone() const override;

  // documentation inherited from the base class
  virtual void
  convert_generalized_support_point_values_to_dof_values(
    const std::vector<Vector<double>> &support_point_values,
    std::vector<double> &              nodal_values) const override;

private:
  /**
   * Only for internal use. Its full name is @p get_dofs_per_object_vector
   * function and it creates the @p dofs_per_object vector that is needed
   * within the constructor to be passed to the constructor of @p
   * FiniteElementData.
   */
  static std::vector<unsigned int>
  get_dpo_vector();

  /**
   * Initialize the FiniteElement<dim>::generalized_support_points and
   * FiniteElement<dim>::generalized_face_support_points fields. Called from
   * the constructor. See the
   * @ref GlossGeneralizedSupport "glossary entry on generalized support points"
   * for more information.
   */
  void
  initialize_support_points();
};

DEAL_II_NAMESPACE_CLOSE

#endif
