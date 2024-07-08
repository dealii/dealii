// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2018 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

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
 * This class implements the non-standard Bernardi-Raugel (BR) element
 * that can be used as one part of a stable velocity/pressure pair for
 * the Stokes equation. The BR element can be seen as either an
 * enriched version of the $Q_1^d$ element with added bubble functions
 * on each edge (in 2d) or face (in 3d), or as a reduced version of
 * the $Q_2^d$ element. It addresses the fact that the $Q_1^d\times
 * Q_0$ combination is not inf-sup stable (requiring a larger velocity
 * space), and that the $Q_2^d\times Q_1$ combination is stable but
 * sub-optimal since the velocity space is too large relative to the
 * pressure space to provide additional accuracy commensurate with the
 * cost of the large number of velocity unknowns.
 *
 * The element was introduced in the paper @cite BR85.
 *
 *
 * <h3>Degrees of freedom</h3>
 *
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
 */
template <int dim>
class FE_BernardiRaugel : public FE_PolyTensor<dim>
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
    std::vector<double>               &nodal_values) const override;

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

  /**
   * Initialize the permutation pattern and the pattern of sign change.
   *
   * @note This function is not fully filled with the correct implementation
   * yet. It needs to be consistently implemented in a future release to work
   * on meshes that contain cells with flipped faces.
   */
  void
  initialize_quad_dof_index_permutation_and_sign_change();
};

DEAL_II_NAMESPACE_CLOSE

#endif
