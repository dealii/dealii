// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2021 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#ifndef dealii_fe_fe_p_bubbles_h
#define dealii_fe_fe_p_bubbles_h

#include <deal.II/base/config.h>

#include <deal.II/base/polynomials_barycentric.h>

#include <deal.II/fe/fe_simplex_p.h>

DEAL_II_NAMESPACE_OPEN

/**
 * @brief Enriched version of FE_SimplexP that can be used with nodal
 * quadrature.
 *
 * Many explicit time integration schemes require solving a @ref GlossMassMatrix "mass matrix" at
 * each time step. There are various ways around this requirement - for
 * example, step-48 replaces the mass matrix with a diagonal approximation,
 * which makes the solution step trivial. In step-48, and also commonly for
 * tensor-product elements, this is done by computing the mass matrix with a
 * lower-order quadrature point based on the nodes of the finite element
 * (i.e., the nodal quadrature rule one obtains by using the shape functions
 * as an interpolatory basis).
 *
 * A major drawback of standard simplex-based finite elements is that they
 * cannot be used with nodal quadrature since some of the quadrature weights
 * end up being either zero or negative, resulting in either an unsolvable or
 * unstable approximation to the mass matrix. For example: the shape functions
 * of FE_SimplexP<2>(2) with support points at vertices have mean values of
 * zero so that element cannot be used with mass lumping.
 *
 * This element avoids this issue by replacing the shape functions of
 * FE_SimplexP with an augmented space amenable to the construction of nodal
 * quadrature rules. For example, on the triangle a single basis function is
 * added corresponding to interpolation at the centroid (and all other basis
 * functions are updated to preserve the partition of unity property). This
 * results in shape functions with positive means (i.e., a valid nodal
 * quadrature formula). Similarly, in 3d, the polynomial space of
 * FE_SimplexP<3>(2) is enriched with five additional degrees of freedom (where
 * four have support points at face centroids and one has a support point at the
 * centroid) to enable construction of valid nodal quadrature rule.
 *
 * Since this FE space includes bubbles (i.e., extra functions which are
 * nonzero only on element interiors), the polynomial degrees of the component
 * basis functions are higher than the actual approximation degree of the
 * element. For example, with a constructor argument <code>degree = 2</code>
 * in 3d, the polynomials are in fact cubic (degree 3) but the order of the
 * approximation is the same as if we were using quadratic (degree 2) finite
 * elements.
 *
 * The 2d quadratic element was first described in @cite fried1975finite. The
 * 3d quadratic element implemented here was first described in
 * @cite Geevers_2018. Higher degree elements amendable to lumping exist but
 * are not yet implemented in this class.
 */
template <int dim, int spacedim = dim>
class FE_SimplexP_Bubbles : public FE_SimplexPoly<dim, spacedim>
{
public:
  /**
   * Constructor, taking the approximation degree as an argument. The
   * polynomial space is typically one degree higher than the approximation
   * space for this element: see the general documentation of this class for
   * more information.
   *
   * @note For <code>degree == 1</code> this element is equivalent to
   * FE_SimplexP(1).
   */
  FE_SimplexP_Bubbles(const unsigned int degree);

  /**
   * @copydoc dealii::FiniteElement::clone()
   */
  virtual std::unique_ptr<FiniteElement<dim, spacedim>>
  clone() const override;

  /**
   * Return a string that uniquely identifies a finite element. This class
   * returns <tt>FE_SimplexP_Bubbles<dim,spacedim>(degree)</tt>, with
   * @p dim, @p spacedim, and @p degree replaced by appropriate values. As
   * usual, @p spacedim is omitted in the codimension zero case.
   */
  virtual std::string
  get_name() const override;

protected:
  /**
   * Degree of the approximation (i.e., the constructor argument).
   */
  unsigned int approximation_degree;
};

DEAL_II_NAMESPACE_CLOSE

#endif
