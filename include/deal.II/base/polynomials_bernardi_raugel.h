// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2018 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


#ifndef dealii_polynomials_bernardi_raugel_h
#define dealii_polynomials_bernardi_raugel_h

#include <deal.II/base/config.h>

#include <deal.II/base/geometry_info.h>
#include <deal.II/base/point.h>
#include <deal.II/base/polynomial.h>
#include <deal.II/base/polynomial_space.h>
#include <deal.II/base/tensor.h>
#include <deal.II/base/tensor_polynomials_base.h>
#include <deal.II/base/tensor_product_polynomials.h>

#include <vector>

DEAL_II_NAMESPACE_OPEN


/**
 * This class implements the Bernardi-Raugel polynomials similarly to the
 * description in the <i>Mathematics of Computation</i> paper from 1985 by
 * Christine Bernardi and Geneviève Raugel.
 *
 * The Bernardi-Raugel polynomials are originally defined as an enrichment
 * of the $(P_1)^d$ elements on simplicial meshes for Stokes problems by the
 * addition of bubble functions, yielding a locking-free finite element which
 * is a subset of $(P_2)^d$ elements. This implementation is an enrichment of
 * $(Q_1)^d$ elements which is a subset of $(Q_2)^d$ elements for
 * quadrilateral and hexahedral meshes.
 *
 * The $BR_1$ bubble functions are defined to have magnitude 1 at the center
 * of face $e_i$ and direction $\mathbf{n}_i$ normal to face $e_i$, and
 * magnitude 0 on all other vertices and faces. Ordering is consistent with
 * the face numbering in GeometryInfo. The vector $\mathbf{n}_i$ points in
 * the positive axis direction and not necessarily normal to the element for
 * consistent orientation across edges.
 *<dl>
 *   <dt> 2d bubble functions (in order)
 *   <dd> $x=0$ edge: $\mathbf{p}_1 = \mathbf{n}_1 (1-x)(y)(1-y)$
 *
 *        $x=1$ edge: $\mathbf{p}_2 = \mathbf{n}_2 (x)(y)(1-y)$
 *
 *        $y=0$ edge: $\mathbf{p}_3 = \mathbf{n}_3 (x)(1-x)(1-y)$
 *
 *        $y=1$ edge: $\mathbf{p}_4 = \mathbf{n}_4 (x)(1-x)(y)$
 *
 *   <dt> 3d bubble functions (in order)
 *   <dd> $x=0$ edge: $\mathbf{p}_1 = \mathbf{n}_1 (1-x)(y)(1-y)(z)(1-z)$
 *
 *        $x=1$ edge: $\mathbf{p}_2 = \mathbf{n}_2 (x)(y)(1-y)(z)(1-z)$
 *
 *        $y=0$ edge: $\mathbf{p}_3 = \mathbf{n}_3 (x)(1-x)(1-y)(z)(1-z)$
 *
 *        $y=1$ edge: $\mathbf{p}_4 = \mathbf{n}_4 (x)(1-x)(y)(z)(1-z)$
 *
 *        $z=0$ edge: $\mathbf{p}_5 = \mathbf{n}_5 (x)(1-x)(y)(1-y)(1-z)$
 *
 *        $z=1$ edge: $\mathbf{p}_6 = \mathbf{n}_6 (x)(1-x)(y)(1-y)(z)$
 *
 *</dl>
 *
 * Then the $BR_1(E)$ polynomials are defined on quadrilaterals and hexahedra
 * by $BR_1(E) = Q_1(E) \oplus \mbox{span}\{\mathbf{p}_i, i=1,...,2d\}$.
 *
 *
 * @ingroup Polynomials
 */
template <int dim>
class PolynomialsBernardiRaugel : public TensorPolynomialsBase<dim>
{
public:
  /**
   * Constructor. Creates all basis functions for Bernardi-Raugel polynomials
   * of given degree.
   *
   * @arg k The degree of the Bernardi-Raugel-space, which is currently
   * limited to the case <tt>k=1</tt>.
   */
  PolynomialsBernardiRaugel(const unsigned int k);

  /**
   * Return the name of the space, which is <tt>BernardiRaugel</tt>.
   */
  std::string
  name() const override;

  /**
   * Compute the value and derivatives of each Bernardi-Raugel
   * polynomial at @p unit_point.
   *
   * The size of the vectors must either be zero or equal <tt>n()</tt>.  In
   * the first case, the function will not compute these values.
   */
  void
  evaluate(const Point<dim>            &unit_point,
           std::vector<Tensor<1, dim>> &values,
           std::vector<Tensor<2, dim>> &grads,
           std::vector<Tensor<3, dim>> &grad_grads,
           std::vector<Tensor<4, dim>> &third_derivatives,
           std::vector<Tensor<5, dim>> &fourth_derivatives) const override;

  /**
   * Return the number of polynomials in the space <tt>BR(degree)</tt> without
   * requiring to build an object of PolynomialsBernardiRaugel. This is
   * required by the FiniteElement classes.
   */
  static unsigned int
  n_polynomials(const unsigned int k);

  /**
   * @copydoc TensorPolynomialsBase::clone()
   */
  virtual std::unique_ptr<TensorPolynomialsBase<dim>>
  clone() const override;

private:
  /**
   * An object representing the polynomial space of Q
   * functions which forms the <tt>BR</tt> polynomials through
   * outer products of these with the corresponding unit ijk
   * vectors.
   */
  const AnisotropicPolynomials<dim> polynomial_space_Q;

  /**
   * An object representing the polynomial space of bubble
   * functions which forms the <tt>BR</tt> polynomials through
   * outer products of these with the corresponding normals.
   */
  const AnisotropicPolynomials<dim> polynomial_space_bubble;

  /**
   * A static member function that creates the polynomial space we use to
   * initialize the #polynomial_space_Q member variable.
   */
  static std::vector<std::vector<Polynomials::Polynomial<double>>>
  create_polynomials_Q();

  /**
   * A static member function that creates the polynomial space we use to
   * initialize the #polynomial_space_bubble member variable.
   */
  static std::vector<std::vector<Polynomials::Polynomial<double>>>
  create_polynomials_bubble();
};


template <int dim>
inline std::string
PolynomialsBernardiRaugel<dim>::name() const
{
  return "BernardiRaugel";
}


DEAL_II_NAMESPACE_CLOSE

#endif
