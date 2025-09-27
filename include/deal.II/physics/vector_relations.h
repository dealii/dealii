// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2021 - 2022 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#ifndef dealii_vector_relations_h
#define dealii_vector_relations_h

#include <deal.II/base/config.h>

#include <deal.II/base/tensor.h>

#include <cmath>
#include <limits>

DEAL_II_NAMESPACE_OPEN


namespace Physics
{
  /**
   * Functions to compute relations between spatial vectors.
   */
  namespace VectorRelations
  {
    /**
     * Calculate the angle $\theta$ between two vectors @p a and @p b. The returned
     * angle will be in the range $[0, \pi]$.
     *
     * This function uses the geometric definition of the scalar product.
     * @f[
     *   \vec{a} \cdot \vec{b} = \|\vec{a}\| \|\vec{b}\| \cos(\theta)
     * @f]
     */
    template <int spacedim, typename Number>
    Number
    angle(const Tensor<1, spacedim, Number> &a,
          const Tensor<1, spacedim, Number> &b);

    /**
     * Calculate the angle $\theta$ between two vectors @p a and @p b, where both
     * vectors are located in a plane described by a normal vector @p axis.
     *
     * The angle computed by this function corresponds to the rotation angle
     * that would transform the vector @p a into the vector @p b around the vector
     * @p axis. Thus, contrary to the function above, we get a @em signed angle
     * which will be in the range $[-\pi, \pi]$.
     *
     * The vector @p axis needs to be a unit vector and be perpendicular to both
     * vectors @p a and @p b.
     *
     * This function uses the geometric definitions of both the scalar and cross
     * product.
     * @f{align*}{
     *   \vec{a} \cdot  \vec{b} &= \|\vec{a}\| \|\vec{b}\| \cos(\theta) \\
     *   \vec{a} \times \vec{b} &= \|\vec{a}\| \|\vec{b}\| \sin(\theta) \vec{n}
     * @f}
     * We can create the tangent of the angle using both products.
     * @f[
     *   \tan{\theta}
     *   = \frac{\sin(\theta)}{\cos(theta)}
     *   = \frac{(\vec{a} \times \vec{b}) \cdot \vec{n}}{\vec{a} \cdot \vec{b}}
     * @f]
     *
     * @note Only applicable for three-dimensional vectors `spacedim == 3`.
     */
    template <int spacedim, typename Number>
    Number
    signed_angle(const Tensor<1, spacedim, Number> &a,
                 const Tensor<1, spacedim, Number> &b,
                 const Tensor<1, spacedim, Number> &axis);
  } // namespace VectorRelations
} // namespace Physics



#ifndef DOXYGEN



template <int spacedim, typename Number>
inline Number
Physics::VectorRelations::angle(const Tensor<1, spacedim, Number> &a,
                                const Tensor<1, spacedim, Number> &b)
{
  const Number a_norm = a.norm();
  const Number b_norm = b.norm();
  Assert(a_norm > 1.e-12 * b_norm && a_norm > 1.e-12 * b_norm,
         ExcMessage("Both vectors need to be non-zero!"));

  Number argument = (a * b) / a_norm / b_norm;

  // std::acos returns nan if argument is out of domain [-1,+1].
  // if argument slightly overshoots these bounds, set it to the bound.
  // allow for 8*eps as a tolerance.
  if ((1. - std::abs(argument)) < 8. * std::numeric_limits<Number>::epsilon())
    argument = std::copysign(1., argument);

  return std::acos(argument);
}



template <int spacedim, typename Number>
inline Number
Physics::VectorRelations::signed_angle(const Tensor<1, spacedim, Number> &a,
                                       const Tensor<1, spacedim, Number> &b,
                                       const Tensor<1, spacedim, Number> &axis)
{
  Assert(spacedim == 3,
         ExcMessage("This function can only be used with spacedim==3!"));

  Assert(std::abs(axis.norm() - 1.) < 1.e-12,
         ExcMessage("The axial vector is not a unit vector."));
  Assert(std::abs(axis * a) < 1.e-12 * a.norm() &&
           std::abs(axis * b) < 1.e-12 * b.norm(),
         ExcMessage("The vectors are not perpendicular to the axial vector."));

  const Number dot = a * b;
  const Number det = axis * cross_product_3d(a, b);

  return std::atan2(det, dot);
}



#endif // DOXYGEN

DEAL_II_NAMESPACE_CLOSE

#endif
