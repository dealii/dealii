// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2016 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#ifndef dealii_geometric_utilities_h
#define dealii_geometric_utilities_h

#include <deal.II/base/config.h>

#include <deal.II/base/point.h>

#include <array>


DEAL_II_NAMESPACE_OPEN


/**
 * A namespace for geometric utility functions that are not particularly
 * specific to finite element computing or numerical programs, but nevertheless
 * are needed in various contexts when writing applications.
 *
 * @ingroup utilities
 */
namespace GeometricUtilities
{
  /**
   * A namespace for coordinate transformations.
   */
  namespace Coordinates
  {
    /**
     * Return spherical coordinates of a Cartesian point @p point.
     * The returned array is filled with radius, azimuth angle $\in [0,2 \pi)$
     * and polar/inclination angle $ \in [0,\pi]$ (omitted in 2d).
     *
     * In 3d the transformation is given by
     * @f{align*}{
     *  r &= \sqrt{x^2+y^2+z^2} \\
     *  \theta &= {\rm atan}(y/x) \\
     *  \phi &= {\rm acos} (z/r)
     * @f}
     *
     * The use of this function is demonstrated in step-75.
     */
    template <int dim>
    std::array<double, dim>
    to_spherical(const Point<dim> &point);

    /**
     * Return the Cartesian coordinates of a spherical point defined by @p scoord
     * which is filled with radius $r \in [0,\infty)$, azimuth angle
     * $\theta \in [0,2 \pi)$ and polar/inclination angle $\phi \in [0,\pi]$
     * (omitted in 2d).
     *
     * In 3d the transformation is given by
     * @f{align*}{
     *  x &= r\, \cos(\theta) \, \sin(\phi) \\
     *  y &= r\, \sin(\theta) \, \sin(\phi) \\
     *  z &= r\, \cos(\phi)
     * @f}
     */
    template <std::size_t dim>
    Point<dim>
    from_spherical(const std::array<double, dim> &scoord);

  } // namespace Coordinates
} // namespace GeometricUtilities

DEAL_II_NAMESPACE_CLOSE

#endif
