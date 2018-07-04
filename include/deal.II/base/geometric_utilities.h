// ---------------------------------------------------------------------
//
// Copyright (C) 2016 - 2018 by the deal.II authors
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

#ifndef dealii_geometric_utilities_h
#define dealii_geometric_utilities_h

#include <deal.II/base/config.h>

#include <deal.II/base/point.h>

#include <array>


DEAL_II_NAMESPACE_OPEN


/**
 * A namespace for geometric utility functions that are not particularly specific
 * to finite element computing or numerical programs, but nevertheless are needed
 * in various contexts when writing applications.
 *
 * @ingroup utilities
 * @author Denis Davydov, 2016
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
     * and polar/inclination angle $ \in [0,\pi]$ (omitted in 2D).
     *
     * In 3D the transformation is given by
     * @f{align*}{
     *  r &= \sqrt{x^2+y^2+z^2} \\
     *  \theta &= {\rm atan}(y/x) \\
     *  \phi &= {\rm acos} (z/r)
     * @f}
     */
    template <int dim>
    std::array<double, dim>
    to_spherical(const Point<dim> &point);

    /**
     * Return the Cartesian coordinates of a spherical point defined by @p scoord
     * which is filled with radius $r \in [0,\infty)$, azimuth angle
     * $\theta \in [0,2 \pi)$ and polar/inclination angle $\phi \in [0,\pi]$
     * (omitted in 2D).
     *
     * In 3D the transformation is given by
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
