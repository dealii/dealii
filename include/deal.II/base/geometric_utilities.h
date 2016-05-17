// ---------------------------------------------------------------------
//
// Copyright (C) 2016 by the deal.II authors
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

#ifndef dealii__geometric_utilities_h
#define dealii__geometric_utilities_h

#include <deal.II/base/config.h>
#include <deal.II/base/point.h>
#include <deal.II/base/std_cxx11/array.h>


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
     * Returns spherical coordinates of a Cartesian point @p point.
     * The returned array is filled with radius, azimuth angle $\in [0,2 \pi)$
     * and polar/inclination angle $ \in [0,\pi]$ (ommited in 2D).
     */
    template <int dim>
    std_cxx11::array<double,dim>
    to_spherical(const Point<dim> &point);

    /**
     * Return the Cartesian coordinates of a spherical point defined by @p scoord
     * which is filled with radius $\in [0,\infty)$, azimuth angle
     * $\in [0,2 \pi)$ and polar/inclination angle $\in [0,\pi]$
     * (ommited in 2D).
     */
    template <std::size_t dim>
    Point<dim>
    from_spherical(const std_cxx11::array<double,dim> &scoord);

  }
}

DEAL_II_NAMESPACE_CLOSE

#endif
