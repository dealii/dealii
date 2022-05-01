// ---------------------------------------------------------------------
//
// Copyright (C) 2022 by the deal.II authors
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

#ifndef dealii_cgal_utilities_h
#define dealii_cgal_utilities_h

#include <deal.II/base/config.h>

#ifdef DEAL_II_WITH_CGAL
#  include <CGAL/Cartesian.h>
#  include <CGAL/Simple_cartesian.h>

DEAL_II_NAMESPACE_OPEN
/**
 * Interface to the Computational Geometry Algorithm Library (CGAL).
 *
 * CGAL is a software project that provides easy access to efficient and
 * reliable geometric algorithms. The library offers data structures and
 * algorithms like triangulations, Voronoi diagrams, Boolean operations on
 * polygons and polyhedra, point set processing, arrangements of curves, surface
 * and volume mesh generation, geometry processing,  alpha shapes, convex hull
 * algorithms, shape reconstruction, AABB and KD trees...
 *
 * You can learn more about the CGAL library at https://www.cgal.org/
 */
namespace CGALWrappers
{
  /**
   * Convert from deal.II Point to any compatible CGAL point.
   *
   * @tparam CGALPointType Any of the CGAL point types
   * @tparam dim Dimension of the point
   * @param [in] p An input deal.II Point<dim>
   * @return CGALPointType A CGAL point
   */
  template <typename CGALPointType, int dim>
  inline CGALPointType
  dealii_point_to_cgal_point(const dealii::Point<dim> &p);

  /**
   * Convert from various CGAL point types to deal.II Point.
   *
   * @tparam dim Dimension of the point
   * @tparam CGALPointType Any of the CGAL point types
   * @param p An input CGAL point type
   * @return dealii::Point<dim> The corresponding deal.II point.
   */
  template <int dim, typename CGALPointType>
  inline dealii::Point<dim>
  cgal_point_to_dealii_point(const CGALPointType &p);
} // namespace CGALWrappers

#  ifndef DOXYGEN
// Template implementations
namespace CGALWrappers
{
  template <typename CGALPointType, int dim>
  inline CGALPointType
  dealii_point_to_cgal_point(const dealii::Point<dim> &p)
  {
    constexpr int cdim = CGALPointType::Ambient_dimension::value;
    static_assert(dim <= cdim, "Only dim <= cdim supported");
    if constexpr (cdim == 1)
      return CGALPointType(p[0]);
    else if constexpr (cdim == 2)
      return CGALPointType(p[0], dim > 1 ? p[1] : 0);
    else if constexpr (cdim == 3)
      return CGALPointType(p[0], dim > 1 ? p[1] : 0, dim > 2 ? p[2] : 0);
    else
      Assert(false, dealii::ExcNotImplemented());
    return CGALPointType();
  }



  template <int dim, typename CGALPointType>
  inline dealii::Point<dim>
  cgal_point_to_dealii_point(const CGALPointType &p)
  {
    constexpr int cdim = CGALPointType::Ambient_dimension::value;
    if constexpr (dim == 1)
      return dealii::Point<dim>(CGAL::to_double(p.x()));
    else if constexpr (dim == 2)
      return dealii::Point<dim>(CGAL::to_double(p.x()),
                                cdim > 1 ? CGAL::to_double(p.y()) : 0);
    else if constexpr (dim == 3)
      return dealii::Point<dim>(CGAL::to_double(p.x()),
                                cdim > 1 ? CGAL::to_double(p.y()) : 0,
                                cdim > 2 ? CGAL::to_double(p.z()) : 0);
    else
      Assert(false, dealii::ExcNotImplemented());
  }
} // namespace CGALWrappers
#  endif

DEAL_II_NAMESPACE_CLOSE

#endif
#endif
