// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2018 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#ifndef dealii_boost_adaptor_bounding_box_h
#define dealii_boost_adaptor_bounding_box_h

#include <deal.II/base/config.h>

#include <deal.II/base/bounding_box.h>

#include <deal.II/boost_adaptors/point.h>

DEAL_II_DISABLE_EXTRA_DIAGNOSTICS
#include <boost/geometry/core/coordinate_dimension.hpp>
#include <boost/geometry/core/coordinate_system.hpp>
#include <boost/geometry/core/coordinate_type.hpp>
#include <boost/geometry/core/cs.hpp>
#include <boost/geometry/core/tag.hpp>
#include <boost/geometry/strategies/strategies.hpp>
DEAL_II_ENABLE_EXTRA_DIAGNOSTICS


DEAL_II_NAMESPACE_OPEN
DEAL_II_NAMESPACE_CLOSE // Do not convert for module purposes

  namespace boost
{
  namespace geometry
  {
    namespace traits
    {
      /**
       * Tag adaptor for dealii::BoundingBox.
       */
      template <int dim, class Number>
      struct tag<dealii::BoundingBox<dim, Number>>
      {
        using type = box_tag;
      };

      /**
       * Point type adaptor for dealii::BoundingBox.
       */
      template <int dim, class Number>
      struct point_type<dealii::BoundingBox<dim, Number>>
      {
        using type = dealii::Point<dim, Number>;
      };

      /**
       * Access to the D-th coordinate of the lower left  corner of a
       * dealii::BoundingBox.
       */
      template <int dim, class Number, std::size_t D>
      struct indexed_access<dealii::BoundingBox<dim, Number>,
#if DEAL_II_BOOST_VERSION_GTE(1, 89, 0)
                            min_corner,
#else
                            // Until Boost 1.88, max_corner was a
                            // static variable in a header file, which
                            // we can't export in the module wrapper
                            // for Boost. Use the variable's numeric
                            // value instead.
                            /*min_corner*/ 0,
#endif
                            D>
      {
        /**
         * Getter function for the D-th coordinate of the lower left corner of
         * a dealii::BoundingBox.
         */
        static inline double
        get(const dealii::BoundingBox<dim, Number> &box)
        {
          return box.get_boundary_points().first[D];
        }

        /**
         * Setter function for the D-th coordinate of the lower left corner of
         * a dealii::BoundingBox.
         */
        static inline void
        set(dealii::BoundingBox<dim, Number> &box, Number value)
        {
          box.get_boundary_points().first[D] = value;
        }
      };

      /**
       * Access to the D-th coordinate of the upper right corner of a
       * dealii::BoundingBox.
       */
      template <int dim, class Number, std::size_t D>
      struct indexed_access<dealii::BoundingBox<dim, Number>,
#if DEAL_II_BOOST_VERSION_GTE(1, 89, 0)
                            max_corner,
#else
                            // Until Boost 1.88, max_corner was a
                            // static variable in a header file, which
                            // we can't export in the module wrapper
                            // for Boost. Use the variable's numeric
                            // value instead.
                            /*max_corner*/ 1,
#endif
                            D>
      {
        /**
         * Getter function for the D-th coordinate of the upper right corner of
         * a dealii::BoundingBox.
         */
        static inline double
        get(const dealii::BoundingBox<dim, Number> &box)
        {
          return box.get_boundary_points().second[D];
        }

        /**
         * Setter function for the D-th coordinate of the upper right corner of
         * a dealii::BoundingBox.
         */
        static inline void
        set(dealii::BoundingBox<dim, Number> &box, Number value)
        {
          box.get_boundary_points().second[D] = value;
        }
      };
    } // namespace traits
  }   // namespace geometry
} // namespace boost

DEAL_II_NAMESPACE_OPEN // Do not convert for module purposes
  DEAL_II_NAMESPACE_CLOSE

#endif
