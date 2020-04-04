// ---------------------------------------------------------------------
//
// Copyright (C) 2018 - 2019 by the deal.II authors
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

#ifndef dealii_boost_adaptor_bounding_box_h
#define dealii_boost_adaptor_bounding_box_h

#include <deal.II/base/config.h>

#include <deal.II/base/bounding_box.h>

DEAL_II_DISABLE_EXTRA_DIAGNOSTICS
#include <deal.II/boost_adaptors/point.h>
DEAL_II_ENABLE_EXTRA_DIAGNOSTICS


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
      struct indexed_access<dealii::BoundingBox<dim, Number>, min_corner, D>
      {
        /**
         * Getter function for the D-th coordinate of the lower left corner of
         * a dealii::BoundingBox.
         */
        static inline double
        get(dealii::BoundingBox<dim, Number> const &box)
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
      struct indexed_access<dealii::BoundingBox<dim, Number>, max_corner, D>
      {
        /**
         * Getter function for the D-th coordinate of the upper right corner of
         * a dealii::BoundingBox.
         */
        static inline double
        get(dealii::BoundingBox<dim, Number> const &box)
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

#endif
