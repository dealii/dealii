// ---------------------------------------------------------------------
//
// Copyright (C) 2018 - 2020 by the deal.II authors
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

#include <deal.II/boost_adaptors/point.h>


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
          std::pair<dealii::Point<dim, Number>, dealii::Point<dim, Number>>
            corner_points = box.get_boundary_points();

          // The caller of this function says that they want the first
          // point updated, but sometimes this creates an invalid
          // bounding box. Check for this, and if necessary make sure
          // that the points remains sorted correctly
          corner_points.first[D] = value;
          if (corner_points.first[D] > corner_points.second[D])
            std::swap(corner_points.first[D], corner_points.second[D]);

          box = dealii::BoundingBox<dim, Number>(corner_points);
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
          std::pair<dealii::Point<dim, Number>, dealii::Point<dim, Number>>
            corner_points = box.get_boundary_points();

          // The caller of this function says that they want the second
          // point updated, but sometimes this creates an invalid
          // bounding box. Check for this, and if necessary make sure
          // that the points remains sorted correctly
          corner_points.second[D] = value;
          if (corner_points.first[D] > corner_points.second[D])
            std::swap(corner_points.first[D], corner_points.second[D]);

          box = dealii::BoundingBox<dim, Number>(corner_points);
        }
      };
    } // namespace traits
  }   // namespace geometry
} // namespace boost

#endif
