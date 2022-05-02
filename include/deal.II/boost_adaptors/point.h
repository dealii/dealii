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

#ifndef dealii_boost_adaptor_point_h
#define dealii_boost_adaptor_point_h

#include <deal.II/base/config.h>

#include <deal.II/base/point.h>

namespace boost
{
  namespace geometry
  {
    namespace traits
    {
      /**
       * Tag adaptor for dealii::Point.
       */
      template <int dim, class Number>
      struct tag<dealii::Point<dim, Number>>
      {
        using type = point_tag;
      };

      /**
       * Coordinate type adaptor for dealii::Point.
       */
      template <int dim, class Number>
      struct coordinate_type<dealii::Point<dim, Number>>
      {
        using type = Number;
      };

      /**
       * Coordinate system adaptor for dealii::Point. By default, we assume
       * that a dealii Point is cartesian point.
       */
      template <int dim, class Number>
      struct coordinate_system<dealii::Point<dim, Number>>
      {
        using type = cs::cartesian;
      };

      /**
       * Dimension adaptor.
       */
      template <int dim, class Number>
      struct dimension<dealii::Point<dim, Number>> : boost::mpl::int_<dim>
      {};

      /**
       * Getter function for D-th coordinate of a dealii Point.
       */
      template <std::size_t D, int dim, class Number>
      struct access<dealii::Point<dim, Number>, D>
      {
        static inline double
        get(dealii::Point<dim, Number> const &p)
        {
          return p[D];
        }

        /**
         * Setter function for D-th coordinate of a dealii Point.
         */
        static inline void
        set(dealii::Point<dim, Number> &p, Number value)
        {
          p[D] = value;
        }
      };
    } // namespace traits
  }   // namespace geometry
} // namespace boost



#endif
