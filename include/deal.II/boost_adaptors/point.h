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

#ifndef dealii_boost_adaptor_point_h
#define dealii_boost_adaptor_point_h

#include <deal.II/base/config.h>

#include <deal.II/base/point.h>

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
        get(const dealii::Point<dim, Number> &p)
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

DEAL_II_NAMESPACE_OPEN // Do not convert for module purposes
  DEAL_II_NAMESPACE_CLOSE


#endif
