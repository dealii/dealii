// ---------------------------------------------------------------------
//
// Copyright (C) 1998 - 2018 by the deal.II authors
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

#ifndef dealii_non_matching_quadrature_generator
#define dealii_non_matching_quadrature_generator

#include <deal.II/base/config.h>

#include <deal.II/base/point.h>

#include <type_traits>

DEAL_II_NAMESPACE_OPEN

namespace NonMatching
{
  namespace internal
  {
    namespace QuadratureGeneratorImplementation
    {
      /**
       * Returns the center $(0.5, ..., 0.5)$ of the unit hypercube:
       * $[0,1]^{dim}$.
       */
      template <int dim>
      Point<dim>
      unit_hypercube_center();


      /**
       * Describes a dim-dimensional hypercube aligned with the coordinate
       * axes. That is, a region
       *
       * @f[
       * [x_0-L/2, x_0+L/2] \times ... \times [x_{dim-1}-L/2, x_{dim-1}+L/2]
       * @f]
       *
       * where L is the side length and $(x_0, ..., x_{dim-1})$ are the
       * coordinates of the center of the hypercube.
       *
       * The purpose of this class is to help us keep track of where we are on
       * the reference cell when generating immersed quadrature rules. There
       * are two main things we want to be able to do.
       *
       * 1. Split a hypercube into its $2^{dim}$ children.
       * 2. Take a cross section orthogonal to a given direction and get this
       * as a Hypercube<dim-1>.
       *
       * These are needed because the algorithm is allowed to recurse both
       * over children and over dimensions.
       *
       * In 3D, the 2 coordinates of the cross section of Hypercube<3> can be
       * ordered in 2 ways. That is, if we take the cross section orthogonal
       * to the y direction we could either order a 3D-coordinate into a
       * 2D-coordinate as $(x,z)$ or as $(z,x)$. This class uses the second
       * convention, corresponding to the coordinates being ordered cyclicly
       * $x \rightarrow y \rightarrow z \rightarrow x \rightarrow ... $
       * To be precise, if we take a cross section:
       *
       * | Orthogonal to | Cross section coordinates ordered as |
       * |:-------------:|:------------------------------------:|
       * |      x        |               (y, z)                 |
       * |      y        |               (z, x)                 |
       * |      z        |               (x, y)                 |
       *
       * This is according to the convention set by the function
       * <code>coordinate_to_one_dim_higher</code>.
       */
      template <int dim>
      class Hypercube
      {
      public:
        /**
         * Constructor, takes the center of the hypercube and its side length.
         * Passing no inarguments creates a hypercube between 0 and 1 in all
         * directions: $[0, 1]^{dim}$.
         */
        Hypercube(
          const Point<dim> &center_coordinate = unit_hypercube_center<dim>(),
          const double      length_of_side    = 1);

        /**
         * Return the lower bound of the hypercube in @p direction.
         */
        double
        lower_bound(const unsigned int direction) const;

        /**
         * Return the upper bound of the hypercube in @p direction.
         */
        double
        upper_bound(const unsigned int direction) const;

        /**
         * Return the bounds of the hypercube in @p direction, as a one-dimensional
         * hypercube.
         */
        Hypercube<1>
        bounds(const unsigned int direction) const;

        /**
         * Returns the cross section of the hypercube orthogonal to @p direction.
         * This is a hypercube in one dimension lower.
         */
        Hypercube<dim - 1>
        cross_section(const unsigned int direction) const;

        /**
         * Returns the point in the center of the hypercube.
         */
        const Point<dim> &
        center() const;

        /**
         * Returns the length of the side of the hypercube.
         */
        double
        side_length() const;

        /**
         * Returns the volume of the hypercube.
         */
        double
        volume() const;

        /**
         * Returns the indexth vertex of the hypercube.
         */
        Point<dim>
        vertex(const unsigned int index) const;

        /**
         * Returns the indexth child of the hypercube. Child is meant in the
         * same way as for a cell.
         */
        Hypercube<dim>
        child(const unsigned int index) const;

      private:
        /**
         * Coordinate in the center of the hypercube.
         */
        const Point<dim> center_coordinate;

        /**
         * Side length of the hypercube.
         */
        const double length_of_side;
      };


      /**
       * This function defines a convention for how coordinates in dim
       * dimensions should translate to the coordinates in dim + 1 dimensions,
       * when one of the coordinates in dim + 1 dimensions is locked to a given
       * value.
       *
       * The convention is the following: Starting from the locked coordinate we
       * store the lower dimensional coordinates consecutively and wrapping
       * around when going over the dimension. This relationship can be
       * described by the following tables:
       *
       *                        2D
       * ------------------------------------------------
       * | locked in 2D | 1D coordinate | 2D coordinate |
       * |--------------------------------------------- |
       * |     x0       |      (a)      |   (x0,  a)    |
       * |     x1       |      (a)      |   (a , x1)    |
       * ------------------------------------------------
       *
       *                       3D
       * --------------------------------------------------
       * | locked in 3D | 2D coordinates | 3D coordinates |
       * |----------------------------------------------- |
       * |     x0       |    (a, b)      | (x0,  a,  b)   |
       * |     x1       |    (a, b)      | ( b, x1,  a)   |
       * |     x2       |    (a, b)      | ( a,  b, x2)   |
       * --------------------------------------------------
       *
       * Given a locked coordinate, this function maps a coordinate index in dim
       * dimension to a coordinate index in dim + 1 dimensions.
       *
       * @param locked_coordinate should be in the range [0, dim+1).
       * @param coordiante_in_dim should be in the range [0, dim).
       * @return A coordinate index in the range [0, dim+1)
       */
      template <int dim>
      inline unsigned int
      coordinate_to_one_dim_higher(const unsigned int locked_coordinate,
                                   const unsigned int coordiante_in_dim)
      {
        AssertIndexRange(locked_coordinate, dim + 1);
        AssertIndexRange(coordiante_in_dim, dim);
        return (locked_coordinate + coordiante_in_dim + 1) % (dim + 1);
      }



      /* -------------- declaration of explicit specializations ------------- */

      template <>
      Hypercube<-1>
      Hypercube<0>::cross_section(const unsigned int) const;

      template <>
      Point<0>
      Hypercube<0>::vertex(const unsigned int) const;

      template <>
      Hypercube<0>
      Hypercube<0>::child(const unsigned int) const;

    } // namespace QuadratureGeneratorImplementation
  }   // namespace internal

} // namespace NonMatching
DEAL_II_NAMESPACE_CLOSE

#endif /* dealii_non_matching_quadrature_generator */
