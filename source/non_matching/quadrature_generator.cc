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

#include <deal.II/base/exceptions.h>
#include <deal.II/base/geometry_info.h>

#include <deal.II/non_matching/quadrature_generator.h>


DEAL_II_NAMESPACE_OPEN

namespace NonMatching
{
  namespace internal
  {
    namespace QuadratureGeneratorImplementation
    {
      template <int dim>
      Point<dim>
      unit_hypercube_center()
      {
        Point<dim> center;
        for (unsigned int i = 0; i < dim; ++i)
          center[i] = .5;

        return center;
      }



      template <int dim>
      Hypercube<dim>::Hypercube(const Point<dim> &center_coordinate,
                                const double      length_of_side)
        : center_coordinate(center_coordinate)
        , length_of_side(length_of_side)
      {
        Assert(0 < dim, ExcImpossibleInDim(dim));
      }



      template <int dim>
      double
      Hypercube<dim>::lower_bound(const unsigned int direction) const
      {
        AssertIndexRange(direction, dim);

        return center_coordinate[direction] - length_of_side / 2;
      }



      template <int dim>
      double
      Hypercube<dim>::upper_bound(const unsigned int direction) const
      {
        AssertIndexRange(direction, dim);

        return center_coordinate[direction] + length_of_side / 2;
      }



      template <int dim>
      Hypercube<1>
      Hypercube<dim>::bounds(const unsigned int direction) const
      {
        AssertIndexRange(direction, dim);

        const Point<1> center(center_coordinate[direction]);
        return Hypercube<1>(center, side_length());
      }



      template <int dim>
      const Point<dim> &
      Hypercube<dim>::center() const
      {
        return center_coordinate;
      }



      template <int dim>
      double
      Hypercube<dim>::side_length() const
      {
        return length_of_side;
      }



      template <int dim>
      double
      Hypercube<dim>::volume() const
      {
        return std::pow(length_of_side, dim);
      }



      template <int dim>
      Hypercube<dim - 1>
      Hypercube<dim>::cross_section(const unsigned int direction) const
      {
        AssertIndexRange(direction, dim);

        Point<dim - 1> center_of_cross_section;
        for (unsigned int d = 0; d < dim - 1; ++d)
          {
            const int index_to_write_from =
              coordinate_to_one_dim_higher<dim - 1>(direction, d);
            center_of_cross_section[d] = center_coordinate[index_to_write_from];
          }

        return Hypercube<dim - 1>(center_of_cross_section, length_of_side);
      }



      template <int dim>
      Point<dim>
      Hypercube<dim>::vertex(const unsigned int index) const
      {
        AssertIndexRange(index, GeometryInfo<dim>::vertices_per_cell);

        Point<dim> bottom_corner;
        for (unsigned int i = 0; i < dim; ++i)
          bottom_corner[i] = center_coordinate[i] - length_of_side / 2;

        const Point<dim> point =
          bottom_corner +
          length_of_side * GeometryInfo<dim>::unit_cell_vertex(index);

        return point;
      }



      template <int dim>
      Hypercube<dim>
      Hypercube<dim>::child(const unsigned int index) const
      {
        AssertIndexRange(index, GeometryInfo<dim>::max_children_per_cell);

        const double child_side_length = length_of_side / 2;

        const Point<dim> unit_cell_vertex =
          GeometryInfo<dim>::unit_cell_vertex(index);

        Point<dim> child_center;
        for (unsigned int i = 0; i < dim; ++i)
          child_center[i] = center_coordinate[i] +
                            child_side_length * (unit_cell_vertex[i] - .5);

        return Hypercube<dim>(child_center, child_side_length);
      }



      template Point<0>
      unit_hypercube_center<0>();
      template Point<1>
      unit_hypercube_center<1>();
      template Point<2>
      unit_hypercube_center<2>();
      template Point<3>
      unit_hypercube_center<3>();

      template class Hypercube<0>;
      template class Hypercube<1>;
      template class Hypercube<2>;
      template class Hypercube<3>;

    } // namespace QuadratureGeneratorImplementation
  }   // namespace internal



} // namespace NonMatching
DEAL_II_NAMESPACE_CLOSE
