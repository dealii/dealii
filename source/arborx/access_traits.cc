// ---------------------------------------------------------------------
//
// Copyright (C) 2021 by the deal.II authors
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

#include <deal.II/arborx/access_traits.h>

#ifdef DEAL_II_WITH_ARBORX

DEAL_II_NAMESPACE_OPEN

namespace ArborXWrappers
{
  // ------------------- PointIntersect ------------------- //
  template <int dim, typename Number>
  PointIntersect::PointIntersect(
    const std::vector<dealii::Point<dim, Number>> &dim_points)
  {
    static_assert(dim != 1, "dim equal to one is not supported.");

    const unsigned int size = dim_points.size();
    points.reserve(size);
    for (unsigned int i = 0; i < size; ++i)
      {
        points.emplace_back(static_cast<float>(dim_points[i][0]),
                            static_cast<float>(dim_points[i][1]),
                            dim == 2 ? 0.f :
                                       static_cast<float>(dim_points[i][2]));
      }
  }



  std::size_t
  PointIntersect::size() const
  {
    return points.size();
  }



  const dealii::Point<3, float> &
  PointIntersect::get(unsigned int i) const
  {
    return points[i];
  }



  // ------------------- BoundingBoxIntersect ------------------- //
  template <int dim, typename Number>
  BoundingBoxIntersect::BoundingBoxIntersect(
    const std::vector<dealii::BoundingBox<dim, Number>> &bb)
  {
    const unsigned int size = bb.size();
    bounding_boxes.reserve(size);
    dealii::Point<3, float> min_corner_arborx(0., 0., 0.);
    dealii::Point<3, float> max_corner_arborx(0., 0., 0.);
    for (unsigned int i = 0; i < size; ++i)
      {
        auto boundary_points                  = bb[i].get_boundary_points();
        dealii::Point<dim, Number> min_corner = boundary_points.first;
        dealii::Point<dim, Number> max_corner = boundary_points.second;
        for (int d = 0; d < dim; ++d)
          {
            min_corner_arborx[d] = static_cast<float>(min_corner[d]);
            max_corner_arborx[d] = static_cast<float>(max_corner[d]);
          }
        bounding_boxes.emplace_back(
          std::make_pair(min_corner_arborx, max_corner_arborx));
      }
  }



  std::size_t
  BoundingBoxIntersect::size() const
  {
    return bounding_boxes.size();
  }



  const dealii::BoundingBox<3, float> &
  BoundingBoxIntersect::get(unsigned int i) const
  {
    return bounding_boxes[i];
  }
} // namespace ArborXWrappers

DEAL_II_NAMESPACE_CLOSE

namespace ArborX
{
  // ------------------- Point Primitives AccessTraits ------------------- //
  template <int dim, typename Number>
  std::size_t
  AccessTraits<std::vector<dealii::Point<dim, Number>>, PrimitivesTag>::size(
    const std::vector<dealii::Point<dim, Number>> &v)
  {
    return v.size();
  }



  template <int dim, typename Number>
  Point
  AccessTraits<std::vector<dealii::Point<dim, Number>>, PrimitivesTag>::get(
    const std::vector<dealii::Point<dim, Number>> &v,
    std::size_t                                    i)
  {
    // ArborX assumes that the point coordinates use float and that the point
    // is 3D
    return {{static_cast<float>(v[i][0]),
             static_cast<float>(v[i][1]),
             dim == 2 ? 0 : static_cast<float>(v[i][2])}};
  }



  // ----------------- BoundingBox Primitives AccessTraits ----------------- //
  template <int dim, typename Number>
  std::size_t
  AccessTraits<std::vector<dealii::BoundingBox<dim, Number>>, PrimitivesTag>::
    size(const std::vector<dealii::BoundingBox<dim, Number>> &v)
  {
    return v.size();
  }



  template <int dim, typename Number>
  Box
  AccessTraits<std::vector<dealii::BoundingBox<dim, Number>>, PrimitivesTag>::
    get(const std::vector<dealii::BoundingBox<dim, Number>> &v, std::size_t i)
  {
    const auto boundary_points                  = v[i].get_boundary_points();
    const dealii::Point<dim, Number> min_corner = boundary_points.first;
    const dealii::Point<dim, Number> max_corner = boundary_points.second;
    // ArborX assumes that the bounding box coordinates use float and that the
    // bounding box is 3D
    return {{static_cast<float>(min_corner[0]),
             static_cast<float>(min_corner[1]),
             dim == 2 ? 0.f : static_cast<float>(min_corner[2])},
            {static_cast<float>(max_corner[0]),
             static_cast<float>(max_corner[1]),
             dim == 2 ? 0.f : static_cast<float>(max_corner[2])}};
  }
} // namespace ArborX

// ----------------------- Instantiations --------------------//
#  include "access_traits.inst"

#endif
