// ---------------------------------------------------------------------
//
// Copyright (C) 2021 - 2022 by the deal.II authors
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
  // ------------------- PointPredicate ------------------- //
  template <int dim, typename Number>
  PointPredicate::PointPredicate(
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
  PointPredicate::size() const
  {
    return points.size();
  }



  const dealii::Point<3, float> &
  PointPredicate::get(unsigned int i) const
  {
    return points[i];
  }



  template <int dim, typename Number>
  PointIntersectPredicate::PointIntersectPredicate(
    const std::vector<dealii::Point<dim, Number>> &points)
    : PointPredicate(points)
  {}



  template <int dim, typename Number>
  PointNearestPredicate::PointNearestPredicate(
    const std::vector<dealii::Point<dim, Number>> &points,
    const unsigned int                             n_nearest_neighbors)
    : PointPredicate(points)
    , n_nearest_neighbors(n_nearest_neighbors)
  {}



  unsigned int
  PointNearestPredicate::get_n_nearest_neighbors() const
  {
    return n_nearest_neighbors;
  }

  // ------------------- BoundingBoxPredicate ------------------- //
  template <int dim, typename Number>
  BoundingBoxPredicate::BoundingBoxPredicate(
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
        for (unsigned int d = 0; d < dim; ++d)
          {
            min_corner_arborx[d] = static_cast<float>(min_corner[d]);
            max_corner_arborx[d] = static_cast<float>(max_corner[d]);
          }
        bounding_boxes.emplace_back(
          std::make_pair(min_corner_arborx, max_corner_arborx));
      }
  }



  std::size_t
  BoundingBoxPredicate::size() const
  {
    return bounding_boxes.size();
  }



  const dealii::BoundingBox<3, float> &
  BoundingBoxPredicate::get(unsigned int i) const
  {
    return bounding_boxes[i];
  }



  template <int dim, typename Number>
  BoundingBoxIntersectPredicate::BoundingBoxIntersectPredicate(
    const std::vector<dealii::BoundingBox<dim, Number>> &bounding_boxes)
    : BoundingBoxPredicate(bounding_boxes)
  {}



  template <int dim, typename Number>
  BoundingBoxNearestPredicate::BoundingBoxNearestPredicate(
    const std::vector<dealii::BoundingBox<dim, Number>> &bounding_boxes,
    const unsigned int                                   n_nearest_neighbors)
    : BoundingBoxPredicate(bounding_boxes)
    , n_nearest_neighbors(n_nearest_neighbors)
  {}



  unsigned int
  BoundingBoxNearestPredicate::get_n_nearest_neighbors() const
  {
    return n_nearest_neighbors;
  }

  // ------------------- SpherePredicate ------------------- //
  template <int dim, typename Number>
  SpherePredicate::SpherePredicate(
    const std::vector<std::pair<dealii::Point<dim, Number>, Number>>
      &dim_spheres)
  {
    static_assert(dim != 1, "dim equal to one is not supported.");

    const unsigned int size = dim_spheres.size();
    spheres.reserve(size);
    for (unsigned int i = 0; i < size; ++i)
      {
        // ArborX assumes that the center coordinates and the radius use float
        // and the sphere is 3D
        spheres.emplace_back(std::make_pair(
          dealii::Point<3, float>(
            static_cast<float>(dim_spheres[i].first[0]),
            static_cast<float>(dim_spheres[i].first[1]),
            dim == 2 ? 0.f : static_cast<float>(dim_spheres[i].first[2])),
          static_cast<float>(dim_spheres[i].second)));
      }
  }



  std::size_t
  SpherePredicate::size() const
  {
    return spheres.size();
  }



  const std::pair<dealii::Point<3, float>, float> &
  SpherePredicate::get(unsigned int i) const
  {
    return spheres[i];
  }



  template <int dim, typename Number>
  SphereIntersectPredicate::SphereIntersectPredicate(
    const std::vector<std::pair<dealii::Point<dim, Number>, Number>> &spheres)
    : SpherePredicate(spheres)
  {}



  template <int dim, typename Number>
  SphereNearestPredicate::SphereNearestPredicate(
    const std::vector<std::pair<dealii::Point<dim, Number>, Number>> &spheres,
    const unsigned int n_nearest_neighbors)
    : SpherePredicate(spheres)
    , n_nearest_neighbors(n_nearest_neighbors)
  {}



  unsigned int
  SphereNearestPredicate::get_n_nearest_neighbors() const
  {
    return n_nearest_neighbors;
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
    return {static_cast<float>(v[i][0]),
            static_cast<float>(v[i][1]),
            dim == 2 ? 0 : static_cast<float>(v[i][2])};
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



  // ---------------------- Sphere Primitives AccessTraits ----------------- //
  template <int dim, typename Number>
  std::size_t
  AccessTraits<std::vector<std::pair<dealii::Point<dim, Number>, Number>>,
               PrimitivesTag>::
    size(const std::vector<std::pair<dealii::Point<dim, Number>, Number>> &v)
  {
    return v.size();
  }



  template <int dim, typename Number>
  Sphere
  AccessTraits<std::vector<std::pair<dealii::Point<dim, Number>, Number>>,
               PrimitivesTag>::
    get(const std::vector<std::pair<dealii::Point<dim, Number>, Number>> &v,
        std::size_t                                                       i)
  {
    // ArborX assumes that the center coordinates and the radius use float and
    // the sphere is 3D
    return {{static_cast<float>(v[i].first[0]),
             static_cast<float>(v[i].first[1]),
             dim == 2 ? 0 : static_cast<float>(v[i].first[2])},
            static_cast<float>(v[i].second)};
  }
} // namespace ArborX

// ----------------------- Instantiations --------------------//
#  include "access_traits.inst"

#endif
