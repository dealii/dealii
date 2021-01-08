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

#ifndef dealii_arborx_access_traits_h
#define dealii_arborx_access_traits_h

#include <deal.II/base/config.h>

#ifdef DEAL_II_WITH_ARBORX
#  include <deal.II/base/bounding_box.h>

#  include <ArborX.hpp>


DEAL_II_NAMESPACE_OPEN

namespace ArborXWrappers
{
  /**
   *  This class is used by ArborXWrappers::BVH to determine which points
   *  intersect the bounding boxes used to build the ArborXWrappers::BVH.
   */
  class PointIntersect
  {
  public:
    /**
     * Constructor. @p dim_points is a list of points which we are interested in
     * knowing if they intersect ArborXWrappers::BVH bounding boxes.
     */
    template <int dim, typename Number>
    PointIntersect(const std::vector<dealii::Point<dim, Number>> &dim_points);

    /**
     * Number of points stored in the structure.
     */
    std::size_t
    size() const;

    /**
     * Return the `i`th Point stored in the object.
     */
    const dealii::Point<3, float> &
    get(unsigned int i) const;

  private:
    std::vector<dealii::Point<3, float>> points;
  };



  /**
   *  This class is used by ArborXWrappers::BVH to determine which bounding
   *  boxes intersect the bounding boxes used to build the ArborXWrappers::BVH.
   */
  class BoundingBoxIntersect
  {
  public:
    /**
     * Constructor. @p bb is a list of bounding boxes which we are interested in
     * knowing if they intersect ArborXWrappers::BVH bounding boxes.
     */
    template <int dim, typename Number>
    BoundingBoxIntersect(
      const std::vector<dealii::BoundingBox<dim, Number>> &bb);

    /**
     * Number of bounding boxes stored in the structure.
     */
    std::size_t
    size() const;

    /**
     * Return the `i`th BoundingBox stored in the object.
     */
    const dealii::BoundingBox<3, float> &
    get(unsigned int i) const;

  private:
    std::vector<dealii::BoundingBox<3, float>> bounding_boxes;
  };
} // namespace ArborXWrappers

DEAL_II_NAMESPACE_CLOSE

/**
 * This namespace contains the implementation of AccessTraits used by ArborX.
 */
namespace ArborX
{
  /**
   * This struct allows ArborX to use std::vector<dealii::Point> as
   * primitive.
   */
  template <int dim, typename Number>
  struct AccessTraits<std::vector<dealii::Point<dim, Number>>, PrimitivesTag>
  {
    using memory_space = Kokkos::HostSpace;

    /**
     * Return the size of the vector @p v.
     */
    static std::size_t
    size(const std::vector<dealii::Point<dim, Number>> &v);

    /**
     * Return an ArborX::Point from the dealii::Point `v[i]`.
     */
    static Point
    get(const std::vector<dealii::Point<dim, Number>> &v, std::size_t i);
  };



  /**
   * This struct allows ArborX to use std::vector<dealii::BoundingBox> as
   * primitive.
   */
  template <int dim, typename Number>
  struct AccessTraits<std::vector<dealii::BoundingBox<dim, Number>>,
                      PrimitivesTag>
  {
    using memory_space = Kokkos::HostSpace;

    /**
     * Return the size of the vector @p v.
     */
    static std::size_t
    size(const std::vector<dealii::BoundingBox<dim, Number>> &v);

    /**
     * Return an ArborX::Box from the dealii::BoundingBox `v[i]`.
     */
    static Box
    get(const std::vector<dealii::BoundingBox<dim, Number>> &v, std::size_t i);
  };



  /**
   * This struct allows ArborX to use PointIntersect in a predicate.
   */
  template <>
  struct AccessTraits<dealii::ArborXWrappers::PointIntersect, PredicatesTag>
  {
    using memory_space = Kokkos::HostSpace;

    /**
     * Number of Point stored in @p pt_intersect.
     */
    static std::size_t
    size(const dealii::ArborXWrappers::PointIntersect &pt_intersect);

    /**
     * Return an Arbox::intersects(ArborX::Point) object constructed from the
     * `i`th dealii::Point stored in @p pt_intersect.
     */
    static auto
    get(const dealii::ArborXWrappers::PointIntersect &pt_intersect,
        std::size_t                                   i);
  };



  /**
   * This struct allows ArborX to use BoundingBoxIntersect in a predicate.
   */
  template <>
  struct AccessTraits<dealii::ArborXWrappers::BoundingBoxIntersect,
                      PredicatesTag>
  {
    using memory_space = Kokkos::HostSpace;

    /**
     * Number of BoundingBox stored in @p bb_intersect.
     */
    static std::size_t
    size(const dealii::ArborXWrappers::BoundingBoxIntersect &bb_intersect);

    /**
     * Return an Arbox::intersects(ArborX::Box) object constructed from the
     * `i`th dealii::BoundingBox stored in @p bb_intersect.
     */
    static auto
    get(const dealii::ArborXWrappers::BoundingBoxIntersect &bb_intersect,
        std::size_t                                         i);
  };

  // ------------------------------- Inline ----------------------------------//

  // The implementation of AccessTraits<..., PredicatesTag> needs to be in the
  // header file otherwise the return type of auto get() cannot be determined.
  // We use auto because ArborX does not expose the type of intersects

  inline std::size_t
  AccessTraits<dealii::ArborXWrappers::PointIntersect, PredicatesTag>::size(
    const dealii::ArborXWrappers::PointIntersect &pt_intersect)
  {
    return pt_intersect.size();
  }



  inline auto
  AccessTraits<dealii::ArborXWrappers::PointIntersect, PredicatesTag>::get(
    const dealii::ArborXWrappers::PointIntersect &pt_intersect,
    std::size_t                                   i)
  {
    const auto dealii_point = pt_intersect.get(i);
    return intersects(Point{dealii_point[0], dealii_point[1], dealii_point[2]});
  }


  inline std::size_t
  AccessTraits<dealii::ArborXWrappers::BoundingBoxIntersect, PredicatesTag>::
    size(const dealii::ArborXWrappers::BoundingBoxIntersect &bb_intersect)
  {
    return bb_intersect.size();
  }



  inline auto
  AccessTraits<dealii::ArborXWrappers::BoundingBoxIntersect, PredicatesTag>::
    get(const dealii::ArborXWrappers::BoundingBoxIntersect &bb_intersect,
        std::size_t                                         i)
  {
    const auto boundary_points = bb_intersect.get(i).get_boundary_points();
    const dealii::Point<3, float> min_corner = boundary_points.first;
    const dealii::Point<3, float> max_corner = boundary_points.second;

    return intersects(Box{{min_corner[0], min_corner[1], min_corner[2]},
                          {max_corner[0], max_corner[1], max_corner[2]}});
  }
} // namespace ArborX

#endif

#endif
