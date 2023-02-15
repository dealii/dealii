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

#ifndef dealii_arborx_bvh_h
#define dealii_arborx_bvh_h

#include <deal.II/base/config.h>

#ifdef DEAL_II_WITH_ARBORX
#  include <deal.II/arborx/access_traits.h>

#  include <ArborX_LinearBVH.hpp>
#  include <Kokkos_Core.hpp>

DEAL_II_NAMESPACE_OPEN

/**
 * This namespace contains wrappers for the ArborX library.
 */
namespace ArborXWrappers
{
  /**
   * This class implements a wrapper around ArborX::BVH. BVH stands for Bounding
   * Volume Hierarchy.
   *
   * From [Wikipedia](https://en.wikipedia.org/wiki/Bounding_Volume_Hierarchy):
   * <blockquote>
   * A bounding volume hierarchy (BVH) is a tree structure on a set of geometric
   * objects. All geometric objects are wrapped in bounding volumes that form
   * the leaf nodes of the tree. These nodes are then grouped as small sets and
   * enclosed within larger bounding volumes. These, in turn, are also grouped
   * and enclosed within other larger bounding volumes in a recursive fashion,
   * eventually resulting in a tree structure with a single bounding volume at
   * the top of the tree. Bounding volume hierarchies are used to support
   * several operations on sets of geometric objects efficiently, such as in
   * collision detection and ray tracing.
   * </blockquote>
   *
   * Because ArborX uses Kokkos, Kokkos needs to be initialized and finalized
   * before using this class.
   */
  class BVH
  {
  public:
    /**
     * Constructor. Use a vector of BoundingBox @p bounding_boxes as primitives.
     */
    template <int dim, typename Number>
    BVH(const std::vector<BoundingBox<dim, Number>> &bounding_boxes);

    /**
     * Constructor. Use a vector of @p points as primitives.
     */
    template <int dim, typename Number>
    BVH(const std::vector<Point<dim, Number>> &points);

    /**
     * Return the indices of those BoundingBox objects that satisfy the @p queries.
     * Because @p queries can contain multiple queries, the function returns a pair
     * of indices and offsets.
     *
     * Below is an example piece of code that does the following: Let us assume
     * that we have a set of bounding boxes for objects -- say, the bounding
     * boxes of each of the cells in a triangulation, or the bounding boxes for
     * each of the parts of a triangulation in a parallel triangulation. We will
     * store those in the `bvh_bounding_boxes` array below.
     *
     * Let us then also assume that we have a set of other bounding boxes, let's
     * say for small objects that are moving around in our domain. We will put
     * these bounding boxes into the `bb_intersect` array. The question we would
     * then like to answer is the following: with which of the BVH bounding
     * box(es) do each of the bb_intersect bounding boxes intersect? In other
     * words, in which cell(s) or partition(s) are the particles?
     *
     * This query is answered by the following piece of code:
     *
     * @code
     * const std::vector<BoundingBox<dim>> query_bounding_boxes = ...
     * ArborXWrappers::BoundingBoxIntersectPredicate
     * bb_intersect(query_bounding_boxes);
     *
     * const std::vector<BoundingBox<dim>> bvh_bounding_boxes = ...
     * ArborxWrappers::BVH bvh(bvh_bounding_boxes);
     *
     * auto [indices, offset] = bvh.query(bb_intersect);
     * @endcode
     *
     * The elements of `bvh_bounding_boxes` that intersect the `j`th BoundingBox
     * of `query_bounding_boxes` are given by:
     *
     * @code
     * std::vector<int> bvh_bounding_box_indices;
     * for (int i = offset[j]; i < offset[j+1]; ++i)
     *   bvh_bounding_box_indices.push_back(indices[i]);
     * @endcode
     *
     * In many other applications, we are interested not only in finding which
     * bounding boxes another bounding box lies in, but in fact which bounding
     * boxes individual points lie in -- say, if instead of objects we have
     * point-like particles moving around. In that case, we would need to answer
     * a query for points, and this can be done as follows:
     *
     * @code
     * const std::vector<Point<dim>> query_points = ...
     * ArborXWrappers::PointIntersectPredicate pt_intersect(query_points);
     *
     * const std::vector<BoundingBox<dim>> bvh_bounding_boxes = ...
     * ArborxWrappers::BVH bvh(bvh_bounding_boxes);
     *
     * auto [indices, offset] = bvh.query(pt_intersect);
     * @endcode
     *
     * As a final example, we want to show how to find the five nearest points
     * of a given set of points. This can done as follows:
     *
     * @code
     * const std::vector<Point<dim>> query_points = ...
     * ArborXWrappers::PointNearestPredicate pt_nearest(query_points, 5);
     *
     * const std::vector<Point<dim>> bvh_points = ...
     * ArborxWrappers::BVH bvh(bvh_points);
     *
     * auto [indices, offset] = bvh.query(pt_nearest);
     * @endcode
     */
    template <typename QueryType>
    std::pair<std::vector<int>, std::vector<int>>
    query(const QueryType &queries);

  private:
    /**
     * Underlying ArborX object.
     */
    ArborX::BVH<Kokkos::HostSpace> bvh;
  };



  template <int dim, typename Number>
  BVH::BVH(const std::vector<BoundingBox<dim, Number>> &bounding_boxes)
    : bvh(Kokkos::DefaultHostExecutionSpace{}, bounding_boxes)
  {}



  template <int dim, typename Number>
  BVH::BVH(const std::vector<Point<dim, Number>> &points)
    : bvh(Kokkos::DefaultHostExecutionSpace{}, points)
  {}



  template <typename QueryType>
  std::pair<std::vector<int>, std::vector<int>>
  BVH::query(const QueryType &queries)
  {
    Kokkos::View<int *, Kokkos::HostSpace> indices("indices", 0);

    Kokkos::View<int *, Kokkos::HostSpace> offset("offset", 0);
    ArborX::query(
      bvh, Kokkos::DefaultHostExecutionSpace{}, queries, indices, offset);
    std::vector<int> indices_vector;
    indices_vector.insert(indices_vector.begin(),
                          indices.data(),
                          indices.data() + indices.extent(0));
    std::vector<int> offset_vector;
    offset_vector.insert(offset_vector.begin(),
                         offset.data(),
                         offset.data() + offset.extent(0));

    return {indices_vector, offset_vector};
  }
} // namespace ArborXWrappers

DEAL_II_NAMESPACE_CLOSE

#endif
#endif
