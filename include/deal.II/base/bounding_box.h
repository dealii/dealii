// ---------------------------------------------------------------------
//
// Copyright (C) 2017 - 2018 by the deal.II authors
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

#ifndef dealii_base_bounding_box_h
#define dealii_base_bounding_box_h


#include <deal.II/base/config.h>

#include <deal.II/base/exceptions.h>
#include <deal.II/base/point.h>
#include <deal.II/base/utilities.h>

DEAL_II_NAMESPACE_OPEN

/**
 * The enumerator NeighborType describes the neighboring relation between
 * two bounding boxes.
 */
enum class NeighborType
{
  /**
   * not neighbours: the intersection is empty
   */
  not_neighbors = 0,

  /**
   * simple neighbors: the boxes intersect with an intersection of dimension at
   * most spacedim - 2
   */
  simple_neighbors = 1,

  /**
   * attached neighbors: neighbors with an intersection of dimension > spacedim
   * - 2
   */
  attached_neighbors = 2,

  /**
   * mergeable neighbors: neighbors which can be expressed with a single
   * Bounding Box, e.g.
   *  @code
   *  .--V--W    .-----V
   *  |  |  | =  |     |
   *  V--W--.    V-----.
   *  @endcode
   * or one is inside the other
   */
  mergeable_neighbors = 3
};

/**
 * A class that represents a bounding box in a space with arbitrary dimension
 * <tt>spacedim</tt>.
 *
 * Objects of this class are used to represent bounding boxes. They are,
 * among other uses, useful in parallel distributed meshes to give a general
 * description of the owners of each portion of the mesh.
 *
 * Bounding boxes are represented by two vertices (bottom left and top right).
 * Geometrically, a bounding box is:
 * - 1 d: a segment (represented by its vertices in the proper order)
 * - 2 d: a rectangle (represented by the vertices V at bottom left, top right)
 * @code
 * .--------V
 * |        |
 * V--------.
 * @endcode
 *
 * - 3 d: a cuboid (in which case the two vertices V follow the convention and
 * are not owned by the same face)
 * @code
 *   .------V
 *  /      /|
 * .------. |
 * |      | /
 * |      |/
 * V------.
 * @endcode
 * Notice that the sides are always parallel to the respective axis.
 *
 * @author Giovanni Alzetta, 2017.
 */
template <int spacedim, typename Number = double>
class BoundingBox
{
public:
  /**
   * Standard constructor. Creates an object that corresponds to an empty box,
   * i.e. a degenerate box with both points being the origin.
   */
  BoundingBox() = default;

  /**
   * Standard constructor for non-empty boxes: it uses a pair of points
   * which describe the box: one for the bottom and one for the top
   * corner.
   */
  BoundingBox(const std::pair<Point<spacedim, Number>, Point<spacedim, Number>>
                &boundary_points);

  /**
   * Construct the bounding box that encloses all the points in the given
   * container.
   *
   * The constructor supports any Container that provides begin() and end()
   * iterators to Point<spacedim, Number> elements.
   */
  template <class Container>
  BoundingBox(const Container &points);

  /**
   * Return a reference to the boundary_points
   */
  std::pair<Point<spacedim, Number>, Point<spacedim, Number>> &
  get_boundary_points();

  /**
   * Check if the current object and @p other_bbox are neighbors, i.e. if the boxes
   * have dimension spacedim, check if their intersection is non empty.
   *
   * Return an enumerator of type NeighborType.
   */
  NeighborType
  get_neighbor_type(const BoundingBox<spacedim, Number> &other_bbox) const;

  /**
   * Enlarge the current object so that it contains @p other_bbox .
   * If the current object already contains @p other_bbox then it is not changed
   * by this function.
   */
  void
  merge_with(const BoundingBox<spacedim, Number> &other_bbox);

  /**
   * Return true if the point is inside the Bounding Box, false otherwise.
   */
  bool
  point_inside(const Point<spacedim, Number> &p) const;

  /**
   * Increase (or decrease) the size of the bounding box by the given amount.
   * After calling this method, the lower left corner of the bounding box will
   * have each coordinate decreased by @p amount, and the upper right corner
   * of the bounding box will have each coordinate increased by @p amount.
   *
   * If you call this method with a negative number, and one of the axes of the
   * original bounding box is smaller than amount/2, the method will trigger
   * an assertion.
   */
  void
  extend(const Number &amount);

  /**
   * Compute the volume (i.e. the dim-dimensional measure) of the BoundingBox.
   */
  double
  volume() const;

  /**
   * Boost serialization function
   */
  template <class Archive>
  void
  serialize(Archive &ar, const unsigned int version);

private:
  std::pair<Point<spacedim, Number>, Point<spacedim, Number>> boundary_points;
};

/*------------------------ Inline functions: BoundingBox --------------------*/

#ifndef DOXYGEN


template <int spacedim, typename Number>
inline BoundingBox<spacedim, Number>::BoundingBox(
  const std::pair<Point<spacedim, Number>, Point<spacedim, Number>>
    &boundary_points)
{
  // We check the Bounding Box is not degenerate
  for (unsigned int i = 0; i < spacedim; ++i)
    Assert(boundary_points.first[i] <= boundary_points.second[i],
           ExcMessage("Bounding Box can't be created: the points' "
                      "order should be bottom left, top right!"));

  this->boundary_points = boundary_points;
}



template <int spacedim, typename Number>
template <class Container>
inline BoundingBox<spacedim, Number>::BoundingBox(const Container &points)
{
  boost::geometry::envelope(
    boost::geometry::model::multi_point<Point<spacedim, Number>>(points.begin(),
                                                                 points.end()),
    *this);
}



template <int spacedim, typename Number>
inline void
BoundingBox<spacedim, Number>::extend(const Number &amount)
{
  for (unsigned int d = 0; d < spacedim; ++d)
    {
      boundary_points.first[d] -= amount;
      boundary_points.second[d] += amount;
      Assert(boundary_points.first[d] <= boundary_points.second[d],
             ExcMessage("Bounding Box can't be shrunk this much: the points' "
                        "order should remain bottom left, top right."));
    }
}


template <int spacedim, typename Number>
template <class Archive>
void
BoundingBox<spacedim, Number>::serialize(Archive &ar,
                                         const unsigned int /*version*/)
{
  ar &boundary_points;
}

#endif // DOXYGEN
DEAL_II_NAMESPACE_CLOSE

#endif
