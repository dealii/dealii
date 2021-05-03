// ---------------------------------------------------------------------
//
// Copyright (C) 2017 - 2020 by the deal.II authors
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
   * Not neighbors: the intersection is empty.
   */
  not_neighbors = 0,

  /**
   * Simple neighbors: the boxes intersect with an intersection of dimension at
   * most `spacedim - 2`. For example, in 2d this means that the two boxes
   * touch at one corner of the each box.
   */
  simple_neighbors = 1,

  /**
   * Attached neighbors: neighbors with an intersection of
   * `dimension > spacedim - 2`. For example, in 2d this means that the two
   * boxes touch along an edge.
   */
  attached_neighbors = 2,

  /**
   * Mergeable neighbors: neighbors which can be expressed with a single
   * BoundingBox, e.g.
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
 * A class that represents a box of arbitrary dimension <tt>spacedim</tt> and
 * with sides parallel to the coordinate axes, that is, a region
 *
 * @f[
 * [x_0^L, x_0^U] \times ... \times [x_{spacedim-1}^L, x_{spacedim-1}^U],
 * @f]
 *
 * where $(x_0^L , ..., x_{spacedim-1}^L) and $(x_0^U , ..., x_{spacedim-1}^U)
 * denote the two vertices (bottom left and top right) which are used to
 * represent the box.
 *
 * Geometrically, a bounding box is thus:
 * - 1D: a segment (represented by its vertices in the proper order)
 * - 2D: a rectangle (represented by the vertices V at bottom left, top right)
 * @code
 * .--------V
 * |        |
 * V--------.
 * @endcode
 *
 * - 3D: a cuboid (in which case the two vertices V follow the convention and
 * are not owned by the same face)
 * @code
 *   .------V
 *  /      /|
 * .------. |
 * |      | /
 * |      |/
 * V------.
 * @endcode
 *
 * Bounding boxes are, for example, useful in parallel distributed meshes to
 * give a general description of the owners of each portion of the mesh.
 *
 * Taking the cross section of a BoundingBox<spacedim> orthogonal to a given
 * direction gives a box in one dimension lower: BoundingBox<spacedim - 1>.
 * In 3D, the 2 coordinates of the cross section of BoundingBox<3> can be
 * ordered in 2 different ways. That is, if we take the cross section orthogonal
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
   * Return a const reference to the boundary_points
   */
  const std::pair<Point<spacedim, Number>, Point<spacedim, Number>> &
  get_boundary_points() const;

  /**
   * Test for equality.
   */
  bool
  operator==(const BoundingBox<spacedim, Number> &box) const;

  /**
   * Test for inequality.
   */
  bool
  operator!=(const BoundingBox<spacedim, Number> &box) const;

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
   * Return true if the point is inside the Bounding Box, false otherwise. The
   * parameter @p tolerance is a factor by which the bounding box is enlarged
   * relative to the dimensions of the bounding box in order to determine in a
   * numerically robust way whether the point is inside.
   */
  bool
  point_inside(
    const Point<spacedim, Number> &p,
    const double tolerance = std::numeric_limits<Number>::epsilon()) const;

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
   * Returns the point in the center of the box.
   */
  Point<spacedim, Number>
  center() const;

  /**
   * Returns the side length of the box in @p direction.
   */
  Number
  side_length(const unsigned int direction) const;

  /**
   * Return the lower bound of the box in @p direction.
   */
  Number
  lower_bound(const unsigned int direction) const;

  /**
   * Return the upper bound of the box in @p direction.
   */
  Number
  upper_bound(const unsigned int direction) const;

  /**
   * Return the bounds of the box in @p direction, as a one-dimensional box.
   */
  BoundingBox<1, Number>
  bounds(const unsigned int direction) const;

  /**
   * Returns the indexth vertex of the box. Vertex is meant in the same way as
   * for a cell, so that @p index $\in [0, 2^{\text{dim}} - 1]$.
   */
  Point<spacedim, Number>
  vertex(const unsigned int index) const;

  /**
   * Returns the indexth child of the box. Child is meant in the same way as for
   * a cell.
   */
  BoundingBox<spacedim, Number>
  child(const unsigned int index) const;

  /**
   * Returns the cross section of the box orthogonal to @p direction.
   * This is a box in one dimension lower.
   *
   * @note Calling this method in 1D will result in an exception since
   * <code>BoundingBox&lt;0&gt;</code> is not implemented.
   */
  BoundingBox<spacedim - 1, Number>
  cross_section(const unsigned int direction) const;

  /**
   * Apply the affine transformation that transforms this BoundingBox to a unit
   * BoundingBox object.
   *
   * If $B$ is this bounding box, and $\hat{B}$ is the unit bounding box,
   * compute the affine mapping that satisfies $G(B) = \hat{B}$ and apply it to
   * @p point.
   */
  Point<spacedim, Number>
  real_to_unit(const Point<spacedim, Number> &point) const;

  /**
   * Apply the affine transformation that transforms the unit BoundingBox object
   * to this object.
   *
   * If $B$ is this bounding box, and $\hat{B}$ is the unit bounding box,
   * compute the affine mapping that satisfies $F(\hat{B}) = B$ and apply it to
   * @p point.
   */
  Point<spacedim, Number>
  unit_to_real(const Point<spacedim, Number> &point) const;

  /**
   * Write or read the data of this object to or from a stream for the
   * purpose of serialization using the [BOOST serialization
   * library](https://www.boost.org/doc/libs/1_74_0/libs/serialization/doc/index.html).
   */
  template <class Archive>
  void
  serialize(Archive &ar, const unsigned int version);

private:
  std::pair<Point<spacedim, Number>, Point<spacedim, Number>> boundary_points;
};

/**
 * Specialization of BoundingBox for spacedim 0. This class exists to enable
 * dimension-independent programming but unconditionally throws an exception
 * in its constructor.
 */
template <typename Number>
class BoundingBox<0, Number>
{
public:
  /**
   * Default constructor. Throws an exception.
   */
  BoundingBox();

  /**
   * Equivalent two-point constructor. Throws an exception.
   */
  BoundingBox(const std::pair<Point<0, Number>, Point<0, Number>> &);

  /**
   * Equivalent container constructor. Throws an exception.
   */
  template <class Container>
  BoundingBox(const Container &);
};


/**
 * Returns the unit box $[0,1]^\text{dim}$.
 *
 * @relates BoundingBox
 */
template <int dim, typename Number = double>
BoundingBox<dim, Number>
create_unit_bounding_box();


namespace internal
{
  /**
   * This function defines a convention for how coordinates in dim
   * dimensions should translate to the coordinates in dim + 1 dimensions,
   * when one of the coordinates in dim + 1 dimensions is locked to a given
   * value.
   *
   * The convention is the following: Starting from the locked coordinate we
   * store the lower dimensional coordinates consecutively and wrap around
   * when going over the dimension. This relationship is, in 2D,
   *
   * | locked in 2D | 1D coordinate | 2D coordinate |
   * |:------------:|:-------------:|:-------------:|
   * |     x0       |      (a)      |   (x0,  a)    |
   * |     x1       |      (a)      |   (a , x1)    |
   *
   * and, in 3D,
   *
   * | locked in 3D | 2D coordinates | 3D coordinates |
   * |:-------------|:--------------:|:--------------:|
   * |     x0       |    (a, b)      | (x0,  a,  b)   |
   * |     x1       |    (a, b)      | ( b, x1,  a)   |
   * |     x2       |    (a, b)      | ( a,  b, x2)   |
   *
   * Given a locked coordinate, this function maps a coordinate index in dim
   * dimension to a coordinate index in dim + 1 dimensions.
   *
   * @param locked_coordinate should be in the range [0, dim+1).
   * @param coordinate_in_dim should be in the range [0, dim).
   * @return A coordinate index in the range [0, dim+1)
   *
   * @relates BoundingBox
   */
  template <int dim>
  inline int
  coordinate_to_one_dim_higher(const int locked_coordinate,
                               const int coordinate_in_dim)
  {
    AssertIndexRange(locked_coordinate, dim + 1);
    AssertIndexRange(coordinate_in_dim, dim);
    return (locked_coordinate + coordinate_in_dim + 1) % (dim + 1);
  }

} // namespace internal

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
  // Use the default constructor in case points is empty instead of setting
  // things to +oo and -oo
  if (points.size() > 0)
    {
      auto &min = boundary_points.first;
      auto &max = boundary_points.second;
      std::fill(min.begin_raw(),
                min.end_raw(),
                std::numeric_limits<Number>::infinity());
      std::fill(max.begin_raw(),
                max.end_raw(),
                -std::numeric_limits<Number>::infinity());

      for (const Point<spacedim, Number> &point : points)
        for (unsigned int d = 0; d < spacedim; ++d)
          {
            min[d] = std::min(min[d], point[d]);
            max[d] = std::max(max[d], point[d]);
          }
    }
}



template <int spacedim, typename Number>
inline std::pair<Point<spacedim, Number>, Point<spacedim, Number>> &
BoundingBox<spacedim, Number>::get_boundary_points()
{
  return this->boundary_points;
}



template <int spacedim, typename Number>
inline const std::pair<Point<spacedim, Number>, Point<spacedim, Number>> &
BoundingBox<spacedim, Number>::get_boundary_points() const
{
  return this->boundary_points;
}



template <int spacedim, typename Number>
inline bool
BoundingBox<spacedim, Number>::
operator==(const BoundingBox<spacedim, Number> &box) const
{
  return boundary_points == box.boundary_points;
}



template <int spacedim, typename Number>
inline bool
BoundingBox<spacedim, Number>::
operator!=(const BoundingBox<spacedim, Number> &box) const
{
  return boundary_points != box.boundary_points;
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



template <typename Number>
inline BoundingBox<0, Number>::BoundingBox()
{
  AssertThrow(false, ExcImpossibleInDim(0));
}



template <typename Number>
inline BoundingBox<0, Number>::BoundingBox(
  const std::pair<Point<0, Number>, Point<0, Number>> &)
{
  AssertThrow(false, ExcImpossibleInDim(0));
}



template <typename Number>
template <class Container>
inline BoundingBox<0, Number>::BoundingBox(const Container &)
{
  AssertThrow(false, ExcImpossibleInDim(0));
}



#endif // DOXYGEN
DEAL_II_NAMESPACE_CLOSE

#endif
