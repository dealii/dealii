// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2017 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#ifndef dealii_base_bounding_box_h
#define dealii_base_bounding_box_h


#include <deal.II/base/config.h>

#include <deal.II/base/exceptions.h>
#include <deal.II/base/point.h>

#include <limits>

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
   * most `spacedim - 2`. For example, in 2d this means that the two rectangles
   * touch at a single point, which must then be a vertex of each box. In 3d,
   * this means that two boxes touch along an edge.
   */
  simple_neighbors = 1,

  /**
   * Attached neighbors: neighbors with an intersection of
   * `dimension > spacedim - 2`. For example, in 2d this means that the two
   * rectangles touch along (parts of) their edges. In 3d, it would mean that
   * two boxes touch along (parts of) their faces.
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
   * or one is inside the other. This is a special case of `attached_neighbors`
   * where the two bounding boxes touch along the entirety of their respective
   * faces, or where they overlap in suitable ways.
   */
  mergeable_neighbors = 3
};

/**
 * A class that represents a box of arbitrary dimension <tt>spacedim</tt> and
 * with sides parallel to the coordinate axes, that is, a region
 * @f[
 * [x_0^L, x_0^U] \times ... \times [x_{spacedim-1}^L, x_{spacedim-1}^U],
 * @f]
 * where $(x_0^L , ..., x_{spacedim-1}^L)$ and $(x_0^U , ..., x_{spacedim-1}^U)$
 * denote the two vertices (bottom left and top right) which are used to
 * represent the box. The quantities $x_k^L$ and $x_k^U$ denote the "lower"
 * and "upper" bounds of values that are within the box for each coordinate
 * direction $k$.
 *
 * Geometrically, a bounding box is thus:
 * - 1d: a segment (represented by its vertices in the proper order)
 * - 2d: a rectangle (represented by the vertices V at bottom left, top right)
 * @code
 * .--------V
 * |        |
 * V--------.
 * @endcode
 *
 * - 3d: a cuboid (in which case the two vertices V follow the convention and
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
 * give a general description of the owners of each portion of the mesh. More
 * generally, bounding boxes are often used to roughly describe a region of
 * space in which an object is contained; if a candidate point is not within
 * the bounding box (a test that is cheap to execute), then it is not necessary
 * to perform an expensive test whether the candidate point is in fact inside
 * the object itself. Bounding boxes are therefore often used as a first,
 * cheap rejection test before more detailed checks. As such, bounding boxes
 * serve many of the same purposes as the
 * [convex hull](https://en.wikipedia.org/wiki/Convex_hull), for which it is
 * also relatively straightforward to compute whether a point is inside or
 * outside, though not quite as cheap as for the bounding box.
 *
 * Taking the cross section of a BoundingBox<spacedim> orthogonal to a given
 * direction gives a box in one dimension lower: BoundingBox<spacedim - 1>.
 * In 3d, the 2 coordinates of the cross section of BoundingBox<3> can be
 * ordered in 2 different ways. That is, if we take the cross section orthogonal
 * to the y direction we could either order a 3d-coordinate into a
 * 2d-coordinate as $(x,z)$ or as $(z,x)$. This class uses the second
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
   * Provide a way to get the dimension of an object without explicit
   * knowledge of it's data type.
   */
  static constexpr unsigned int dimension = spacedim;

  /**
   * Standard constructor. Creates an object that corresponds to an empty box,
   * i.e. a degenerate box with both points being the origin.
   */
  BoundingBox() = default;

  /**
   * Standard copy constructor operator.
   */
  BoundingBox(const BoundingBox<spacedim, Number> &box) = default;

  /**
   * Standard copy assignment operator.
   */
  BoundingBox<spacedim, Number> &
  operator=(const BoundingBox<spacedim, Number> &t) = default;

  /**
   * Standard constructor for an empty box around a point @p point.
   */
  BoundingBox(const Point<spacedim, Number> &point);

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
   */
  bool
  has_overlap_with(
    const BoundingBox<spacedim, Number> &other_bbox,
    const double tolerance = std::numeric_limits<Number>::epsilon()) const;

  /**
   * Check which NeighborType @p other_bbox is to the current object.
   */
  NeighborType
  get_neighbor_type(
    const BoundingBox<spacedim, Number> &other_bbox,
    const double tolerance = std::numeric_limits<Number>::epsilon()) const;

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
  extend(const Number amount);

  /**
   * The same as above with the difference that a new BoundingBox instance is
   * created without changing the current object.
   */
  BoundingBox<spacedim, Number>
  create_extended(const Number amount) const;

  /**
   * Increase (or decrease) each side of the bounding box by the given
   * @p relative_amount.
   *
   * After calling this method, the lower left corner of the bounding box will
   * have each coordinate decreased by @p relative_amount * side_length(direction),
   * and the upper right corner of the bounding box will have each coordinate
   * increased by @p relative_amount * side_length(direction).
   *
   * If you call this method with a negative number, and one of the axes of the
   * original bounding box is smaller than relative_amount *
   * side_length(direction) / 2, the method will trigger an assertion.
   */
  BoundingBox<spacedim, Number>
  create_extended_relative(const Number relative_amount) const;

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
   * @note Calling this method in 1d will result in an exception since
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
   * Returns the signed distance from a @p point orthogonal to the bounds of the
   * box in @p direction. The signed distance is negative for points inside the
   * interval described by the bounds of the rectangle in the respective
   * direction, zero for points on the interval boundary and positive for points
   * outside.
   */
  Number
  signed_distance(const Point<spacedim, Number> &point,
                  const unsigned int             direction) const;

  /**
   * Returns the signed distance from a @p point to the bounds of the box. The
   * signed distance is negative for points inside the rectangle, zero for
   * points on the rectangle and positive for points outside the rectangle.
   */
  Number
  signed_distance(const Point<spacedim, Number> &point) const;

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
   * when going over the dimension. This relationship is, in 2d,
   *
   * | locked in 2D | 1d coordinate | 2d coordinate |
   * |:------------:|:-------------:|:-------------:|
   * |     x0       |      (a)      |   (x0,  a)    |
   * |     x1       |      (a)      |   (a , x1)    |
   *
   * and, in 3d,
   *
   * | locked in 3D | 2d coordinates | 3d coordinates |
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
  const Point<spacedim, Number> &p)
  : BoundingBox({p, p})
{}



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
      for (unsigned int d = 0; d < spacedim; ++d)
        {
          min[d] = std::numeric_limits<Number>::infinity();
          max[d] = -std::numeric_limits<Number>::infinity();
        }

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
BoundingBox<spacedim, Number>::operator==(
  const BoundingBox<spacedim, Number> &box) const
{
  return boundary_points == box.boundary_points;
}



template <int spacedim, typename Number>
inline bool
BoundingBox<spacedim, Number>::operator!=(
  const BoundingBox<spacedim, Number> &box) const
{
  return boundary_points != box.boundary_points;
}



template <int spacedim, typename Number>
inline void
BoundingBox<spacedim, Number>::extend(const Number amount)
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
inline BoundingBox<spacedim, Number>
BoundingBox<spacedim, Number>::create_extended(const Number amount) const
{
  // create and modify copy
  auto bb = *this;
  bb.extend(amount);

  return bb;
}



template <int spacedim, typename Number>
inline BoundingBox<spacedim, Number>
BoundingBox<spacedim, Number>::create_extended_relative(
  const Number relative_amount) const
{
  // create and modify copy
  auto bb = *this;

  for (unsigned int d = 0; d < spacedim; ++d)
    {
      bb.boundary_points.first[d] -= relative_amount * side_length(d);
      bb.boundary_points.second[d] += relative_amount * side_length(d);
      Assert(bb.boundary_points.first[d] <= bb.boundary_points.second[d],
             ExcMessage("Bounding Box can't be shrunk this much: the points' "
                        "order should remain bottom left, top right."));
    }

  return bb;
}



template <int spacedim, typename Number>
template <class Archive>
void
BoundingBox<spacedim, Number>::serialize(Archive &ar,
                                         const unsigned int /*version*/)
{
  // Avoid including boost/serialization/utility.hpp by unpacking the pair
  // ourselves
  ar &boundary_points.first;
  ar &boundary_points.second;
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
