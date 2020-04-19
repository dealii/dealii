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

#ifndef dealii_numerics_rtree_h
#define dealii_numerics_rtree_h

#include <deal.II/base/config.h>

#include <deal.II/base/point.h>

#include <deal.II/boost_adaptors/bounding_box.h>
#include <deal.II/boost_adaptors/point.h>
#include <deal.II/boost_adaptors/segment.h>

DEAL_II_DISABLE_EXTRA_DIAGNOSTICS
#include <boost/geometry/index/rtree.hpp>
#include <boost/geometry/strategies/strategies.hpp>
DEAL_II_ENABLE_EXTRA_DIAGNOSTICS

#include <memory>


DEAL_II_NAMESPACE_OPEN

/**
 * A wrapper for the boost::geometry::index::rtree class, implementing a
 * self-balancing spatial index (the R-tree) capable of storing various types of
 * values, using different balancing algorithms.
 *
 * From [Wikipedia](https://en.wikipedia.org/wiki/R-tree):
 * <blockquote>
 * R-trees are tree data structures used for spatial access methods, i.e., for
 * indexing multi-dimensional information such as geographical coordinates,
 * rectangles or polygons. The R-tree was proposed by Antonin Guttman in 1984
 * and has found significant use in both theoretical and applied contexts. A
 * common real-world usage for an R-tree might be to store spatial objects such
 * as restaurant locations or the polygons that typical maps are made of:
 * streets, buildings, outlines of lakes, coastlines, etc. and then find answers
 * quickly to queries such as "Find all museums within 2 km of my current
 * location", "retrieve all road segments within 2 km of my location" (to
 * display them in a navigation system) or "find the nearest gas station"
 * (although not taking roads into account). The R-tree can also accelerate
 * nearest neighbor search for various distance metrics, including
 * great-circle distance.
 *
 * The key idea of the data structure is to group nearby objects and represent
 * them with their minimum bounding rectangle in the next higher level of the
 * tree; the "R" in R-tree is for rectangle. Since all objects lie within this
 * bounding rectangle, a query that does not intersect the bounding rectangle
 * also cannot intersect any of the contained objects. At the leaf level, each
 * rectangle describes a single object; at higher levels the aggregation of an
 * increasing number of objects. This can also be seen as an increasingly coarse
 * approximation of the data set.
 *
 * The key difficulty of R-tree is to build an efficient tree that on one hand
 * is balanced (so the leaf nodes are at the same height) on the other hand the
 * rectangles do not cover too much empty space and do not overlap too much (so
 * that during search, fewer subtrees need to be processed). For example, the
 * original idea for inserting elements to obtain an efficient tree is to always
 * insert into the subtree that requires least enlargement of its bounding box.
 * Once that page is full, the data is split into two sets that should cover the
 * minimal area each. Most of the research and improvements for R-trees aims at
 * improving the way the tree is built and can be grouped into two objectives:
 * building an efficient tree from scratch (known as bulk-loading) and
 * performing changes on an existing tree (insertion and deletion).
 * </blockquote>
 *
 * An RTree may store any type of @p LeafType as long as it is possible to extract
 * an @p Indexable that the RTree can handle and compare values. An @p Indexable
 * is a type adapted to the Point, BoundingBox or Segment concept, for which
 * distance and equality comparison are implemented. The deal.II Point, Segment,
 * and BoundingBox classes satisfy this requirement, but you can mix in any
 * geometry object that boost::geometry accepts as indexable.
 *
 * In particular, given an @p Indexable type (for example a Point,  a BoundingBox,
 * or a Segment), @p LeafType can by any of @p Indexable, `std::pair<Indexable, T>`,
 * `boost::tuple<Indexable, ...>` or `std::tuple<Indexable, ...>`.
 *
 * The optional argument @p IndexType is used only when adding elements to the
 * tree one by one. If a range insertion is used, then the tree is built using
 * the packing algorithm.
 *
 * Linear, quadratic, and rstar algorithms are available if one wants to
 * construct the tree sequentially. However, none of these is very efficient,
 * and users should use the packing algorithm when possible.
 *
 * The packing algorithm constructs the tree all at once, and may be used when
 * you have all the leaves at your disposal.
 *
 * This class is usually used in combination with one of the two helper
 * functions pack_rtree(), that takes a container or a range of iterators to
 * construct the RTree using the packing algorithm.
 *
 * An example usage is the following:
 *
 * @code
 * std::vector<Point<2>> points = generate_some_points();
 * auto tree = pack_rtree(points.begin(), points.end());
 * // or, equivalently:
 * // auto tree = pack_rtree(points);
 * @endcode
 *
 * The tree is accessed by using [`boost::geometry::index`
 * queries](https://www.boost.org/doc/libs/1_68_0/libs/geometry/doc/html/geometry/spatial_indexes/queries.html).
 * For example, after constructing the tree with the snippet above, one can ask
 * for the closest points to a segment in the following way:
 *
 * @code
 * namespace bgi = boost::geometry::index;
 *
 * Segment<2> segment(Point<2>(0,0), Point<2>(1,1));
 *
 * std::vector<Point<2>> nearest;
 * tree.query(bgi::nearest(segment,3), std::back_inserter(intersection));
 * // Returns the 3 closest points to the Segment defined above.
 * @endcode
 *
 * @author Luca Heltai, 2018.
 */
template <typename LeafType,
          typename IndexType = boost::geometry::index::linear<16>>
using RTree = boost::geometry::index::rtree<LeafType, IndexType>;

/**
 * Construct the correct RTree object by passing an iterator range.
 *
 * Notice that the order of the parameters is the opposite with respect to the
 * RTree class, since we can automatically infer the @p LeafType from the
 * arguments, and we only need to specify the @p IndexType if the default is not
 * adequate.
 *
 * @author Luca Heltai, 2018.
 */
template <typename IndexType = boost::geometry::index::linear<16>,
          typename LeafTypeIterator>
RTree<typename LeafTypeIterator::value_type, IndexType>
pack_rtree(const LeafTypeIterator &begin, const LeafTypeIterator &end);

/**
 * Construct the correct RTree object by passing an STL container type.
 *
 * Notice that the order of the parameters is the opposite with respect to the
 * RTree class, since we can automatically infer the @p LeafType from the
 * arguments, and we only need to specify the @p IndexType if the default is not
 * adequate.
 *
 * @author Luca Heltai, 2018.
 */
template <typename IndexType = boost::geometry::index::linear<16>,
          typename ContainerType>
RTree<typename ContainerType::value_type, IndexType>
pack_rtree(const ContainerType &container);



// Inline and template functions
#ifndef DOXYGEN
template <typename IndexType, typename LeafTypeIterator>
RTree<typename LeafTypeIterator::value_type, IndexType>
pack_rtree(const LeafTypeIterator &begin, const LeafTypeIterator &end)
{
  return RTree<typename LeafTypeIterator::value_type, IndexType>(begin, end);
}



template <typename IndexType, typename ContainerType>
RTree<typename ContainerType::value_type, IndexType>
pack_rtree(const ContainerType &container)
{
  return pack_rtree<IndexType>(container.begin(), container.end());
}

#endif

DEAL_II_NAMESPACE_CLOSE

#endif
