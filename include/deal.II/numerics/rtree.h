// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2018 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#ifndef dealii_numerics_rtree_h
#define dealii_numerics_rtree_h

#include <deal.II/base/config.h>

#include <deal.II/base/point.h>
#include <deal.II/base/std_cxx20/iota_view.h>

#include <deal.II/boost_adaptors/bounding_box.h>
#include <deal.II/boost_adaptors/point.h>
#include <deal.II/boost_adaptors/segment.h>

DEAL_II_DISABLE_EXTRA_DIAGNOSTICS
#include <boost/geometry/algorithms/distance.hpp>
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
 * std::vector<Point<2>> intersection;
 * tree.query(bgi::nearest(segment,3), std::back_inserter(intersection));
 * // Returns the 3 closest points to the Segment defined above.
 * @endcode
 *
 * In general, a tree does not need to store the actual objects, as long as it
 * knows how to access a const reference to an indexable type. This is
 * achieved by passing the optional template argument @p IndexableGetter, that
 * extracts a const reference to one of the possible indexable types given an
 * object of type @p LeafType. As an example, one may store points in a container,
 * and only create a tree of the indices within the container, using the
 * IndexableGetterFromIndices class defined below, and the function
 * pack_rtree_of_indices().
 */
template <typename LeafType,
          typename IndexType = boost::geometry::index::linear<16>,
          typename IndexableGetter =
            boost::geometry::index::indexable<LeafType>>
using RTree =
  boost::geometry::index::rtree<LeafType, IndexType, IndexableGetter>;

/**
 * Construct the correct RTree object by passing an iterator range.
 *
 * Notice that the order of the parameters is the opposite with respect to the
 * RTree class, since we can automatically infer the @p LeafType from the
 * arguments.
 */
template <typename IndexType = boost::geometry::index::linear<16>,
          typename LeafTypeIterator,
          typename IndexableGetter = boost::geometry::index::indexable<
            typename LeafTypeIterator::value_type>>
RTree<typename LeafTypeIterator::value_type, IndexType, IndexableGetter>
pack_rtree(const LeafTypeIterator &begin, const LeafTypeIterator &end);

/**
 * Construct an RTree object by passing an STL container type.
 * This function is used in step-70.
 *
 * Notice that the order of the template parameters is the opposite with respect
 * to the
 * RTree class, since we can automatically infer the @p LeafType from the
 * arguments, and we only need to specify the @p IndexType if the default is not
 * adequate.
 */
template <typename IndexType = boost::geometry::index::linear<16>,
          typename ContainerType,
          typename IndexableGetter = boost::geometry::index::indexable<
            typename ContainerType::value_type>>
RTree<typename ContainerType::value_type, IndexType, IndexableGetter>
pack_rtree(const ContainerType &container);

/**
 * A class that may be used as an @p IndexableGetter template argument for an
 * RTree that stores indices to entries in a @p Container type.
 *
 * This class may be used as a proxy to extract an indexable type from
 * compatible containers. For example:
 * @code
 * std::vector<std::pair<Point<dim>, double> > point_temperature = fill();
 * IndexableGetterFromIndices<decltype(point_temperature)>
 *    getter(point_temperature);
 *
 * const Point<dim> &p = getter(i); // returns point_temperature[i].first;
 * @endcode
 *
 * This class is used by the pack_rtree_of_indices() function to construct an
 * RTree where the leaves are indices pointing to the entries of the container
 * passed to this class.
 */
template <typename Container>
class IndexableGetterFromIndices
{
public:
  /**
   * An alias for the boost type that is used to extract a Point, Segment, or
   * BoundingBox from compatible types (pairs, tuples, etc.).
   */
  using IndexableGetter =
    typename boost::geometry::index::indexable<typename Container::value_type>;

  /**
   * An alias to the actual geometrical type.
   */
  using result_type = typename IndexableGetter::result_type;

  /**
   * An alias to the index type.
   */
  using size_t = typename Container::size_type;

  /**
   * Constructor. Store a const reference to a container.
   */
  explicit IndexableGetterFromIndices(const Container &c)
    : container(c)
  {}

  /**
   * Implements the @p IndexableGetter requirements of the rtree class.
   */
  result_type
  operator()(size_t i) const
  {
    return getter(container[i]);
  }

private:
  /**
   * A const reference to the container.
   */
  const Container &container;

  /**
   * An instantiation of the getter that allows easy translation from the
   * container @p value_type and the actual indexable type.
   */
  IndexableGetter getter;
};

/**
 * Construct a RTree object that stores the indices of an existing container of
 * indexable types. The only requirement on the container is that it supports
 * operator[] for any index between 0 and the size of the container (i.e., a
 * std::vector, or an std::array will do, however an std::map won't).
 *
 * Differently from the object created by the pack_rtree() function, in this
 * case we don't store the actual geometrical types, but just their indices:
 *
 * @code
 * namespace bgi = boost::geometry::index;
 * std::vector<Point<dim>> some_points = fill();
 * auto tree = pack_rtree(points); // the tree contains a copy of the points
 * auto index_tree = pack_rtree_of_indices(points); // the tree contains only
 *                                                  // the indices of the
 *                                                  // points
 * BoundingBox<dim> box = build_a_box();
 *
 * for(const auto &p: tree       | bgi::adaptors::queried(bgi::intersects(box)))
 *   std::cout << "Point p: " << p << std::endl;
 *
 * for(const auto &i: index_tree | bgi::adaptors::queried(bgi::intersects(box)))
 *   std::cout << "Point p: " << some_points[i] << std::endl;
 * @endcode
 *
 * The leaves stored in the tree are the indices of the corresponding entry in
 * the container. A reference to the external container is stored internally,
 * but keep in mind that if you change the container, you should rebuild the
 * tree.
 *
 * @warning This function does not work on Windows. As an alternative,
 * build a tree with pairs of point and id. For instance to
 * create an intersection between `BoundingBox<dim> box` and
 * `std::vector<BoundingBox<dim>> boxes`, don't write:
 * @code
 * const auto tree = pack_rtree_of_indices(boxes);
 *
 * for (const auto &i : tree | bgi::adaptors::queried(bgi::intersects(box)))
 *   std::cout << i << " " << boxes[i] << std::endl;
 * @endcode
 * but instead:
 * @code
 * std::vector<std::pair<BoundingBox<dim>, unsigned int>> boxes_and_ids;
 *
 * for(unsigned int i = 0; i < boxes.size(); ++i)
 *   boxes_and_ids.emplace_back(boxes[i], i);
 *
 * const auto tree = pack_rtree(boxes_and_ids);
 *
 * for (const auto &i : tree | bgi::adaptors::queried(bgi::intersects(box)))
 *   std::cout << i->second << " " << boxes[i->second] << std::endl;
 * @endcode
 *
 */
template <typename IndexType = boost::geometry::index::linear<16>,
          typename ContainerType>
RTree<typename ContainerType::size_type,
      IndexType,
      IndexableGetterFromIndices<ContainerType>>
pack_rtree_of_indices(const ContainerType &container);

/**
 * Helper structure that allows one to extract a level from an RTree as a vector
 * of BoundingBox objects.
 *
 * This structure implements a boost::geometry::index::detail::rtree::visitor
 * object, allowing one to visit any existing RTree object, and return the
 * vector of bounding boxes associated to a specific target level of the tree.
 *
 * Although possible, direct usage of this structure is cumbersome. The
 * suggested usage of this class is through the helper function
 * extract_rtree_level().
 */
template <typename Value,
          typename Options,
          typename Translator,
          typename Box,
          typename Allocators>
struct ExtractLevelVisitor
  : public boost::geometry::index::detail::rtree::visitor<
      Value,
      typename Options::parameters_type,
      Box,
      Allocators,
      typename Options::node_tag,
      true>::type
{
  /**
   * Construct a vector @p boxes of BoundingBox objects corresponding to the
   * @p target_level of the tree.
   */
  inline ExtractLevelVisitor(
    const Translator  &translator,
    const unsigned int target_level,
    std::vector<BoundingBox<boost::geometry::dimension<Box>::value>> &boxes);

  /**
   * An alias that identifies an InternalNode of the tree.
   */
  using InternalNode =
    typename boost::geometry::index::detail::rtree::internal_node<
      Value,
      typename Options::parameters_type,
      Box,
      Allocators,
      typename Options::node_tag>::type;

  /**
   * An alias that identifies a Leaf of the tree.
   */
  using Leaf = typename boost::geometry::index::detail::rtree::leaf<
    Value,
    typename Options::parameters_type,
    Box,
    Allocators,
    typename Options::node_tag>::type;

  /**
   * Implements the visitor interface for InternalNode objects. If the node
   * belongs to the @p target_level, then fill the bounding box vector.
   */
  inline void
  operator()(const InternalNode &node);

  /**
   * Implements the visitor interface for Leaf objects.
   */
  inline void
  operator()(const Leaf &);

  /**
   * Translator interface, required by the boost implementation of the rtree.
   */
  const Translator &translator;

  /**
   * Store the level we are currently visiting.
   */
  std::size_t level;

  /**
   * The level we want to extract from the RTree object.
   */
  const std::size_t target_level;

  /**
   * A reference to the input vector of BoundingBox objects.
   */
  std::vector<BoundingBox<boost::geometry::dimension<Box>::value>> &boxes;
};

/**
 * Given a RTree object @p rtree, and a target level @p level, return a vector
 * of BoundingBox objects containing all the bounding boxes that make the given
 * @p level of the @p rtree. This function is a convenient wrapper around the
 * ExtractLevelVisitor class. It is used in step-70.
 *
 * Since an RTree object is a balanced tree, you can expect each entry of the
 * resulting vector to contain roughly the same number of children, and
 * ultimately, the same number of leaf objects. If you request for a level that
 * is not present in the RTree, the last level is returned.
 *
 * A typical usage of this function is in the context of
 * parallel::distributed::Triangulation objects, where one would like to
 * construct a rough representation of the area which is covered by the locally
 * owned cells of the active process, and exchange this information with other
 * processes. The finest level of information is given by the leaves, which in
 * this context would be the collection of all the bounding boxes associated
 * to the locally owned cells of the triangulation. Exchanging this information
 * with all participating processes would defeat the purpose of parallel
 * computations. If however one constructs an RTree containing these bounding
 * boxes (for example, by calling
 * GridTools::Cache::get_cell_bounding_boxes_rtree()), and then extracts one of
 * the first levels of the RTree, only a handful of BoundingBox objects would be
 * returned, allowing the user to have a very efficient description of the
 * geometry of the domain, and of its distribution among processes.
 *
 * An example usage is given by the following snippet of code:
 * @code
 * parallel::distributed::Triangulation<2> tria(MPI_COMM_WORLD);
 * GridGenerator::hyper_ball(tria);
 * tria.refine_global(4);
 *
 * std::vector<BoundingBox<2>> all_boxes(tria.n_locally_owned_active_cells());
 * unsigned int                i = 0;
 * for (const auto &cell : tria.active_cell_iterators())
 *   if (cell->is_locally_owned())
 *     all_boxes[i++] = cell->bounding_box();
 *
 * const auto tree  = pack_rtree(all_boxes);
 * const auto boxes = extract_rtree_level(tree, 1);
 * @endcode
 *
 * When run on three processes, the complete set of the BoundingBox objects
 * surrounding only the locally owned cells and the second level of the rtree
 * constructed with those boxes would look like in the following pictures (one
 * image per process):
 *
 * @image html rtree-process-0.png
 * @image html rtree-process-1.png
 * @image html rtree-process-2.png
 */
template <typename Rtree>
inline std::vector<BoundingBox<
  boost::geometry::dimension<typename Rtree::indexable_type>::value>>
extract_rtree_level(const Rtree &tree, const unsigned int level);



/**
 * Given an Rtree object @p tree and a target level @p level, this function returns
 * the bounding boxes associated to the children on level l+1 and stores them in
 * a vector. The resulting type is hence a vector of vectors of BoundingBox
 * objects.
 * If @p v is such a vector, then @p v has the following property:
 * @p v[i] is a vector with all of the bounding boxes associated to the children
 * of the i-th node of the tree.
 */
template <typename Rtree>
std::vector<std::vector<BoundingBox<
  boost::geometry::dimension<typename Rtree::indexable_type>::value>>>
extract_children_of_level(const Rtree &tree, const unsigned int level);



// Inline and template functions
#ifndef DOXYGEN

template <typename IndexType,
          typename LeafTypeIterator,
          typename IndexableGetter>
RTree<typename LeafTypeIterator::value_type, IndexType, IndexableGetter>
pack_rtree(const LeafTypeIterator &begin, const LeafTypeIterator &end)
{
  return RTree<typename LeafTypeIterator::value_type,
               IndexType,
               IndexableGetter>(begin, end);
}



template <typename IndexType, typename ContainerType, typename IndexableGetter>
RTree<typename ContainerType::value_type, IndexType, IndexableGetter>
pack_rtree(const ContainerType &container)
{
  return pack_rtree<IndexType, decltype(container.begin()), IndexableGetter>(
    container.begin(), container.end());
}



template <typename IndexType, typename ContainerType>
RTree<typename ContainerType::size_type,
      IndexType,
      IndexableGetterFromIndices<ContainerType>>
pack_rtree_of_indices(const ContainerType &container)
{
  // We need an array that holds the indices we want to pack. The rtree
  // implementation in BOOST, for reasons not entirely clear, insists
  // on using a reference to the elements of the range. This is fine if
  // the indices are stored in a container, so that's what we do.
  // (It would be nice if we could just pass a std::ranges::iota_view
  // instead, but that has no nested 'reference' type, and this then
  // trips up BOOST rtree.)
  std::vector<typename ContainerType::size_type> indices(container.size());
  for (typename ContainerType::size_type i = 0; i < container.size(); ++i)
    indices[i] = i;

  return RTree<typename ContainerType::size_type,
               IndexType,
               IndexableGetterFromIndices<ContainerType>>(
    indices.begin(),
    indices.end(),
    IndexType(),
    IndexableGetterFromIndices<ContainerType>(container));
}



template <typename Value,
          typename Options,
          typename Translator,
          typename Box,
          typename Allocators>
ExtractLevelVisitor<Value, Options, Translator, Box, Allocators>::
  ExtractLevelVisitor(
    const Translator  &translator,
    const unsigned int target_level,
    std::vector<BoundingBox<boost::geometry::dimension<Box>::value>> &boxes)
  : translator(translator)
  , level(0)
  , target_level(target_level)
  , boxes(boxes)
{}



template <typename Value,
          typename Options,
          typename Translator,
          typename Box,
          typename Allocators>
void
ExtractLevelVisitor<Value, Options, Translator, Box, Allocators>::operator()(
  const ExtractLevelVisitor::InternalNode &node)
{
  using ElementsType =
    typename boost::geometry::index::detail::rtree::elements_type<
      InternalNode>::type;

  const auto &elements = boost::geometry::index::detail::rtree::elements(node);

  if (level == target_level)
    {
      const auto offset = boxes.size();
      boxes.resize(offset + elements.size());

      unsigned int i = offset;
      for (typename ElementsType::const_iterator it = elements.begin();
           it != elements.end();
           ++it)
        {
          boost::geometry::convert(it->first, boxes[i]);
          ++i;
        }
      return;
    }

  const std::size_t level_backup = level;
  ++level;

  for (typename ElementsType::const_iterator it = elements.begin();
       it != elements.end();
       ++it)
    {
      boost::geometry::index::detail::rtree::apply_visitor(*this, *it->second);
    }

  level = level_backup;
}

template <typename Value,
          typename Options,
          typename Translator,
          typename Box,
          typename Allocators>
void
ExtractLevelVisitor<Value, Options, Translator, Box, Allocators>::operator()(
  const ExtractLevelVisitor::Leaf &)
{}



template <typename Rtree>
inline std::vector<BoundingBox<
  boost::geometry::dimension<typename Rtree::indexable_type>::value>>
extract_rtree_level(const Rtree &tree, const unsigned int level)
{
  constexpr unsigned int dim =
    boost::geometry::dimension<typename Rtree::indexable_type>::value;

  using RtreeView =
    boost::geometry::index::detail::rtree::utilities::view<Rtree>;
  RtreeView rtv(tree);

  std::vector<BoundingBox<dim>> boxes;

  if (rtv.depth() == 0)
    {
      // The below algorithm does not work for `rtv.depth()==0`, which might
      // happen if the number entries in the tree is too small.
      // In this case, simply return a single bounding box.
      boxes.resize(1);
      boost::geometry::convert(tree.bounds(), boxes[0]);
    }
  else
    {
      const unsigned int target_level =
        std::min<unsigned int>(level, rtv.depth() - 1);

      ExtractLevelVisitor<typename RtreeView::value_type,
                          typename RtreeView::options_type,
                          typename RtreeView::translator_type,
                          typename RtreeView::box_type,
                          typename RtreeView::allocators_type>
        extract_level_visitor(rtv.translator(), target_level, boxes);
      rtv.apply_visitor(extract_level_visitor);
    }

  return boxes;
}



template <class Rtree>
unsigned int
n_levels(const Rtree &tree)
{
  boost::geometry::index::detail::rtree::utilities::view<Rtree> rtv(tree);
  return rtv.depth();
}


template <typename Value,
          typename Options,
          typename Translator,
          typename Box,
          typename Allocators>
struct NodeVisitor : public boost::geometry::index::detail::rtree::visitor<
                       Value,
                       typename Options::parameters_type,
                       Box,
                       Allocators,
                       typename Options::node_tag,
                       true>::type
{
  inline NodeVisitor(
    const Translator &translator,
    unsigned int      target_level,
    std::vector<
      std::vector<BoundingBox<boost::geometry::dimension<Box>::value>>> &boxes);


  /**
   * An alias that identifies an InternalNode of the tree.
   */
  using InternalNode =
    typename boost::geometry::index::detail::rtree::internal_node<
      Value,
      typename Options::parameters_type,
      Box,
      Allocators,
      typename Options::node_tag>::type;

  /**
   * An alias that identifies a Leaf of the tree.
   */
  using Leaf = typename boost::geometry::index::detail::rtree::leaf<
    Value,
    typename Options::parameters_type,
    Box,
    Allocators,
    typename Options::node_tag>::type;

  /**
   * Implements the visitor interface for InternalNode objects. If the node
   * belongs to the level next to @p target_level, then fill the bounding box vector for that node.
   */
  inline void
  operator()(const InternalNode &node);

  /**
   * Implements the visitor interface for Leaf objects.
   */
  inline void
  operator()(const Leaf &);

  /**
   * Translator interface, required by the boost implementation of the rtree.
   */
  const Translator &translator;

  /**
   * Store the level we are currently visiting.
   */
  std::size_t level;

  /**
   * Index used to keep track of the number of different visited nodes during
   * recursion/
   */
  std::size_t node_counter;

  /**
   * The level where children are living.
   * Before: "we want to extract from the RTree object."
   */
  const std::size_t target_level;

  /**
   * A reference to the input vector of vector of BoundingBox objects. This
   * vector v has the following property: v[i] = vector with all
   * of the BoundingBox bounded by the i-th node of the Rtree.
   */
  std::vector<std::vector<BoundingBox<boost::geometry::dimension<Box>::value>>>
    &boxes_in_boxes;
};



template <typename Value,
          typename Options,
          typename Translator,
          typename Box,
          typename Allocators>
NodeVisitor<Value, Options, Translator, Box, Allocators>::NodeVisitor(
  const Translator  &translator,
  const unsigned int target_level,
  std::vector<std::vector<BoundingBox<boost::geometry::dimension<Box>::value>>>
    &bb_in_boxes)
  : translator(translator)
  , level(0)
  , node_counter(0)
  , target_level(target_level)
  , boxes_in_boxes(bb_in_boxes)
{}



template <typename Value,
          typename Options,
          typename Translator,
          typename Box,
          typename Allocators>
void
NodeVisitor<Value, Options, Translator, Box, Allocators>::operator()(
  const NodeVisitor::InternalNode &node)
{
  using elements_type =
    typename boost::geometry::index::detail::rtree::elements_type<
      InternalNode>::type; //  pairs of bounding box and pointer to child node
  const elements_type &elements =
    boost::geometry::index::detail::rtree::elements(node);

  if (level == target_level)
    {
      const unsigned int n_children = elements.size();
      const auto         offset     = boxes_in_boxes.size();
      boxes_in_boxes.resize(offset + n_children);
    }

  if (level == target_level + 1)
    {
      // I have now access to children of level target_level
      boxes_in_boxes[node_counter].resize(elements.size()); // number of bboxes
      unsigned int i = 0;
      for (typename elements_type::const_iterator it = elements.begin();
           it != elements.end();
           ++it)
        {
          boost::geometry::convert(it->first, boxes_in_boxes[node_counter][i]);
          ++i;
        }
      // Children have been stored, go to the next parent.
      ++node_counter;
      return;
    }

  std::size_t level_backup = level;
  ++level;

  for (typename elements_type::const_iterator it = elements.begin();
       it != elements.end();
       ++it)
    {
      boost::geometry::index::detail::rtree::apply_visitor(*this, *it->second);
    }

  level = level_backup;
}



template <typename Value,
          typename Options,
          typename Translator,
          typename Box,
          typename Allocators>
void
NodeVisitor<Value, Options, Translator, Box, Allocators>::operator()(
  const NodeVisitor::Leaf &)
{
  // No children for leaf nodes.
  boxes_in_boxes.clear();
}

template <typename Rtree>
inline std::vector<std::vector<BoundingBox<
  boost::geometry::dimension<typename Rtree::indexable_type>::value>>>
extract_children_of_level(const Rtree &tree, const unsigned int level)
{
  constexpr unsigned int dim =
    boost::geometry::dimension<typename Rtree::indexable_type>::value;

  using RtreeView =
    boost::geometry::index::detail::rtree::utilities::view<Rtree>;
  RtreeView rtv(tree);

  std::vector<std::vector<BoundingBox<dim>>> boxes_in_boxes;

  if (rtv.depth() == 0)
    {
      // The below algorithm does not work for `rtv.depth()==0`, which might
      // happen if the number entries in the tree is too small.
      // In this case, simply return a single bounding box.
      boxes_in_boxes.resize(1);
      boxes_in_boxes[0].resize(1);
      boost::geometry::convert(tree.bounds(), boxes_in_boxes[0][0]);
    }
  else
    {
      const unsigned int target_level =
        std::min<unsigned int>(level, rtv.depth() - 1);

      NodeVisitor<typename RtreeView::value_type,
                  typename RtreeView::options_type,
                  typename RtreeView::translator_type,
                  typename RtreeView::box_type,
                  typename RtreeView::allocators_type>
        node_visitor(rtv.translator(), target_level, boxes_in_boxes);
      rtv.apply_visitor(node_visitor);
    }

  return boxes_in_boxes;
}



#endif

DEAL_II_NAMESPACE_CLOSE

#endif
