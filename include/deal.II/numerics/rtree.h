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
#include <boost/geometry/algorithms/buffer.hpp>
#include <boost/geometry/algorithms/distance.hpp>
#include <boost/geometry/index/rtree.hpp>
#include <boost/geometry/strategies/strategies.hpp>
DEAL_II_ENABLE_EXTRA_DIAGNOSTICS

#include <memory>
#include <type_traits>
#include <utility>


DEAL_II_NAMESPACE_OPEN

namespace internal
{
  namespace RTreeImplementation
  {
    /**
     * Sentinel type used by RTreeFunctionalVisitor to indicate that no callback
     * was provided.
     */
    struct EmptyVisitor
    {};
  } // namespace RTreeImplementation
} // namespace internal

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
 * An alias to a view of an RTree.
 *
 * This is a convenience alias to refer to the view type of an RTree. This is
 * the only way to access the internal structure of the RTree in a read-only
 * fashion.
 *
 * @tparam RTreeType
 */
template <typename RTreeType>
using RTreeView =
  boost::geometry::index::detail::rtree::utilities::view<RTreeType>;

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
 * A depth-first visitor for RTree that invokes user-provided lambdas.
 *
 * For every RTree object, we distinguish between the root node, the internal
 * nodes (i.e., nodes that have other nodes as children) and leaf nodes (i.e.,
 * nodes that store the actual leaf values). During the traversal of the tree,
 * we may want to hook into the traversal at different levels of the hierarchy.
 * For example, we may want to print the bounding box of each internal node,
 * then print the bounding box of each leaf node, and then finally print each
 * element stored in the leaf nodes.
 *
 * This is a convenience visitor for the RTree data structure that forwards each
 * traversal event to four (optional) independent callables:
 * - an internal node visitor (called once per internal node, before
 *   descending),
 * - a leaf node visitor (called once per leaf node, before visiting its leaf
 *   elements),
 * - an element visitor (called once per leaf value stored in a leaf node),
 * - an indexable visitor (called once per leaf value, after translating Value
 *   -> Indexable with the tree's indexable getter).
 *
 * The traversal is depth-first: starting from the root, internal nodes are
 * visited, their children are recursed immediately, and leaf elements are
 * delivered in the order stored inside each leaf node. This allows users to
 * plug lambdas directly without defining a specific visitor class.
 *
 * Example (taken from one of the tests): build a tree from 8 random points,
 * then hook all four callbacks:
 *
 * @code
 * std::vector<Point<2>> pts(8); for (auto &p : pts) p = random_point<2>();
 *
 * auto tree = pack_rtree_of_indices<bgi::linear<2>>(pts);
 *
 * // internal node visitor auto internal = [](std::size_t level, auto &parent,
 * const BoundingBox<2> &box, auto &child) { std::cout << "Internal node (level
 * " << level << ") box " <<
 * Patterns::Tools::to_string(box.get_boundary_points()) << std::endl;
 * };
 *
 * // leaf visitor auto leaf = [](std::size_t level, auto &, const
 * BoundingBox<2> &box, auto &) { std::cout << "    Leaf      (level " << level
 * << ") box " << Patterns::Tools::to_string(box.get_boundary_points()) <<
 * std::endl;
 * };
 *
 * // element visitor (receives the stored value) auto element = [](auto &,
 * const auto &value) { std::cout << "        element " << value;
 * };
 *
 * // indexable visitor (receives the indexable geometry) auto indexable =
 * [](auto &, const auto &indexable) { std::cout << ",  indexable " << indexable
 * << std::endl;
 * };
 *
 * std::cout << "Visiting R-tree with " << pts.size() << " points, and " <<
 *           n_levels(tree) << " levels." << std::endl;
 *
 * visit_rtree(tree, internal, leaf, element, indexable);
 * @endcode
 *
 * A possible output (indices and boxes depend on the random points) is:
 * @code
 * Internal node (level 0) box 0.277775, 0.197551 : 0.911647, 0.55397 Leaf
 *     (level 1) box 0.277775, 0.513401 : 0.364784, 0.55397 element 4, indexable
 *     0.277775 0.55397 element 6,  indexable 0.364784 0.513401 Leaf (level 1)
 *     box 0.840188, 0.197551 : 0.911647, 0.394383 element 0, indexable 0.840188
 *     0.394383 element 2,  indexable 0.911647 0.197551 Internal node (level 0)
 *     box 0.335223, 0.628871 : 0.95223, 0.916195 Leaf (level 1) box 0.335223,
 *     0.628871 : 0.477397, 0.76823 element 3, indexable 0.335223 0.76823
 *     element 5,  indexable 0.477397 0.628871 Leaf (level 1) box 0.783099,
 *     0.79844 : 0.95223, 0.916195 element 1,  indexable 0.783099 0.79844
 *     element 7,  indexable 0.95223 0.916195
 * @endcode
 *
 * This structure exists to make it easy to instrument or analyze an rtree
 * (bounding boxes, counts, aggregations, debug prints, etc.) without writing
 * the boilerplate of a custom visitor type every time.
 *
 * Users can call directly the function visit_rtree() (see below) instead of
 * instantiating this class.
 */
template <
  typename RTreeType,
  typename InternalVisitor  = internal::RTreeImplementation::EmptyVisitor,
  typename LeafVisitor      = internal::RTreeImplementation::EmptyVisitor,
  typename ElementVisitor   = internal::RTreeImplementation::EmptyVisitor,
  typename IndexableVisitor = internal::RTreeImplementation::EmptyVisitor>
struct RTreeFunctionalVisitor
  : public boost::geometry::index::detail::rtree::visitor<
      typename RTreeView<RTreeType>::value_type,
      typename RTreeView<RTreeType>::options_type::parameters_type,
      typename RTreeView<RTreeType>::box_type,
      typename RTreeView<RTreeType>::allocators_type,
      typename RTreeView<RTreeType>::options_type::node_tag,
      true>
{
  /**
   * Spatial dimension of the indexable geometry stored in the tree.
   */
  static constexpr unsigned int dim =
    boost::geometry::dimension<typename RTreeType::indexable_type>::value;

  /**
   * Internal node type exposed by Boost for this R-tree configuration.
   */
  using internal_node_type =
    typename boost::geometry::index::detail::rtree::internal_node<
      typename RTreeView<RTreeType>::value_type,
      typename RTreeView<RTreeType>::options_type::parameters_type,
      typename RTreeView<RTreeType>::box_type,
      typename RTreeView<RTreeType>::allocators_type,
      typename RTreeView<RTreeType>::options_type::node_tag>::type;

  /**
   * Leaf type exposed by Boost for this R-tree configuration
   */
  using leaf_type = typename boost::geometry::index::detail::rtree::leaf<
    typename RTreeView<RTreeType>::value_type,
    typename RTreeView<RTreeType>::options_type::parameters_type,
    typename RTreeView<RTreeType>::box_type,
    typename RTreeView<RTreeType>::allocators_type,
    typename RTreeView<RTreeType>::options_type::node_tag>::type;

  /**
   * Value stored in leaves (matches RTreeType::value_type).
   */
  using value_type = typename RTreeType::value_type;

  /**
   * Indexable geometry obtained from a value via the tree's getter.
   */
  using indexable_type = typename RTreeType::indexable_type;

  /**
   * Translator (indexable getter) associated with the tree view.
   */
  using translator_type = typename RTreeView<RTreeType>::translator_type;

  /**
   * Sentinel type used to detect empty visitors.
   */
  using empty_visitor = internal::RTreeImplementation::EmptyVisitor;

  /**
   * Compile-time flags indicating if internal visitor was provided.
   */
  static constexpr bool has_internal_visitor =
    !std::is_same_v<InternalVisitor, empty_visitor>;

  /**
   * Compile-time flags indicating if leaf visitor was provided.
   */
  static constexpr bool has_leaf_visitor =
    !std::is_same_v<LeafVisitor, empty_visitor>;

  /**
   * Compile-time flags indicating if element visitor was provided.
   */
  static constexpr bool has_element_visitor =
    !std::is_same_v<ElementVisitor, empty_visitor>;

  /**
   * Compile-time flags indicating if indexable visitor was provided.
   */
  static constexpr bool has_indexable_visitor =
    !std::is_same_v<IndexableVisitor, empty_visitor>;

  /**
   * Optional callback invoked once per internal node, before descending into
   * its children.
   *
   * This function is only called when the child is an internal_node_type. If
   * the child is a leaf_type, the leaf_visitor is called instead.
   *
   * Arguments:
   * - level: depth excluding the root (root's children are at level 0; the
   *   maximum value is @p n_levels(tree)-1),
   * - parent_node: the internal node being visited (the current parent),
   * - node_bounding_box: bounding box of the child about to be traversed,
   * - node: the child internal node
   */
  InternalVisitor internal_node_visitor;

  /**
   * Optional callback invoked once per leaf node, before iterating its
   * elements.
   *
   * Arguments mirror the internal callback, but the last argument is the child
   * leaf node.
   */
  LeafVisitor leaf_visitor;

  /**
   * Optional callback invoked for each value stored in a leaf.
   *
   * Arguments:
   * - leaf: the leaf holding the value,
   * - value: the stored value_type.
   */
  ElementVisitor element_visitor;

  /**
   * Optional callback invoked for each value after translating it to the
   * corresponding indexable type via the tree's indexable getter. If the tree
   * is built without an indexable getter, this function is in fact identical to
   * the element_visitor, and therefore both can be used interchangeably.
   *
   * Arguments:
   * - leaf: the leaf holding the value,
   * - indexable: the Indexable associated with the value.
   */
  IndexableVisitor indexable_visitor;

  /**
   * Implements the internal-node branch of the
   * `boost::geometry::index::detail::rtree::visitor` interface so this
   * functor can be passed to apply_visitor() when traversing an rtree.
   */
  void
  operator()(internal_node_type &node);

  /**
   * Implements the leaf branch of the
   * `boost::geometry::index::detail::rtree::visitor` interface so this
   * functor can be passed to apply_visitor() when traversing an rtree.
   */
  void
  operator()(leaf_type &leaf);

  /**
   * Constructor wiring optional callbacks and starting a traversal of the
   * provided rtree through apply_visitor().
   *
   * You may call this constructor directly, but the preferred way is to use
   * the visit_rtree() function below.
   */
  RTreeFunctionalVisitor(const RTreeType &tree,
                         InternalVisitor  internal_node_visitor = {},
                         LeafVisitor      leaf_visitor          = {},
                         ElementVisitor   element_visitor       = {},
                         IndexableVisitor indexable_visitor     = {});

  std::size_t     current_level;
  translator_type translator;
};

/**
 * Depth-first traversal helper for any R-tree instance.
 *
 * This function instantiates RTreeFunctionalVisitor<RTreeType> and walks the
 * tree depth-first, invoking the provided callbacks at three granularities:
 * internal nodes, leaf nodes, and leaf elements (values and/or indexables).
 *
 * @tparam RTreeType The concrete boost::geometry::index::rtree type.
 * @param tree The R-tree instance to traverse (passed by non-const reference so
 * the visitor signature matches Boost's visitor API; traversal itself is
 * read-only).
 * @param internal_node_visitor Optional callback invoked once per internal node
 * before its children are visited. Signature:
 * @code
 * void(std::size_t level,
 *      const internal_node_type &parent,
 *      const BoundingBox<dim> &child_box,
 *      const internal_node_type &child);
 * @endcode
 * The @p level is the depth of the child without counting the root (root's
 * children are level 0, and the maximum value is @p n_levels(tree)-1). The
 * @p child_box is the bounding box stored alongside the child pointer in the
 * parent.
 * @param leaf_visitor Optional callback invoked once per leaf node before its
 * elements are iterated. Signature:
 * @code
 * void(std::size_t level,
 *      const internal_node_type &parent,
 *      const BoundingBox<dim> &leaf_box,
 *      const leaf_type &leaf);
 * @endcode
 * @param element_visitor Optional callback invoked for each value stored in a
 * leaf. Signature:
 * @code
 * void(const leaf_type &leaf,
 *      const value_type &value);
 * @endcode
 * @param indexable_visitor Optional callback invoked for each value after it is
 * translated to its indexable counterpart via the tree's indexable getter.
 * Signature:
 * @code
 * void(const leaf_type &leaf,
 *      const indexable_type &indexable);
 * @endcode
 *
 * The traversal order is depth-first: an internal node triggers its callback
 * and immediately recurses into each child (internal or leaf); leaves trigger
 * their callback and then the element/indexable callbacks for each stored
 * value.
 *
 * For an end-to-end usage example and sample output, see the documentation of
 * RTreeFunctionalVisitor.
 */
template <
  typename RTreeType,
  typename InternalVisitor  = internal::RTreeImplementation::EmptyVisitor,
  typename LeafVisitor      = internal::RTreeImplementation::EmptyVisitor,
  typename ElementVisitor   = internal::RTreeImplementation::EmptyVisitor,
  typename IndexableVisitor = internal::RTreeImplementation::EmptyVisitor>
void
visit_rtree(const RTreeType &tree,
            InternalVisitor  internal_node_visitor = {},
            LeafVisitor      leaf_visitor          = {},
            ElementVisitor   element_visitor       = {},
            IndexableVisitor indexable_visitor     = {});


/**
 * Given a RTree object @p rtree, and a target level @p level, return a vector
 * of BoundingBox objects containing all the bounding boxes that make the given
 * @p level of the @p rtree. This function is used in step-70.
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



template <typename RTreeType,
          typename InternalVisitor,
          typename LeafVisitor,
          typename ElementVisitor,
          typename IndexableVisitor>
void
RTreeFunctionalVisitor<RTreeType,
                       InternalVisitor,
                       LeafVisitor,
                       ElementVisitor,
                       IndexableVisitor>::operator()(internal_node_type &node)
{
  const std::size_t level_backup = current_level;
  for (const auto &element :
       boost::geometry::index::detail::rtree::elements(node))
    {
      BoundingBox<dim> box;
      boost::geometry::convert(element.first, box);

      const std::size_t child_level = level_backup;

      if (auto *leaf_ptr = boost::get<leaf_type>(&*element.second))
        {
          if constexpr (has_leaf_visitor)
            leaf_visitor(child_level, node, box, *leaf_ptr);
        }
      else if (auto *internal_ptr =
                 boost::get<internal_node_type>(&*element.second))
        {
          if constexpr (has_internal_visitor)
            internal_node_visitor(child_level, node, box, *internal_ptr);
        }

      current_level = child_level + 1;
      boost::geometry::index::detail::rtree::apply_visitor(*this,
                                                           *element.second);
    }
  current_level = level_backup;
}



template <typename RTreeType,
          typename InternalVisitor,
          typename LeafVisitor,
          typename ElementVisitor,
          typename IndexableVisitor>
void
RTreeFunctionalVisitor<RTreeType,
                       InternalVisitor,
                       LeafVisitor,
                       ElementVisitor,
                       IndexableVisitor>::operator()(leaf_type &leaf)
{
  for (const auto &element :
       boost::geometry::index::detail::rtree::elements(leaf))
    {
      if constexpr (has_element_visitor)
        element_visitor(leaf, element);
      if constexpr (has_indexable_visitor)
        {
          auto idx = translator(element);
          indexable_visitor(leaf, idx);
        }
    }
}



template <typename RTreeType,
          typename InternalVisitor,
          typename LeafVisitor,
          typename ElementVisitor,
          typename IndexableVisitor>
RTreeFunctionalVisitor<RTreeType,
                       InternalVisitor,
                       LeafVisitor,
                       ElementVisitor,
                       IndexableVisitor>::
  RTreeFunctionalVisitor(const RTreeType &tree,
                         InternalVisitor  internal_node_visitor,
                         LeafVisitor      leaf_visitor,
                         ElementVisitor   element_visitor,
                         IndexableVisitor indexable_visitor)
  : internal_node_visitor(std::move(internal_node_visitor))
  , leaf_visitor(std::move(leaf_visitor))
  , element_visitor(std::move(element_visitor))
  , indexable_visitor(std::move(indexable_visitor))
  , current_level(0)
  , translator(RTreeView<RTreeType>(tree).translator())
{
  static_assert(!has_internal_visitor ||
                  std::is_invocable_r_v<void,
                                        InternalVisitor,
                                        const std::size_t &,
                                        const internal_node_type &,
                                        const BoundingBox<dim> &,
                                        const internal_node_type &>,
                "internal_node_visitor must be callable with (std::size_t, "
                "internal_node_type, BoundingBox<dim>, internal_node_type)");

  static_assert(
    !has_leaf_visitor || std::is_invocable_r_v<void,
                                               LeafVisitor,
                                               const std::size_t &,
                                               const internal_node_type &,
                                               const BoundingBox<dim> &,
                                               const leaf_type &>,
    "leaf_visitor must be callable with (std::size_t, internal_node_type, "
    "BoundingBox<dim>, leaf_type)");

  static_assert(
    !has_element_visitor || std::is_invocable_r_v<void,
                                                  ElementVisitor,
                                                  const leaf_type &,
                                                  const value_type &>,
    "element_visitor must be callable with (leaf_type, value_type)");

  static_assert(
    !has_indexable_visitor || std::is_invocable_r_v<void,
                                                    IndexableVisitor,
                                                    const leaf_type &,
                                                    const indexable_type &>,
    "indexable_visitor must be callable with (leaf_type, indexable_type)");

  RTreeView<RTreeType> rtv(tree);
  rtv.apply_visitor(*this);
}



template <typename RTreeType,
          typename InternalVisitor,
          typename LeafVisitor,
          typename ElementVisitor,
          typename IndexableVisitor>
void
visit_rtree(const RTreeType &tree,
            InternalVisitor  internal_node_visitor,
            LeafVisitor      leaf_visitor,
            ElementVisitor   element_visitor,
            IndexableVisitor indexable_visitor)
{
  RTreeFunctionalVisitor<RTreeType,
                         std::decay_t<InternalVisitor>,
                         std::decay_t<LeafVisitor>,
                         std::decay_t<ElementVisitor>,
                         std::decay_t<IndexableVisitor>>
    visitor(tree,
            std::move(internal_node_visitor),
            std::move(leaf_visitor),
            std::move(element_visitor),
            std::move(indexable_visitor));
}



template <typename Rtree>
inline std::vector<BoundingBox<
  boost::geometry::dimension<typename Rtree::indexable_type>::value>>
extract_rtree_level(const Rtree &tree, const unsigned int level)
{
  constexpr unsigned int dim =
    boost::geometry::dimension<typename Rtree::indexable_type>::value;

  std::vector<BoundingBox<dim>> boxes;

  unsigned int n_levels_in_tree = n_levels(tree);

  if (n_levels_in_tree == 0)
    {
      // The algorithm below does not work for `n_levels_in_tree==0`, which
      // might happen if the number entries in the tree is too small. In this
      // case, simply return a single bounding box.
      boxes.resize(1);
      boost::geometry::convert(tree.bounds(), boxes[0]);
    }
  else
    {
      const unsigned int target_level =
        std::min<unsigned int>(level, n_levels_in_tree - 1);

      using Visitor           = RTreeFunctionalVisitor<Rtree>;
      using InternalNode      = typename Visitor::internal_node_type;
      using LeafNode          = typename Visitor::leaf_type;
      const auto node_visitor = [&](const std::size_t &current_level,
                                    const InternalNode &,
                                    const BoundingBox<dim> &box,
                                    const InternalNode &) {
        if (current_level == target_level)
          boxes.push_back(box);
      };
      const auto leaf_visitor = [&](const std::size_t &current_level,
                                    const InternalNode &,
                                    const BoundingBox<dim> &box,
                                    const LeafNode &) {
        if (current_level == target_level)
          boxes.push_back(box);
      };
      visit_rtree(tree, node_visitor, leaf_visitor);
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


template <typename Rtree>
inline std::vector<std::vector<BoundingBox<
  boost::geometry::dimension<typename Rtree::indexable_type>::value>>>
extract_children_of_level(const Rtree &tree, const unsigned int level)
{
  constexpr unsigned int dim =
    boost::geometry::dimension<typename Rtree::indexable_type>::value;

  unsigned int n_levels_in_tree = n_levels(tree);

  std::vector<std::vector<BoundingBox<dim>>> boxes_in_boxes;

  if (n_levels_in_tree == 0)
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
        std::min<unsigned int>(level, n_levels_in_tree - 1);

      // Track all nodes on target_level and collect their children (level
      // target_level+1) grouped per parent. Depth-first visitation guarantees
      // children immediately follow their parent, so a simple counter suffices.
      std::size_t current_parent_index = 0;
      std::size_t node_counter         = 0;

      const auto internal = [&](const std::size_t current_level,
                                const auto &,
                                const BoundingBox<dim> &child_box,
                                const auto &) {
        // Register nodes that live on the target level so their children can
        // be grouped.
        if (current_level == target_level)
          {
            current_parent_index = node_counter++;
            boxes_in_boxes.emplace_back();
            return;
          }

        // Record children of nodes on the target level.
        if (current_level == target_level + 1)
          boxes_in_boxes[current_parent_index].push_back(child_box);
      };

      const auto leaf = [&](const std::size_t current_level,
                            const auto &,
                            const BoundingBox<dim> &child_box,
                            const auto &) {
        if (current_level == target_level + 1)
          boxes_in_boxes[current_parent_index].push_back(child_box);
      };

      visit_rtree(tree, internal, leaf);
    }

  return boxes_in_boxes;
}



#endif

DEAL_II_NAMESPACE_CLOSE

#endif
