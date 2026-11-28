// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
// Copyright (C) 2001 - 2026 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Detailed license information governing the source code and contributions
// can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
//
// -----------------------------------------------------------------------------


#ifndef agglomerator_h
#define agglomerator_h


#include <deal.II/base/config.h>

#include <deal.II/base/bounding_box.h>

#include <boost/geometry/algorithms/distance.hpp>
#include <boost/geometry/index/rtree.hpp>
#include <boost/geometry/strategies/strategies.hpp>

namespace dealii
{
  template <int dim, int spacedim>
  class AgglomerationHandler;

  namespace internal
  {
    template <typename Value,
              typename Options,
              typename Translator,
              typename Box,
              typename Allocators>
    struct Rtree_visitor
      : public boost::geometry::index::detail::rtree::visitor<
          Value,
          typename Options::parameters_type,
          Box,
          Allocators,
          typename Options::node_tag,
          true>::type
    {
      inline Rtree_visitor(
        const Translator  &translator,
        const unsigned int target_level,
        std::vector<std::vector<typename Triangulation<
          boost::geometry::dimension<Box>::value>::active_cell_iterator>>
                                                        &agglomerates_,
        std::vector<types::global_cell_index>           &n_nodes_per_level,
        std::map<std::pair<types::global_cell_index, types::global_cell_index>,
                 std::vector<types::global_cell_index>> &parent_to_children);

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
       * belongs to the level next to @p target_level, then fill the bounding
       * box vector for that node.
       */
      inline void
      operator()(const InternalNode &node);

      /**
       * Implements the visitor interface for Leaf objects.
       */
      inline void
      operator()(const Leaf &);

      /**
       * Translator interface, required by the boost implementation of the
       * rtree.
       */
      const Translator &translator;

      /**
       * Store the level we are currently visiting.
       */
      size_t level;

      /**
       * Index used to keep track of the number of different visited nodes
       * during recursion/
       */
      size_t node_counter;

      /**
       * The level where children are living.
       * Before: "we want to extract from the RTree object."
       */
      const size_t target_level;

      /**
       * A reference to the input vector of vector of BoundingBox objects. This
       * vector v has the following property: v[i] = vector with all the mesh
       * iterators composing the i-th agglomerate.
       */
      std::vector<std::vector<typename Triangulation<
        boost::geometry::dimension<Box>::value>::active_cell_iterator>>
        &agglomerates;

      /**
       * Store the total number of nodes on each level.
       */
      std::vector<types::global_cell_index> &n_nodes_per_level;

      /**
       * Map that associates to a given node on level l its children, identified
       * by their integer index.
       */
      std::map<std::pair<types::global_cell_index, types::global_cell_index>,
               std::vector<types::global_cell_index>>
        &parent_node_to_children_nodes;
    };



    template <typename Value,
              typename Options,
              typename Translator,
              typename Box,
              typename Allocators>
    Rtree_visitor<Value, Options, Translator, Box, Allocators>::Rtree_visitor(
      const Translator  &translator,
      const unsigned int target_level,
      std::vector<std::vector<typename Triangulation<
        boost::geometry::dimension<Box>::value>::active_cell_iterator>>
                                                      &agglomerates_,
      std::vector<types::global_cell_index>           &n_nodes_per_level_,
      std::map<std::pair<types::global_cell_index, types::global_cell_index>,
               std::vector<types::global_cell_index>> &parent_to_children)
      : translator(translator)
      , level(0)
      , node_counter(0)
      , target_level(target_level)
      , agglomerates(agglomerates_)
      , n_nodes_per_level(n_nodes_per_level_)
      , parent_node_to_children_nodes(parent_to_children)
    {}



    template <typename Value,
              typename Options,
              typename Translator,
              typename Box,
              typename Allocators>
    void
    Rtree_visitor<Value, Options, Translator, Box, Allocators>::operator()(
      const Rtree_visitor::InternalNode &node)
    {
      using elements_type =
        typename boost::geometry::index::detail::rtree::elements_type<
          InternalNode>::type; //  pairs of bounding box and pointer to child
                               //  node
      const elements_type &elements =
        boost::geometry::index::detail::rtree::elements(node);

      if (level < target_level)
        {
          size_t level_backup = level;
          ++level;

          for (typename elements_type::const_iterator it = elements.begin();
               it != elements.end();
               ++it)
            {
              boost::geometry::index::detail::rtree::apply_visitor(*this,
                                                                   *it->second);
            }

          level = level_backup;
        }
      else if (level == target_level)
        {
          const auto offset = agglomerates.size();
          agglomerates.resize(offset + 1);
          size_t level_backup = level;

          ++level;
          for (const auto &entry : elements)
            {
              boost::geometry::index::detail::rtree::apply_visitor(
                *this, *entry.second);
            }
          // Done with node number 'node_counter' on level target_level.

          ++node_counter; // visited all children of an internal node
          n_nodes_per_level[target_level]++;

          level = level_backup;
        }
      else if (level > target_level)
        {
          // I am on a child (internal) node on a deeper level.

          // Keep visiting until you go to the leafs.
          size_t level_backup = level;

          ++level;

          // looping through entries of node
          for (const auto &entry : elements)
            {
              boost::geometry::index::detail::rtree::apply_visitor(
                *this, *entry.second);
            }
          // done with node on level l > target_level (not just
          // "target_level+1).
          n_nodes_per_level[level_backup]++;
          const types::global_cell_index node_idx =
            n_nodes_per_level[level_backup] - 1; // so to start from 0

          parent_node_to_children_nodes[{n_nodes_per_level[level_backup - 1],
                                         level_backup - 1}]
            .push_back(node_idx);

          level = level_backup;
        }
    }



    template <typename Value,
              typename Options,
              typename Translator,
              typename Box,
              typename Allocators>
    void
    Rtree_visitor<Value, Options, Translator, Box, Allocators>::operator()(
      const Rtree_visitor::Leaf &leaf)
    {
      using elements_type =
        typename boost::geometry::index::detail::rtree::elements_type<
          Leaf>::type; //  pairs of bounding box and pointer to child node
      const elements_type &elements =
        boost::geometry::index::detail::rtree::elements(leaf);

      if (level == target_level)
        {
          // If I want to extract from leaf node, i.e. the target_level is the
          // last one where leafs are grouped together.
          const auto offset = agglomerates.size();
          agglomerates.resize(offset + 1);

          for (const auto &it : elements)
            agglomerates[node_counter].push_back(it.second);

          ++node_counter;
          n_nodes_per_level[target_level]++;
        }
      else
        {
          for (const auto &it : elements)
            agglomerates[node_counter].push_back(it.second);


          if (level == target_level + 1)
            {
              const unsigned int node_idx = n_nodes_per_level[level];

              parent_node_to_children_nodes[{n_nodes_per_level[level - 1],
                                             level - 1}]
                .push_back(node_idx);
              n_nodes_per_level[level]++;
            }
        }
    }
  } // namespace internal



  /**
   * Helper class which handles agglomeration based on the R-tree data
   * structure. Notice that the R-tree type is assumed to be an R-star-tree.
   */
  template <int dim, typename RtreeType>
  class CellsAgglomerator
  {
  public:
    template <int, int>
    friend class AgglomerationHandler;

    /**
     * Constructor. It takes a given rtree and an integer representing the
     * index of the level to be extracted.
     */
    CellsAgglomerator(const RtreeType   &rtree,
                      const unsigned int extraction_level);

    /**
     * Extract agglomerates based on the current tree and the extraction level.
     * This function returns a reference to
     */
    const std::vector<
      std::vector<typename Triangulation<dim>::active_cell_iterator>> &
    extract_agglomerates();

    /**
     * Get total number of levels.
     */
    inline unsigned int
    get_n_levels() const;

    /**
     * Return the number of nodes present in level @p level.
     */
    inline types::global_cell_index
    get_n_nodes_per_level(const unsigned int level) const;

    /**
     * This function returns a map which associates to each node on level
     * @p extraction_level a list of children.
     */
    inline const std::map<
      std::pair<types::global_cell_index, types::global_cell_index>,
      std::vector<types::global_cell_index>> &
    get_hierarchy() const;

  private:
    /**
     * Raw pointer to the actual R-tree.
     */
    RtreeType *rtree;

    /**
     * Extraction level.
     */
    const unsigned int extraction_level;

    /**
     * Store agglomerates obtained after recursive extraction on nodes of
     * level @p extraction_level.
     */
    std::vector<std::vector<typename Triangulation<dim>::active_cell_iterator>>
      agglomerates_on_level;

    /**
     * Vector storing the number of nodes (and, ultimately, agglomerates) for
     * each level.
     */
    std::vector<types::global_cell_index> n_nodes_per_level;

    /**
     * Map which maps a node parent @n on level @p l to a vector of integers
     * which stores the index of children.
     */
    std::map<std::pair<types::global_cell_index, types::global_cell_index>,
             std::vector<types::global_cell_index>>
      parent_node_to_children_nodes;
  };



  template <int dim, typename RtreeType>
  CellsAgglomerator<dim, RtreeType>::CellsAgglomerator(
    const RtreeType   &tree,
    const unsigned int extraction_level_)
    : extraction_level(extraction_level_)
  {
    rtree = const_cast<RtreeType *>(&tree);
    Assert(n_levels(*rtree), ExcMessage("At least two levels are needed."));
  }



  template <int dim, typename RtreeType>
  const std::vector<
    std::vector<typename Triangulation<dim>::active_cell_iterator>> &
  CellsAgglomerator<dim, RtreeType>::extract_agglomerates()
  {
    AssertThrow(extraction_level <= n_levels(*rtree),
                ExcInternalError("You are trying to extract level " +
                                 std::to_string(extraction_level) +
                                 " of the tree, but it only has a total of " +
                                 std::to_string(n_levels(*rtree)) +
                                 " levels."));
    using RtreeView =
      boost::geometry::index::detail::rtree::utilities::view<RtreeType>;
    RtreeView rtv(*rtree);

    n_nodes_per_level.resize(rtv.depth() +
                             1); // store how many nodes we have for each level.

    if (rtv.depth() == 0)
      {
        // The below algorithm does not work for `rtv.depth()==0`, which might
        // happen if the number entries in the tree is too small.
        agglomerates_on_level.resize(1);
        agglomerates_on_level[0].resize(1);
      }
    else
      {
        const unsigned int target_level =
          std::min<unsigned int>(extraction_level, rtv.depth());

        internal::Rtree_visitor<typename RtreeView::value_type,
                                typename RtreeView::options_type,
                                typename RtreeView::translator_type,
                                typename RtreeView::box_type,
                                typename RtreeView::allocators_type>
          extractor_visitor(rtv.translator(),
                            target_level,
                            agglomerates_on_level,
                            n_nodes_per_level,
                            parent_node_to_children_nodes);


        rtv.apply_visitor(extractor_visitor);
      }
    return agglomerates_on_level;
  }



  // ------------------------------ inline functions -------------------------


  template <int dim, typename RtreeType>
  inline unsigned int
  CellsAgglomerator<dim, RtreeType>::get_n_levels() const
  {
    return n_levels(*rtree);
  }



  template <int dim, typename RtreeType>
  inline types::global_cell_index
  CellsAgglomerator<dim, RtreeType>::get_n_nodes_per_level(
    const unsigned int level) const
  {
    return n_nodes_per_level[level];
  }



  template <int dim, typename RtreeType>
  inline const std::map<
    std::pair<types::global_cell_index, types::global_cell_index>,
    std::vector<types::global_cell_index>> &
  CellsAgglomerator<dim, RtreeType>::get_hierarchy() const
  {
    Assert(parent_node_to_children_nodes.size(),
           ExcMessage(
             "The hierarchy has not been computed. Did you forget to call"
             " extract_agglomerates() first?"));
    return parent_node_to_children_nodes;
  }
} // namespace dealii
#endif