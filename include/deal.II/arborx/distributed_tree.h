// ---------------------------------------------------------------------
//
// Copyright (C) 2022 by the deal.II authors
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

#ifndef dealii_arborx_distributed_tree_h
#define dealii_arborx_distributed_tree_h

#include <deal.II/base/config.h>

#if defined(DEAL_II_ARBORX_WITH_MPI) && defined(DEAL_II_WITH_MPI)
#  include <deal.II/arborx/access_traits.h>

#  include <ArborX_DistributedTree.hpp>
#  include <Kokkos_Core.hpp>

DEAL_II_NAMESPACE_OPEN

namespace ArborXWrappers
{
  /**
   * This class implements a wrapper around ArborX::DistributedTree, the
   * distributed version of ArborX:BVH.
   *
   * Because ArborX uses Kokkos, Kokkos needs to be initialized and finalized
   * before using this class.
   */
  class DistributedTree
  {
  public:
    /**
     * Constructor. Use a vector of BoundingBox @p bounding_boxes as primitives.
     * The BoundingBox objects in @p bounding_boxes are local to the MPI process.
     */
    template <int dim, typename Number>
    DistributedTree(
      const MPI_Comm &                             comm,
      const std::vector<BoundingBox<dim, Number>> &bounding_boxes);

    /**
     * Constructor. Use a vector of @p points as primitives. The Point objects
     * in @p points are local to the MPI process.
     */
    template <int dim, typename Number>
    DistributedTree(const MPI_Comm &                       comm,
                    const std::vector<Point<dim, Number>> &points);

    /**
     * Return the indices and the MPI ranks of those BoundingBox objects that
     * satisfy the @p queries.
     * Because @p queries can contain multiple queries, the function returns a pair
     * of indices and ranks and the associated offsets.
     *
     * Valid `QueryType` classes are
     * ArborXWrappers::BoundingBoxIntersectPredicate,
     * ArborXWrappers::BoundingBoxNearestPredicate,
     * ArborXWrappers::PointIntersectPredicate, and
     * ArborXWrappers::PointNearestPredicate,
     */
    template <typename QueryType>
    std::pair<std::vector<std::pair<int, int>>, std::vector<int>>
    query(const QueryType &queries);

  private:
    /**
     * Underlying ArborX object.
     */
    ArborX::DistributedTree<Kokkos::HostSpace> distributed_tree;
  };



  template <int dim, typename Number>
  DistributedTree::DistributedTree(
    const MPI_Comm &                             comm,
    const std::vector<BoundingBox<dim, Number>> &bounding_boxes)
    : distributed_tree(comm,
                       Kokkos::DefaultHostExecutionSpace{},
                       bounding_boxes)
  {}



  template <int dim, typename Number>
  DistributedTree::DistributedTree(
    const MPI_Comm &                       comm,
    const std::vector<Point<dim, Number>> &points)
    : distributed_tree(comm, Kokkos::DefaultHostExecutionSpace{}, points)
  {}



  template <typename QueryType>
  std::pair<std::vector<std::pair<int, int>>, std::vector<int>>
  DistributedTree::query(const QueryType &queries)
  {
    Kokkos::View<int *, Kokkos::HostSpace> offsets("offsets", 0);
    Kokkos::View<Kokkos::pair<int, int> *, Kokkos::HostSpace> indices_ranks(
      "indices_ranks", 0);
    distributed_tree.query(Kokkos::DefaultHostExecutionSpace{},
                           queries,
                           indices_ranks,
                           offsets);

    std::vector<std::pair<int, int>> indices_ranks_vector;
    for (unsigned int i = 0; i < indices_ranks.extent(0); ++i)
      {
        indices_ranks_vector.emplace_back(indices_ranks(i).first,
                                          indices_ranks(i).second);
      }

    std::vector<int> offsets_vector;
    offsets_vector.insert(offsets_vector.begin(),
                          offsets.data(),
                          offsets.data() + offsets.extent(0));

    return {indices_ranks_vector, offsets_vector};
  }
} // namespace ArborXWrappers

DEAL_II_NAMESPACE_CLOSE

#endif
#endif
