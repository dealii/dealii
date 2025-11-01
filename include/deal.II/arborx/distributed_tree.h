// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2022 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#ifndef dealii_arborx_distributed_tree_h
#define dealii_arborx_distributed_tree_h

#include <deal.II/base/config.h>

#if defined(DEAL_II_ARBORX_WITH_MPI) && defined(DEAL_II_WITH_MPI)
#  include <deal.II/arborx/access_traits.h>

#  include <deal.II/base/mpi.h>

#  include <ArborX_DistributedTree.hpp>
#  include <Kokkos_Core.hpp>

DEAL_II_NAMESPACE_OPEN

namespace ArborXWrappers
{
  /**
   * This class implements a wrapper around ArborX::DistributedTree, the
   * distributed version of ArborX::BVH.
   */
#  if ARBORX_VERSION_MAJOR < 2
  class DistributedTree
  {
  public:
    /**
     * Constructor. Use a vector of BoundingBox @p bounding_boxes as primitives.
     * The BoundingBox objects in @p bounding_boxes are local to the MPI process.
     */
    template <int dim, typename Number>
    DistributedTree(
      const MPI_Comm                               comm,
      const std::vector<BoundingBox<dim, Number>> &bounding_boxes);

    /**
     * Constructor. Use a vector of @p points as primitives. The Point objects
     * in @p points are local to the MPI process.
     */
    template <int dim, typename Number>
    DistributedTree(const MPI_Comm                         comm,
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
     * ArborXWrappers::PointIntersectPredicate,
     * ArborXWrappers::PointNearestPredicate,
     * ArborXWrappers::SphereIntersectPredicate, and
     * ArborXWrappers::SphereNearestPredicate.
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
    const MPI_Comm                               comm,
    const std::vector<BoundingBox<dim, Number>> &bounding_boxes)
    : distributed_tree(comm,
                       Kokkos::DefaultHostExecutionSpace{},
                       bounding_boxes)
  {}



  template <int dim, typename Number>
  DistributedTree::DistributedTree(
    const MPI_Comm                         comm,
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
#  else
  template <typename Value>
  class DistributedTree
  {
  public:
    /**
     * Constructor. Use a vector of Point, BoundingBox, or sphere
     * (std::pair<Point,Number>) @p values as primitives. The objects in @p values
     * are local to the MPI process.
     */
    DistributedTree(const MPI_Comm comm, const std::vector<Value> &values);

    /**
     * Return the indices and the MPI ranks of those BoundingBox objects that
     * satisfy the @p queries.
     * Because @p queries can contain multiple queries, the function returns a pair
     * of indices and ranks and the associated offsets.
     *
     * Valid `QueryType` classes are
     * ArborXWrappers::BoundingBoxIntersectPredicate,
     * ArborXWrappers::BoundingBoxNearestPredicate,
     * ArborXWrappers::PointIntersectPredicate,
     * ArborXWrappers::PointNearestPredicate,
     * ArborXWrappers::SphereIntersectPredicate, and
     * ArborXWrappers::SphereNearestPredicate.
     */
    template <typename QueryType>
    std::pair<std::vector<std::pair<int, int>>, std::vector<int>>
    query(const QueryType &queries);

  private:
    /**
     * Underlying ArborX object.
     */
    ArborX::DistributedTree<Kokkos::HostSpace,
                            ArborX::PairValueIndex<Value, unsigned int>,
                            internal::IndexableGetter>
      distributed_tree;
    /**
     * Callback to extract the index and the rank of each primitive that
     * satisfies a query.
     */
    internal::ExtractIndexRank callback;
  };



  template <typename Value>
  DistributedTree<Value>::DistributedTree(const MPI_Comm            comm,
                                          const std::vector<Value> &values)
    : distributed_tree(comm,
                       Kokkos::DefaultHostExecutionSpace{},
                       ArborX::Experimental::attach_indices(values),
                       internal::IndexableGetter{})
    , callback{Utilities::MPI::this_mpi_process(comm)}
  {}

  template <typename Value>
  template <typename QueryType>
  std::pair<std::vector<std::pair<int, int>>, std::vector<int>>
  DistributedTree<Value>::query(const QueryType &queries)
  {
    Kokkos::View<int *, Kokkos::HostSpace> offsets("offsets", 0);
    Kokkos::View<Kokkos::pair<unsigned int, unsigned int> *, Kokkos::HostSpace>
      indices_ranks("indices_ranks", 0);
    if constexpr (QueryType::is_nearest)
      {
        distributed_tree.query(
          Kokkos::DefaultHostExecutionSpace{},
          queries,
          ArborX::Experimental::declare_callback_constrained(callback),
          indices_ranks,
          offsets);
      }
    else
      {
        distributed_tree.query(Kokkos::DefaultHostExecutionSpace{},
                               queries,
                               callback,
                               indices_ranks,
                               offsets);
      }

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

#  endif
} // namespace ArborXWrappers

DEAL_II_NAMESPACE_CLOSE

#else

// Make sure the scripts that create the C++20 module input files have
// something to latch on if the preprocessor #ifdef above would
// otherwise lead to an empty content of the file.
DEAL_II_NAMESPACE_OPEN
DEAL_II_NAMESPACE_CLOSE

#endif
#endif
