// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2020 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#ifndef dealii_mg_transfer_global_coarsening_h
#define dealii_mg_transfer_global_coarsening_h

#include <deal.II/base/config.h>

#include <deal.II/grid/tria.h>

#include <deal.II/multigrid/mg_transfer_matrix_free.h>


DEAL_II_NAMESPACE_OPEN

// Forward declarations
#ifndef DOXYGEN
namespace RepartitioningPolicyTools
{
  template <int dim, int spacedim>
  class Base;
}
#endif

/**
 * Global coarsening utility functions.
 */
namespace MGTransferGlobalCoarseningTools
{
  /**
   * Common polynomial coarsening sequences.
   *
   * @note These polynomial coarsening sequences up to a degree of 9 are
   *   precompiled in MGTwoLevelTransfer. See also:
   *   MGTwoLevelTransfer::fast_polynomial_transfer_supported()
   */
  enum class PolynomialCoarseningSequenceType
  {
    /**
     * Half polynomial degree by integer division. For example, for degree=7
     * the following sequence would be obtained:: 7 -> 3 -> 1
     */
    bisect,
    /**
     * Decrease the polynomial degree by one. E.g., for degree=7 following
     * sequence would result: 7 -> 6 -> 5 -> 4 -> 3 -> 2 -> 1
     */
    decrease_by_one,
    /**
     * Decrease the polynomial degree to one. E.g., for degree=7 following
     * sequence would result: 7 -> 1
     */
    go_to_one
  };

  /**
   * For a given @p degree and polynomial coarsening sequence @p p_sequence,
   * determine the next coarser degree.
   */
  unsigned int
  create_next_polynomial_coarsening_degree(
    const unsigned int                      degree,
    const PolynomialCoarseningSequenceType &p_sequence);

  /**
   * For a given @p max_degree and polynomial coarsening sequence @p p_sequence,
   * determine the full sequence of polynomial degrees, sorted in ascending
   * order.
   */
  std::vector<unsigned int>
  create_polynomial_coarsening_sequence(
    const unsigned int                      max_degree,
    const PolynomialCoarseningSequenceType &p_sequence);

  /**
   * For a given triangulation @p tria, determine the geometric coarsening
   * sequence by repeated global coarsening of the provided triangulation.
   *
   * @note For convenience, a reference to the input triangulation is stored in
   *   the last entry of the return vector.
   * @note Currently, not implemented for parallel::fullydistributed::Triangulation.
   * @note The type of the returned triangulations is the same as of the input
   *   triangulation.
   */
  template <int dim, int spacedim>
  std::vector<std::shared_ptr<const Triangulation<dim, spacedim>>>
  create_geometric_coarsening_sequence(
    const Triangulation<dim, spacedim> &tria);

  /**
   * Similar to the above function but also taking a @p policy for
   * repartitioning the triangulations on the coarser levels. If
   * @p preserve_fine_triangulation is set, the input triangulation is not
   * altered,
   * else the triangulation is coarsened. If @p repartition_fine_triangulation
   * is set, the triangulation on the finest level is repartitioned as well. If
   * the flags are set to true/false, the input triangulation is simply used as
   * the finest triangulation.
   *
   * @note For convenience, a reference to the input triangulation is stored in
   *   the last entry of the return vector.
   * @note The type of the returned triangulations is
   *   parallel::fullydistributed::Triangulation.
   * @note Currently, only implemented for parallel::distributed::Triangulation.
   */
  template <int dim, int spacedim>
  std::vector<std::shared_ptr<const Triangulation<dim, spacedim>>>
  create_geometric_coarsening_sequence(
    Triangulation<dim, spacedim>                         &tria,
    const RepartitioningPolicyTools::Base<dim, spacedim> &policy,
    const bool preserve_fine_triangulation,
    const bool repartition_fine_triangulation);

  /**
   * Similar to the above function but taking in a constant version of
   * @p tria. As a consequence, it can not be used for coarsening directly,
   * so a temporary copy will be created internally.
   */
  template <int dim, int spacedim>
  std::vector<std::shared_ptr<const Triangulation<dim, spacedim>>>
  create_geometric_coarsening_sequence(
    const Triangulation<dim, spacedim>                   &tria,
    const RepartitioningPolicyTools::Base<dim, spacedim> &policy,
    const bool repartition_fine_triangulation = false);

} // namespace MGTransferGlobalCoarseningTools


template <int dim, typename Number>
using MGTransferMF DEAL_II_DEPRECATED_EARLY = MGTransferMatrixFree<dim, Number>;

template <int dim, typename Number>
using MGTransferBlockMF DEAL_II_DEPRECATED_EARLY =
  MGTransferBlockMatrixFree<dim, Number>;

template <int dim, typename VectorType>
using MGTransferGlobalCoarsening DEAL_II_DEPRECATED_EARLY =
  MGTransferMatrixFree<dim,
                       typename VectorType::value_type,
                       ::dealii::MemorySpace::Host>;

template <int dim, typename VectorType>
using MGTransferBlockGlobalCoarsening DEAL_II_DEPRECATED_EARLY =
  MGTransferBlockMatrixFree<dim, typename VectorType::value_type>;


DEAL_II_NAMESPACE_CLOSE

#endif
