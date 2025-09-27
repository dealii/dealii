// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2021 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#ifndef dealii_distributed_repartitioning_policy_tools_h
#define dealii_distributed_repartitioning_policy_tools_h

#include <deal.II/grid/tria.h>

#include <deal.II/lac/la_parallel_vector.h>

DEAL_II_NAMESPACE_OPEN

/**
 * A namespace with repartitioning policies. These classes return vectors
 * of the new owners of the active locally owned and ghost cells of a
 * Triangulation object. The returned vectors can be used, e.g., in
 * TriangulationDescription::Utilities::create_description_from_triangulation()
 * to create a TriangulationDescription::Description based on a given
 * Triangulation and the predescribed partition, which can be used to
 * set up a parallel::fullydistributed::Triangulation objects.
 *
 * These policies can be also used in context of
 * MGTransferGlobalCoarseningTools::create_geometric_coarsening_sequence() to
 * prescribe arbitrary partitioning in multgrid levels of global coarsening
 * multigrid schmeme.
 */
namespace RepartitioningPolicyTools
{
  /**
   * The base class for repartitioning policies.
   *
   * Used in
   * MGTransferGlobalCoarseningTools::create_geometric_coarsening_sequence().
   * See the description of RepartitioningPolicyTools for more information.
   */
  template <int dim, int spacedim = dim>
  class Base : public EnableObserverPointer
  {
  public:
    /**
     * Return a vector of the new owners of the active locally owned and ghost
     * cells.
     */
    virtual LinearAlgebra::distributed::Vector<double>
    partition(const Triangulation<dim, spacedim> &tria_coarse_in) const = 0;

    /**
     * Destructor.
     */
    virtual ~Base() = default;
  };

  /**
   * A dummy policy that simply returns an empty vector, which is interpreted
   * in MGTransferGlobalCoarseningTools::create_geometric_coarsening_sequence()
   * in a way that the triangulation is not repartitioned.
   */
  template <int dim, int spacedim = dim>
  class DefaultPolicy : public Base<dim, spacedim>
  {
  public:
    /**
     * Constructor.
     *
     * @param tighten allows to renumber of subdomains so that empty ranks are
     *   positioned at the end.
     */
    DefaultPolicy(const bool tighten = false);

    virtual LinearAlgebra::distributed::Vector<double>
    partition(
      const Triangulation<dim, spacedim> &tria_coarse_in) const override;

  private:
    const bool tighten;
  };

  /**
   * A policy that partitions coarse grids based on a base triangulation
   * according to a first-child policy. The triangulation to be partitioned
   * should be able to be obtained by a sequence of (global) coarsening steps.
   */
  template <int dim, int spacedim = dim>
  class FirstChildPolicy : public Base<dim, spacedim>
  {
  public:
    /**
     * Constructor taking the base (fine) triangulation.
     */
    FirstChildPolicy(const Triangulation<dim, spacedim> &tria_fine);

    virtual LinearAlgebra::distributed::Vector<double>
    partition(
      const Triangulation<dim, spacedim> &tria_coarse_in) const override;

  private:
    /**
     * Number of coarse cells.
     */
    const unsigned int n_coarse_cells;

    /**
     * Number of global levels.
     */
    const unsigned int n_global_levels;

    /**
     * Index set constructed from the triangulation passed to the constructor.
     * It contains all the cells that would be owned by the current process
     * if the levels would be partitioned according to a first-child policy.
     */
    IndexSet is_level_partitions;
  };

  /**
   * A policy that allows to specify a minimal number of cells per
   * process. If a threshold is reached, processes might be left
   * without cells. The cells will be distributed evenly among the
   * remaining processes.
   */
  template <int dim, int spacedim = dim>
  class MinimalGranularityPolicy : public Base<dim, spacedim>
  {
  public:
    /**
     * Constructor taking the minimum number of cells per process.
     */
    MinimalGranularityPolicy(const unsigned int n_min_cells);

    virtual LinearAlgebra::distributed::Vector<double>
    partition(const Triangulation<dim, spacedim> &tria_in) const override;

  private:
    /**
     * Minimum number of cells per process.
     */
    const unsigned int n_min_cells;
  };

  /**
   * A policy that allows to specify a weight of each cell. The underlying
   * algorithm will try to distribute the weights equally among the processes.
   */
  template <int dim, int spacedim = dim>
  class CellWeightPolicy : public Base<dim, spacedim>
  {
  public:
    /**
     * Constructor taking a function that gives a weight to each cell.
     */
    CellWeightPolicy(
      const std::function<unsigned int(
        const typename Triangulation<dim, spacedim>::cell_iterator &,
        const CellStatus)> &weighting_function);

    virtual LinearAlgebra::distributed::Vector<double>
    partition(const Triangulation<dim, spacedim> &tria_in) const override;

  private:
    /**
     * A function that gives a weight to each cell.
     */
    const std::function<
      unsigned int(const typename Triangulation<dim, spacedim>::cell_iterator &,
                   const CellStatus)>
      weighting_function;
  };

} // namespace RepartitioningPolicyTools

DEAL_II_NAMESPACE_CLOSE

#endif
