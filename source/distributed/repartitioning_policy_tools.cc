// ---------------------------------------------------------------------
//
// Copyright (C) 2021 by the deal.II authors
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


#include <deal.II/base/mpi_compute_index_owner_internal.h>
#include <deal.II/base/mpi_consensus_algorithms.h>

#include <deal.II/distributed/repartitioning_policy_tools.h>

#include <deal.II/multigrid/mg_transfer_global_coarsening.templates.h>

DEAL_II_NAMESPACE_OPEN


namespace RepartitioningPolicyTools
{
  namespace
  {
    template <int dim, int spacedim>
    void
    add_indices_recursevly_for_first_child_policy(
      const TriaIterator<CellAccessor<dim, spacedim>> &cell,
      const internal::CellIDTranslator<dim> &          cell_id_translator,
      IndexSet &                                       is_fine)
    {
      is_fine.add_index(cell_id_translator.translate(cell));

      if (cell->level() > 0 &&
          (cell->index() % GeometryInfo<dim>::max_children_per_cell) == 0)
        add_indices_recursevly_for_first_child_policy(cell->parent(),
                                                      cell_id_translator,
                                                      is_fine);
    }
  } // namespace

  template <int dim, int spacedim>
  LinearAlgebra::distributed::Vector<double>
  DefaultPolicy<dim, spacedim>::partition(
    const Triangulation<dim, spacedim> &) const
  {
    return {}; // nothing to do
  }



  template <int dim, int spacedim>
  FirstChildPolicy<dim, spacedim>::FirstChildPolicy(
    const Triangulation<dim, spacedim> &tria_fine)
    : n_coarse_cells(tria_fine.n_global_coarse_cells())
    , n_global_levels(tria_fine.n_global_levels())
  {
    Assert(
      tria_fine.all_reference_cells_are_hyper_cube(),
      ExcMessage(
        "FirstChildPolicy is only working for pure hex meshes at the moment."))

      const internal::CellIDTranslator<dim>
        cell_id_translator(n_coarse_cells, n_global_levels);
    is_level_partitions.set_size(cell_id_translator.size());

    for (const auto &cell : tria_fine.active_cell_iterators())
      if (cell->is_locally_owned())
        add_indices_recursevly_for_first_child_policy(cell,
                                                      cell_id_translator,
                                                      is_level_partitions);
  }



  template <int dim, int spacedim>
  LinearAlgebra::distributed::Vector<double>
  FirstChildPolicy<dim, spacedim>::partition(
    const Triangulation<dim, spacedim> &tria_coarse_in) const
  {
    const auto communicator = tria_coarse_in.get_communicator();

    const internal::CellIDTranslator<dim> cell_id_translator(n_coarse_cells,
                                                             n_global_levels);

    IndexSet is_coarse(cell_id_translator.size());

    for (const auto &cell : tria_coarse_in.active_cell_iterators())
      if (cell->is_locally_owned())
        is_coarse.add_index(cell_id_translator.translate(cell));

    std::vector<unsigned int> owning_ranks_of_coarse_cells(
      is_coarse.n_elements());
    {
      Utilities::MPI::internal::ComputeIndexOwner::ConsensusAlgorithmsPayload
        process(is_level_partitions,
                is_coarse,
                communicator,
                owning_ranks_of_coarse_cells,
                false);

      Utilities::MPI::ConsensusAlgorithms::Selector<
        std::pair<types::global_cell_index, types::global_cell_index>,
        unsigned int>
        consensus_algorithm(process, communicator);
      consensus_algorithm.run();
    }

    const auto tria =
      dynamic_cast<const parallel::TriangulationBase<dim, spacedim> *>(
        &tria_coarse_in);

    Assert(tria, ExcNotImplemented());

    LinearAlgebra::distributed::Vector<double> partition(
      tria->global_active_cell_index_partitioner().lock());

    for (const auto &cell : tria_coarse_in.active_cell_iterators())
      if (cell->is_locally_owned())
        partition[cell->global_active_cell_index()] =
          owning_ranks_of_coarse_cells[is_coarse.index_within_set(
            cell_id_translator.translate(cell))];

    return partition;
  }


  template <int dim, int spacedim>
  MinimalGranularityPolicy<dim, spacedim>::MinimalGranularityPolicy(
    const unsigned int n_min_cells)
    : n_min_cells(n_min_cells)
  {}



  template <int dim, int spacedim>
  LinearAlgebra::distributed::Vector<double>
  MinimalGranularityPolicy<dim, spacedim>::partition(
    const Triangulation<dim, spacedim> &tria_in) const
  {
    const auto tria =
      dynamic_cast<const parallel::TriangulationBase<dim, spacedim> *>(
        &tria_in);

    Assert(tria, ExcNotImplemented());

    LinearAlgebra::distributed::Vector<double> partition(
      tria->global_active_cell_index_partitioner().lock());

    // step 1) check if all processes have enough cells

    unsigned int n_locally_owned_active_cells = 0;
    for (const auto &cell : tria_in.active_cell_iterators())
      if (cell->is_locally_owned())
        ++n_locally_owned_active_cells;

    const auto comm = tria_in.get_communicator();

    if (Utilities::MPI::min(n_locally_owned_active_cells, comm) >= n_min_cells)
      return {}; // all processes have enough cells

    // step 2) there are processes which do not have enough cells so that
    // a repartitioning kicks in with the aim that all processes that own
    // cells have at least the specified number of cells

    const unsigned int n_global_active_cells = tria_in.n_global_active_cells();

    const unsigned int n_partitions =
      std::max<unsigned int>(1,
                             std::min(n_global_active_cells / n_min_cells,
                                      Utilities::MPI::n_mpi_processes(comm)));

    const unsigned int min_cells = n_global_active_cells / n_partitions;

    const auto convert = [&](const unsigned int i) {
      // determine the owner of a given index: we try to assign each process
      // the same number of cells; if there is a remainder, we assign the
      // first processes an additional cell (so that the difference of number of
      // locally owned cells is never larger than one between processes).
      const unsigned int n_partitions_with_additional_cell =
        n_global_active_cells - min_cells * n_partitions;

      const unsigned int rank =
        (i < (min_cells + 1) * n_partitions_with_additional_cell) ?
          (i / (min_cells + 1)) :
          ((i - n_partitions_with_additional_cell) / min_cells);

      AssertIndexRange(rank, n_partitions);

      return rank;
    };

    for (const auto i : partition.locally_owned_elements())
      partition[i] = convert(i);

    return partition;
  }


} // namespace RepartitioningPolicyTools



/*-------------- Explicit Instantiations -------------------------------*/
#include "repartitioning_policy_tools.inst"

DEAL_II_NAMESPACE_CLOSE
