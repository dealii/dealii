// ---------------------------------------------------------------------
//
// Copyright (C) 2021 - 2022 by the deal.II authors
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
#include <deal.II/distributed/tria_base.h>

#include <deal.II/grid/cell_id_translator.h>
#include <deal.II/grid/filtered_iterator.h>

DEAL_II_NAMESPACE_OPEN


namespace RepartitioningPolicyTools
{
  namespace
  {
    template <int dim, int spacedim>
    void
    add_indices_recursively_for_first_child_policy(
      const TriaIterator<CellAccessor<dim, spacedim>> &cell,
      const internal::CellIDTranslator<dim> &          cell_id_translator,
      IndexSet &                                       is_fine)
    {
      is_fine.add_index(cell_id_translator.translate(cell));

      if (cell->level() > 0 && cell->parent()->child(0) == cell)
        add_indices_recursively_for_first_child_policy(cell->parent(),
                                                       cell_id_translator,
                                                       is_fine);
    }
  } // namespace


  template <int dim, int spacedim>
  DefaultPolicy<dim, spacedim>::DefaultPolicy(const bool tighten)
    : tighten(tighten)
  {}

  template <int dim, int spacedim>
  LinearAlgebra::distributed::Vector<double>
  DefaultPolicy<dim, spacedim>::partition(
    const Triangulation<dim, spacedim> &tria_in) const
  {
    if (tighten == false)
      return {}; // nothing to do

    const auto tria =
      dynamic_cast<const parallel::TriangulationBase<dim, spacedim> *>(
        &tria_in);

    if (tria == nullptr)
      return {}; // nothing to do, since serial triangulation

#ifndef DEAL_II_WITH_MPI
    Assert(false, ExcNeedsMPI());
    return {};
#else

    const auto comm = tria->get_communicator();

    const unsigned int process_has_active_locally_owned_cells =
      tria->n_locally_owned_active_cells() > 0;
    const unsigned int n_processes_with_active_locally_owned_cells =
      Utilities::MPI::sum(process_has_active_locally_owned_cells, comm);

    if (n_processes_with_active_locally_owned_cells ==
        Utilities::MPI::n_mpi_processes(comm))
      return {}; // nothing to do, since all processes have cells

    unsigned int offset = 0;

    const int ierr = MPI_Exscan(&process_has_active_locally_owned_cells,
                                &offset,
                                1,
                                Utilities::MPI::mpi_type_id_for_type<decltype(
                                  process_has_active_locally_owned_cells)>,
                                MPI_SUM,
                                comm);
    AssertThrowMPI(ierr);

    LinearAlgebra::distributed::Vector<double> partition(
      tria->global_active_cell_index_partitioner().lock());

    partition = offset;

    return partition;
#endif
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
        add_indices_recursively_for_first_child_policy(cell,
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
        std::vector<
          std::pair<types::global_cell_index, types::global_cell_index>>,
        std::vector<unsigned int>>
        consensus_algorithm(process, communicator);
      consensus_algorithm.run();
    }

    const auto tria =
      dynamic_cast<const parallel::TriangulationBase<dim, spacedim> *>(
        &tria_coarse_in);

    Assert(tria, ExcNotImplemented());

    LinearAlgebra::distributed::Vector<double> partition(
      tria->global_active_cell_index_partitioner().lock());

    for (const auto &cell : tria_coarse_in.active_cell_iterators() |
                              IteratorFilters::LocallyOwnedCell())
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

    const unsigned int n_locally_owned_active_cells =
      std::count_if(tria_in.begin_active(),
                    typename Triangulation<dim, spacedim>::active_cell_iterator(
                      tria_in.end()),
                    [](const auto &cell) { return cell.is_locally_owned(); });

    const auto comm = tria_in.get_communicator();

    if (Utilities::MPI::min(n_locally_owned_active_cells, comm) >= n_min_cells)
      return {}; // all processes have enough cells

    // step 2) there are processes which do not have enough cells so that
    // a repartitioning kicks in with the aim that all processes that own
    // cells have at least the specified number of cells

    const types::global_cell_index n_global_active_cells =
      tria_in.n_global_active_cells();

    const unsigned int n_partitions =
      std::max<unsigned int>(1,
                             std::min<types::global_cell_index>(
                               n_global_active_cells / n_min_cells,
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



  template <int dim, int spacedim>
  CellWeightPolicy<dim, spacedim>::CellWeightPolicy(
    const std::function<
      unsigned int(const typename Triangulation<dim, spacedim>::cell_iterator &,
                   const typename Triangulation<dim, spacedim>::CellStatus)>
      &weighting_function)
    : weighting_function(weighting_function)
  {}



  template <int dim, int spacedim>
  LinearAlgebra::distributed::Vector<double>
  CellWeightPolicy<dim, spacedim>::partition(
    const Triangulation<dim, spacedim> &tria_in) const
  {
#ifndef DEAL_II_WITH_MPI
    (void)tria_in;
    return {};
#else

    const auto tria =
      dynamic_cast<const parallel::TriangulationBase<dim, spacedim> *>(
        &tria_in);

    Assert(tria, ExcNotImplemented());

    const auto partitioner =
      tria->global_active_cell_index_partitioner().lock();

    std::vector<unsigned int> weights(partitioner->local_size());

    const auto mpi_communicator = tria_in.get_communicator();
    const auto n_subdomains = Utilities::MPI::n_mpi_processes(mpi_communicator);

    // determine weight of each cell
    for (const auto &cell :
         tria->active_cell_iterators() | IteratorFilters::LocallyOwnedCell())
      weights[partitioner->global_to_local(cell->global_active_cell_index())] =
        weighting_function(
          cell, Triangulation<dim, spacedim>::CellStatus::CELL_PERSIST);

    // determine weight of all the cells locally owned by this process
    std::uint64_t process_local_weight = 0;
    for (const auto &weight : weights)
      process_local_weight += weight;

    // determine partial sum of weights of this process
    std::uint64_t process_local_weight_offset = 0;

    int ierr = MPI_Exscan(
      &process_local_weight,
      &process_local_weight_offset,
      1,
      Utilities::MPI::mpi_type_id_for_type<decltype(process_local_weight)>,
      MPI_SUM,
      tria->get_communicator());
    AssertThrowMPI(ierr);

    // total weight of all processes
    std::uint64_t total_weight =
      process_local_weight_offset + process_local_weight;

    ierr =
      MPI_Bcast(&total_weight,
                1,
                Utilities::MPI::mpi_type_id_for_type<decltype(total_weight)>,
                n_subdomains - 1,
                mpi_communicator);
    AssertThrowMPI(ierr);

    // setup partition
    LinearAlgebra::distributed::Vector<double> partition(partitioner);

    for (std::uint64_t i = 0, weight = process_local_weight_offset;
         i < partition.locally_owned_size();
         weight += weights[i], ++i)
      partition.local_element(i) =
        static_cast<double>(weight * n_subdomains / total_weight);

    return partition;
#endif
  }


} // namespace RepartitioningPolicyTools



/*-------------- Explicit Instantiations -------------------------------*/
#include "repartitioning_policy_tools.inst"

DEAL_II_NAMESPACE_CLOSE
