// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2021 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


#include <deal.II/base/mpi.h>

#include <deal.II/distributed/repartitioning_policy_tools.h>
#include <deal.II/distributed/tria_base.h>

#include <deal.II/fe/mapping_q1.h>

#include <deal.II/grid/cell_id_translator.h>
#include <deal.II/grid/filtered_iterator.h>

#include <deal.II/multigrid/mg_transfer_global_coarsening.templates.h>

#include <deal.II/numerics/vector_tools.h>

DEAL_II_NAMESPACE_OPEN


namespace RepartitioningPolicyTools
{
  namespace
  {
    template <int dim, int spacedim>
    void
    add_indices_recursively_for_first_child_policy(
      const TriaIterator<CellAccessor<dim, spacedim>> &cell,
      const internal::CellIDTranslator<dim>           &cell_id_translator,
      IndexSet                                        &is_fine)
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

    const auto comm = tria->get_mpi_communicator();

    const unsigned int process_has_active_locally_owned_cells =
      tria->n_locally_owned_active_cells() > 0;
    const unsigned int n_processes_with_active_locally_owned_cells =
      Utilities::MPI::sum(process_has_active_locally_owned_cells, comm);

    if (n_processes_with_active_locally_owned_cells ==
        Utilities::MPI::n_mpi_processes(comm))
      return {}; // nothing to do, since all processes have cells

    unsigned int offset = 0;

    const int ierr =
      MPI_Exscan(&process_has_active_locally_owned_cells,
                 &offset,
                 1,
                 Utilities::MPI::mpi_type_id_for_type<
                   decltype(process_has_active_locally_owned_cells)>,
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
        "FirstChildPolicy is only working for pure hex meshes at the moment."));

    const internal::CellIDTranslator<dim> cell_id_translator(n_coarse_cells,
                                                             n_global_levels);
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
    const auto communicator = tria_coarse_in.get_mpi_communicator();

    const internal::CellIDTranslator<dim> cell_id_translator(n_coarse_cells,
                                                             n_global_levels);

    IndexSet is_coarse(cell_id_translator.size());

    for (const auto &cell : tria_coarse_in.active_cell_iterators())
      if (cell->is_locally_owned())
        is_coarse.add_index(cell_id_translator.translate(cell));

    const std::vector<unsigned int> owning_ranks_of_coarse_cells =
      Utilities::MPI::compute_index_owner(is_level_partitions,
                                          is_coarse,
                                          communicator);

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

    const auto comm = tria_in.get_mpi_communicator();

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
                   const CellStatus)> &weighting_function)
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

    std::vector<unsigned int> weights(partitioner->locally_owned_size());

    const auto mpi_communicator = tria_in.get_mpi_communicator();
    const auto n_subdomains = Utilities::MPI::n_mpi_processes(mpi_communicator);

    // determine weight of each cell
    for (const auto &cell :
         tria->active_cell_iterators() | IteratorFilters::LocallyOwnedCell())
      weights[partitioner->global_to_local(cell->global_active_cell_index())] =
        weighting_function(cell, CellStatus::cell_will_persist);

    // determine weight of all the cells locally owned by this process
    std::uint64_t process_local_weight = 0;
    for (const auto &weight : weights)
      process_local_weight += weight;

    // determine partial sum of weights of this process, as well as the total
    // weight
    const auto [process_local_weight_offset, total_weight] =
      Utilities::MPI::partial_and_total_sum(process_local_weight,
                                            tria->get_mpi_communicator());

    // set up partition
    LinearAlgebra::distributed::Vector<double> partition(partitioner);

    for (std::uint64_t i = 0, weight = process_local_weight_offset;
         i < partition.locally_owned_size();
         weight += weights[i], ++i)
      partition.local_element(i) =
        static_cast<double>(weight * n_subdomains / total_weight);

    return partition;
#endif
  }


  template <int dim, int spacedim>
  ImmersedMeshPolicy<dim, spacedim>::AdditionalData::AdditionalData()
    : reduction_type(ReductionType::highest_count)
    , immersed_identification(true)
    , n_samples(1)
  {}



  template <int dim, int spacedim>
  ImmersedMeshPolicy<dim, spacedim>::ImmersedMeshPolicy(
    const Triangulation<dim, spacedim> &tria_background,
    const AdditionalData               &data)
    : tria_background(&tria_background)
    , dof_handler_background(nullptr)
    , data(data)
  {}



  template <int dim, int spacedim>
  ImmersedMeshPolicy<dim, spacedim>::ImmersedMeshPolicy(
    const DoFHandler<dim, spacedim> &dof_handler_background,
    const AdditionalData            &data)
    : tria_background(&dof_handler_background.get_triangulation())
    , dof_handler_background(&dof_handler_background)
    , data(data)
  {}



  template <int dim, int spacedim>
  LinearAlgebra::distributed::Vector<double>
  ImmersedMeshPolicy<dim, spacedim>::partition(
    const Triangulation<dim, spacedim> &tria_immersed) const
  {
    std::vector<std::vector<unsigned int>> cell_ranks(
      tria_immersed.n_active_cells());

    MappingQ1<dim, spacedim> mapping_default;

    const Mapping<dim, spacedim> &mapping_background =
      data.mapping_background ? *data.mapping_background : mapping_default;
    const Mapping<dim, spacedim> &mapping_immersed =
      data.mapping_immersed ? *data.mapping_immersed : mapping_default;

    if ((dof_handler_background == nullptr) && data.immersed_identification)
      {
        std::vector<Point<spacedim>> points;

        Quadrature<dim> quadrature;

        if (data.n_samples == 1)
          quadrature = QGauss<dim>(1);
        else
          quadrature = QGaussLobatto<dim>(data.n_samples);

        for (const auto &cell : tria_immersed.active_cell_iterators())
          if (cell->is_locally_owned())
            for (const auto &p : quadrature.get_points())
              points.push_back(
                mapping_immersed.transform_unit_to_real_cell(cell, p));

        Utilities::MPI::RemotePointEvaluation<dim, spacedim> rpe;
        rpe.reinit(points, *tria_background, mapping_background);

        const auto evaluate_function = [&](const ArrayView<double> &values,
                                           const auto              &cell_data) {
          for (const auto cell : cell_data.cell_indices())
            {
              const auto unit_points = cell_data.get_unit_points(cell);
              const auto local_value = cell_data.get_data_view(cell, values);

              for (unsigned int q = 0; q < unit_points.size(); ++q)
                local_value[q] = Utilities::MPI::this_mpi_process(
                  tria_background->get_communicator());
            }
        };

        const std::vector<double> point_ranks =
          rpe.template evaluate_and_process<double>(evaluate_function);

        unsigned int counter = 0;
        for (const auto &cell : tria_immersed.active_cell_iterators())
          if (cell->is_locally_owned())
            {
              unsigned int start =
                rpe.get_point_ptrs()[counter * quadrature.size()];
              unsigned int end =
                rpe.get_point_ptrs()[(counter + 1) * quadrature.size()];

              for (unsigned int i = start; i < end; ++i)
                cell_ranks[cell->active_cell_index()].push_back(
                  static_cast<unsigned int>(point_ranks[i]));

              counter++;
            }
      }
    else
      {
        std::vector<Point<spacedim>> points;

        if (dof_handler_background == nullptr)
          {
            Quadrature<dim> quadrature;

            if (data.n_samples == 1)
              quadrature = QGauss<dim>(1);
            else
              quadrature = QGaussLobatto<dim>(data.n_samples);

            for (const auto &cell : tria_background->active_cell_iterators())
              if (cell->is_locally_owned())
                for (const auto &p : quadrature.get_points())
                  points.push_back(
                    mapping_background.transform_unit_to_real_cell(cell, p));
          }
        else
          {
            std::tie(points, std::ignore, std::ignore) =
              internal::collect_unconstrained_unique_support_points(
                *dof_handler_background,
                mapping_background,
                AffineConstraints<double>());
          }

        Utilities::MPI::RemotePointEvaluation<dim, spacedim> rpe;
        rpe.reinit(points, tria_immersed, mapping_immersed);

        std::vector<double> integration_values(
          points.size(),
          Utilities::MPI::this_mpi_process(
            tria_background->get_communicator()));

        const auto integration_function = [&](const auto &values,
                                              const auto &cell_data) {
          for (const auto cell : cell_data.cell_indices())
            {
              const auto unit_points = cell_data.get_unit_points(cell);
              const auto local_value = cell_data.get_data_view(cell, values);

              for (unsigned int q = 0; q < unit_points.size(); ++q)
                cell_ranks[cell_data.get_active_cell_iterator(cell)
                             ->active_cell_index()]
                  .push_back(static_cast<unsigned int>(local_value[q]));
            }
        };

        rpe.template process_and_evaluate<double>(integration_values,
                                                  integration_function);
      }

    const auto tria =
      dynamic_cast<const parallel::TriangulationBase<dim, spacedim> *>(
        &tria_immersed);

    Assert(tria, ExcNotImplemented());

    // 3) set partitioning
    LinearAlgebra::distributed::Vector<double> partition(
      tria->global_active_cell_index_partitioner().lock());


    const auto reduce = [this](const auto &vector) -> unsigned int {
      AssertThrow(!vector.empty(), ExcInternalError());

      if (data.reduction_type == ReductionType::smallest_rank)
        {
          unsigned int rank = numbers::invalid_unsigned_int;

          for (const auto rank_i : vector)
            rank = std::min<unsigned int>(rank, rank_i);

          return rank;
        }
      else if (data.reduction_type == ReductionType::highest_count)
        {
          std::map<unsigned int, unsigned int> rank_counter;

          for (const auto rank_i : vector)
            rank_counter[rank_i] = 0;
          for (const auto rank_i : vector)
            rank_counter[rank_i]++;

          const auto pr =
            std::max_element(rank_counter.begin(),
                             rank_counter.end(),
                             [](const auto &p1, const auto &p2) {
                               if (p1.second != p2.second)
                                 return p1.second < p2.second;

                               return p1.first > p2.first; // stable search
                             });

          AssertThrow(pr != rank_counter.end(), ExcInternalError());

          return pr->first;
        }
      else
        {
          AssertThrow(false, ExcNotImplemented());

          return 0;
        }
    };

    for (const auto &cell : tria_immersed.active_cell_iterators())
      if (cell->is_locally_owned())
        partition[cell->global_active_cell_index()] =
          reduce(cell_ranks[cell->active_cell_index()]);

    partition.update_ghost_values();

    return partition;
  }


} // namespace RepartitioningPolicyTools



/*-------------- Explicit Instantiations -------------------------------*/
#include "distributed/repartitioning_policy_tools.inst"

DEAL_II_NAMESPACE_CLOSE
