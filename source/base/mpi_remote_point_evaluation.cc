// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2021 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#include <deal.II/base/config.h>

#include <deal.II/base/bounding_box.h>
#include <deal.II/base/mpi_consensus_algorithms.h>
#include <deal.II/base/mpi_remote_point_evaluation.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/mapping.h>

#include <deal.II/grid/filtered_iterator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/grid_tools_cache.h>
#include <deal.II/grid/tria.h>

DEAL_II_NAMESPACE_OPEN


namespace Utilities
{
  namespace MPI
  {
    template <int dim, int spacedim>
    RemotePointEvaluation<dim, spacedim>::AdditionalData::AdditionalData(
      const double                              tolerance,
      const bool                                enforce_unique_mapping,
      const unsigned int                        rtree_level,
      const std::function<std::vector<bool>()> &marked_vertices)
      : tolerance(tolerance)
      , enforce_unique_mapping(enforce_unique_mapping)
      , rtree_level(rtree_level)
      , marked_vertices(marked_vertices)
    {}



    template <int dim, int spacedim>
    RemotePointEvaluation<dim, spacedim>::RemotePointEvaluation(
      const AdditionalData &additional_data)
      : additional_data(additional_data)
      , ready_flag(false)
    {}



    template <int dim, int spacedim>
    RemotePointEvaluation<dim, spacedim>::RemotePointEvaluation(
      const double                              tolerance,
      const bool                                enforce_unique_mapping,
      const unsigned int                        rtree_level,
      const std::function<std::vector<bool>()> &marked_vertices)
      : additional_data(tolerance,
                        enforce_unique_mapping,
                        rtree_level,
                        marked_vertices)
      , ready_flag(false)
    {}



    template <int dim, int spacedim>
    RemotePointEvaluation<dim, spacedim>::~RemotePointEvaluation()
    {
      if (tria_signal.connected())
        tria_signal.disconnect();
    }



    template <int dim, int spacedim>
    void
    RemotePointEvaluation<dim, spacedim>::reinit(
      const std::vector<Point<spacedim>> &points,
      const Triangulation<dim, spacedim> &tria,
      const Mapping<dim, spacedim>       &mapping)
    {
      const GridTools::Cache<dim, spacedim> cache(tria, mapping);

      this->reinit(cache, points);
    }



    template <int dim, int spacedim>
    void
    RemotePointEvaluation<dim, spacedim>::reinit(
      const GridTools::Cache<dim, spacedim> &cache,
      const std::vector<Point<spacedim>>    &points)
    {
#ifndef DEAL_II_WITH_MPI
      Assert(false, ExcNeedsMPI());
      (void)cache;
      (void)points;
#else
      if (tria_signal.connected())
        tria_signal.disconnect();

      tria_signal = cache.get_triangulation().signals.any_change.connect(
        [&]() { this->ready_flag = false; });

      // compress r-tree to a minimal set of bounding boxes
      std::vector<std::vector<BoundingBox<spacedim>>> global_bboxes;
      global_bboxes.emplace_back(
        extract_rtree_level(cache.get_locally_owned_cell_bounding_boxes_rtree(),
                            additional_data.rtree_level));

      const auto data =
        GridTools::internal::distributed_compute_point_locations(
          cache,
          points,
          global_bboxes,
          additional_data.marked_vertices ? additional_data.marked_vertices() :
                                            std::vector<bool>(),
          additional_data.tolerance,
          true,
          additional_data.enforce_unique_mapping);

      this->reinit(data, cache.get_triangulation(), cache.get_mapping());
#endif
    }



    template <int dim, int spacedim>
    void
    RemotePointEvaluation<dim, spacedim>::reinit(
      const GridTools::internal::
        DistributedComputePointLocationsInternal<dim, spacedim> &data,
      const Triangulation<dim, spacedim>                        &tria,
      const Mapping<dim, spacedim>                              &mapping)
    {
      this->tria    = &tria;
      this->mapping = &mapping;

      this->recv_ranks = data.recv_ranks;
      this->recv_ptrs  = data.recv_ptrs;

      this->send_ranks = data.send_ranks;
      this->send_ptrs  = data.send_ptrs;

      this->recv_permutation = {};
      this->recv_permutation.resize(data.recv_components.size());
      this->point_ptrs.assign(data.n_searched_points + 1, 0);
      for (unsigned int i = 0; i < data.recv_components.size(); ++i)
        {
          AssertIndexRange(std::get<2>(data.recv_components[i]),
                           this->recv_permutation.size());
          this->recv_permutation[std::get<2>(data.recv_components[i])] = i;

          AssertIndexRange(std::get<1>(data.recv_components[i]) + 1,
                           this->point_ptrs.size());
          this->point_ptrs[std::get<1>(data.recv_components[i]) + 1]++;
        }

      std::pair<unsigned int, unsigned int> n_owning_processes_default{
        numbers::invalid_unsigned_int, 0};
      std::pair<unsigned int, unsigned int> n_owning_processes_local =
        n_owning_processes_default;

      for (unsigned int i = 0; i < data.n_searched_points; ++i)
        {
          std::get<0>(n_owning_processes_local) =
            std::min(std::get<0>(n_owning_processes_local),
                     this->point_ptrs[i + 1]);
          std::get<1>(n_owning_processes_local) =
            std::max(std::get<1>(n_owning_processes_local),
                     this->point_ptrs[i + 1]);

          this->point_ptrs[i + 1] += this->point_ptrs[i];
        }

      const auto n_owning_processes_global =
        Utilities::MPI::all_reduce<std::pair<unsigned int, unsigned int>>(
          n_owning_processes_local,
          tria.get_mpi_communicator(),
          [&](const auto &a,
              const auto &b) -> std::pair<unsigned int, unsigned int> {
            if (a == n_owning_processes_default)
              return b;

            if (b == n_owning_processes_default)
              return a;

            return std::pair<unsigned int, unsigned int>{
              std::min(std::get<0>(a), std::get<0>(b)),
              std::max(std::get<1>(a), std::get<1>(b))};
          });

      if (n_owning_processes_global == n_owning_processes_default)
        {
          unique_mapping        = true;
          all_points_found_flag = true;
        }
      else
        {
          unique_mapping = (std::get<0>(n_owning_processes_global) == 1) &&
                           (std::get<1>(n_owning_processes_global) == 1);
          all_points_found_flag = std::get<0>(n_owning_processes_global) > 0;
        }

      Assert(additional_data.enforce_unique_mapping == false || unique_mapping,
             ExcInternalError());

      cell_data        = std::make_unique<CellData>(tria);
      send_permutation = {};

      std::pair<int, int> dummy{-1, -1};
      for (const auto &i : data.send_components)
        {
          if (dummy != std::get<0>(i))
            {
              dummy = std::get<0>(i);
              cell_data->cells.emplace_back(dummy);
              cell_data->reference_point_ptrs.emplace_back(
                cell_data->reference_point_values.size());
            }

          cell_data->reference_point_values.emplace_back(std::get<3>(i));
          send_permutation.emplace_back(std::get<5>(i));
        }

      cell_data->reference_point_ptrs.emplace_back(
        cell_data->reference_point_values.size());

      unsigned int max_size_recv = 0;
      for (unsigned int i = 0; i < recv_ranks.size(); ++i)
        max_size_recv =
          std::max(max_size_recv, recv_ptrs[i + 1] - recv_ptrs[i]);

      unsigned int max_size_send = 0;
      for (unsigned int i = 0; i < send_ranks.size(); ++i)
        max_size_send =
          std::max(max_size_send, send_ptrs[i + 1] - send_ptrs[i]);

      this->buffer_size_with_sorting =
        std::max(send_permutation.size() * 2 + max_size_recv,
                 point_ptrs.back() + send_permutation.size() + max_size_send);

      this->buffer_size_without_sorting = send_permutation.size();

      // invert permutation matrices
      recv_permutation_inv.resize(recv_permutation.size());
      for (unsigned int c = 0; c < recv_permutation.size(); ++c)
        recv_permutation_inv[recv_permutation[c]] = c;

      send_permutation_inv.resize(send_permutation.size());
      for (unsigned int c = 0; c < send_permutation.size(); ++c)
        send_permutation_inv[send_permutation[c]] = c;

      this->ready_flag = true;
    }



    template <int dim, int spacedim>
    RemotePointEvaluation<dim, spacedim>::CellData::CellData(
      const Triangulation<dim, spacedim> &triangulation)
      : triangulation(triangulation)
    {}



    template <int dim, int spacedim>
    std_cxx20::ranges::iota_view<unsigned int, unsigned int>
    RemotePointEvaluation<dim, spacedim>::CellData::cell_indices() const
    {
      return std_cxx20::ranges::iota_view<unsigned int, unsigned int>(
        0, static_cast<unsigned int>(cells.size()));
    }



    template <int dim, int spacedim>
    typename Triangulation<dim, spacedim>::active_cell_iterator
    RemotePointEvaluation<dim, spacedim>::CellData::get_active_cell_iterator(
      const unsigned int cell) const
    {
      AssertIndexRange(cell, cells.size());
      return {&triangulation, cells[cell].first, cells[cell].second};
    }



    template <int dim, int spacedim>
    ArrayView<const Point<dim>>
    RemotePointEvaluation<dim, spacedim>::CellData::get_unit_points(
      const unsigned int cell) const
    {
      AssertIndexRange(cell, cells.size());
      return {reference_point_values.data() + reference_point_ptrs[cell],
              reference_point_ptrs[cell + 1] - reference_point_ptrs[cell]};
    }



    template <int dim, int spacedim>
    const typename RemotePointEvaluation<dim, spacedim>::CellData &
    RemotePointEvaluation<dim, spacedim>::get_cell_data() const
    {
      return *cell_data;
    }



    template <int dim, int spacedim>
    const std::vector<unsigned int> &
    RemotePointEvaluation<dim, spacedim>::get_point_ptrs() const
    {
      return point_ptrs;
    }



    template <int dim, int spacedim>
    bool
    RemotePointEvaluation<dim, spacedim>::is_map_unique() const
    {
      return unique_mapping;
    }



    template <int dim, int spacedim>
    bool
    RemotePointEvaluation<dim, spacedim>::all_points_found() const
    {
      return all_points_found_flag;
    }



    template <int dim, int spacedim>
    bool
    RemotePointEvaluation<dim, spacedim>::point_found(
      const unsigned int i) const
    {
      AssertIndexRange(i, point_ptrs.size() - 1);

      if (all_points_found_flag)
        return true;
      else
        return (point_ptrs[i + 1] - point_ptrs[i]) > 0;
    }



    template <int dim, int spacedim>
    const Triangulation<dim, spacedim> &
    RemotePointEvaluation<dim, spacedim>::get_triangulation() const
    {
      return *tria;
    }



    template <int dim, int spacedim>
    const Mapping<dim, spacedim> &
    RemotePointEvaluation<dim, spacedim>::get_mapping() const
    {
      return *mapping;
    }



    template <int dim, int spacedim>
    bool
    RemotePointEvaluation<dim, spacedim>::is_ready() const
    {
      return ready_flag;
    }



    template <int dim, int spacedim>
    const std::vector<unsigned int> &
    RemotePointEvaluation<dim, spacedim>::get_send_permutation() const
    {
      return send_permutation;
    }



    template <int dim, int spacedim>
    const std::vector<unsigned int> &
    RemotePointEvaluation<dim, spacedim>::get_inverse_recv_permutation() const
    {
      return recv_permutation_inv;
    }

  } // end of namespace MPI
} // end of namespace Utilities

#include "base/mpi_remote_point_evaluation.inst"

DEAL_II_NAMESPACE_CLOSE
