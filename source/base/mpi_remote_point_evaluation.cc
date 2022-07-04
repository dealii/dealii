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
    RemotePointEvaluation<dim, spacedim>::RemotePointEvaluation(
      const double                              tolerance,
      const bool                                enforce_unique_mapping,
      const unsigned int                        rtree_level,
      const std::function<std::vector<bool>()> &marked_vertices)
      : tolerance(tolerance)
      , enforce_unique_mapping(enforce_unique_mapping)
      , rtree_level(rtree_level)
      , marked_vertices(marked_vertices)
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
      const Mapping<dim, spacedim> &      mapping)
    {
#ifndef DEAL_II_WITH_MPI
      Assert(false, ExcNeedsMPI());
      (void)points;
      (void)tria;
      (void)mapping;
#else
      if (tria_signal.connected())
        tria_signal.disconnect();

      tria_signal =
        tria.signals.any_change.connect([&]() { this->ready_flag = false; });

      this->tria    = &tria;
      this->mapping = &mapping;

      std::vector<BoundingBox<spacedim>> local_boxes;
      for (const auto &cell :
           tria.active_cell_iterators() | IteratorFilters::LocallyOwnedCell())
        local_boxes.push_back(mapping.get_bounding_box(cell));

      // create r-tree of bounding boxes
      const auto local_tree = pack_rtree(local_boxes);

      // compress r-tree to a minimal set of bounding boxes
      const auto local_reduced_box =
        extract_rtree_level(local_tree, rtree_level);

      // gather bounding boxes of other processes
      const auto global_bboxes =
        Utilities::MPI::all_gather(tria.get_communicator(), local_reduced_box);

      const GridTools::Cache<dim, spacedim> cache(tria, mapping);

      const auto data =
        GridTools::internal::distributed_compute_point_locations(
          cache,
          points,
          global_bboxes,
          marked_vertices ? marked_vertices() : std::vector<bool>(),
          tolerance,
          true,
          enforce_unique_mapping);

      this->recv_ranks = data.recv_ranks;
      this->recv_ptrs  = data.recv_ptrs;

      this->send_ranks = data.send_ranks;
      this->send_ptrs  = data.send_ptrs;

      this->recv_permutation = {};
      this->recv_permutation.resize(data.recv_components.size());
      this->point_ptrs.assign(points.size() + 1, 0);
      for (unsigned int i = 0; i < data.recv_components.size(); ++i)
        {
          AssertIndexRange(std::get<2>(data.recv_components[i]),
                           this->recv_permutation.size());
          this->recv_permutation[std::get<2>(data.recv_components[i])] = i;

          AssertIndexRange(std::get<1>(data.recv_components[i]) + 1,
                           this->point_ptrs.size());
          this->point_ptrs[std::get<1>(data.recv_components[i]) + 1]++;
        }

      std::tuple<unsigned int, unsigned int> n_owning_processes_default{
        numbers::invalid_unsigned_int, 0};
      std::tuple<unsigned int, unsigned int> n_owning_processes_local =
        n_owning_processes_default;

      for (unsigned int i = 0; i < points.size(); ++i)
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
        Utilities::MPI::all_reduce<std::tuple<unsigned int, unsigned int>>(
          n_owning_processes_local,
          tria.get_communicator(),
          [&](const auto &a,
              const auto &b) -> std::tuple<unsigned int, unsigned int> {
            if (a == n_owning_processes_default)
              return b;

            if (b == n_owning_processes_default)
              return a;

            return std::tuple<unsigned int, unsigned int>{
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

      Assert(enforce_unique_mapping == false || unique_mapping,
             ExcInternalError());

      cell_data        = {};
      send_permutation = {};

      std::pair<int, int> dummy{-1, -1};
      for (const auto &i : data.send_components)
        {
          if (dummy != std::get<0>(i))
            {
              dummy = std::get<0>(i);
              cell_data.cells.emplace_back(dummy);
              cell_data.reference_point_ptrs.emplace_back(
                cell_data.reference_point_values.size());
            }

          cell_data.reference_point_values.emplace_back(std::get<3>(i));
          send_permutation.emplace_back(std::get<5>(i));
        }

      cell_data.reference_point_ptrs.emplace_back(
        cell_data.reference_point_values.size());

      this->ready_flag = true;
#endif
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

  } // end of namespace MPI
} // end of namespace Utilities

#include "mpi_remote_point_evaluation.inst"

DEAL_II_NAMESPACE_CLOSE
