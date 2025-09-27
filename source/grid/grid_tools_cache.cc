// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2017 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#include <deal.II/base/bounding_box.h>
#include <deal.II/base/mpi_stub.h>

#include <deal.II/grid/filtered_iterator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/grid_tools_cache.h>

DEAL_II_NAMESPACE_OPEN

namespace GridTools
{
  template <int dim, int spacedim>
  Cache<dim, spacedim>::Cache(const Triangulation<dim, spacedim> &tria,
                              const Mapping<dim, spacedim>       &mapping)
    : update_flags(update_all)
    , tria(&tria)
    , mapping(&mapping)
  {
    tria_change_signal =
      tria.signals.any_change.connect([&]() { mark_for_update(update_all); });
  }



  template <int dim, int spacedim>
  Cache<dim, spacedim>::Cache(const Triangulation<dim, spacedim> &tria)
    : update_flags(update_all)
    , tria(&tria)
  {
    tria_change_signal =
      tria.signals.any_change.connect([&]() { mark_for_update(update_all); });

    // Allow users to set this class up with an empty Triangulation and no
    // Mapping argument by deferring Mapping assignment until after the
    // Triangulation exists.
    if (tria.get_reference_cells().size() == 0)
      {
        tria_create_signal = tria.signals.create.connect([&]() {
          Assert(tria.get_reference_cells().size() > 0, ExcInternalError());
          Assert(tria.get_reference_cells().size() == 1, ExcNotImplemented());
          mapping = &tria.get_reference_cells()[0]
                       .template get_default_linear_mapping<dim, spacedim>();
        });
      }
    else
      {
        Assert(tria.get_reference_cells().size() == 1, ExcNotImplemented());
        mapping = &tria.get_reference_cells()[0]
                     .template get_default_linear_mapping<dim, spacedim>();
      }
  }

  template <int dim, int spacedim>
  Cache<dim, spacedim>::~Cache()
  {
    if (tria_change_signal.connected())
      tria_change_signal.disconnect();
    if (tria_create_signal.connected())
      tria_create_signal.disconnect();
  }



  template <int dim, int spacedim>
  void
  Cache<dim, spacedim>::mark_for_update(const CacheUpdateFlags &flags)
  {
    update_flags |= flags;
  }



  template <int dim, int spacedim>
  const std::vector<
    std::set<typename Triangulation<dim, spacedim>::active_cell_iterator>> &
  Cache<dim, spacedim>::get_vertex_to_cell_map() const
  {
    // In the following, we will first check whether the data structure
    // in question needs to be updated (in which case we update it, and
    // reset the flag that indices that this needs to happen to zero), and
    // then return it. Make this thread-safe by using a mutex to guard
    // all of this:
    std::lock_guard<std::mutex> lock(vertex_to_cells_mutex);

    if (update_flags & update_vertex_to_cell_map)
      {
        vertex_to_cells = GridTools::vertex_to_cell_map(*tria);

        // Atomically clear the flag that indicates that this data member
        // needs to be updated:
        update_flags &= ~update_vertex_to_cell_map;
      }
    return vertex_to_cells;
  }



  template <int dim, int spacedim>
  const std::vector<std::vector<Tensor<1, spacedim>>> &
  Cache<dim, spacedim>::get_vertex_to_cell_centers_directions() const
  {
    // In the following, we will first check whether the data structure
    // in question needs to be updated (in which case we update it, and
    // reset the flag that indices that this needs to happen to zero), and
    // then return it. Make this thread-safe by using a mutex to guard
    // all of this:
    std::lock_guard<std::mutex> lock(vertex_to_cell_centers_mutex);

    if (update_flags & update_vertex_to_cell_centers_directions)
      {
        vertex_to_cell_centers = GridTools::vertex_to_cell_centers_directions(
          *tria, get_vertex_to_cell_map());

        // Atomically clear the flag that indicates that this data member
        // needs to be updated:
        update_flags &= ~update_vertex_to_cell_centers_directions;
      }
    return vertex_to_cell_centers;
  }



  template <int dim, int spacedim>
  const std::map<unsigned int, Point<spacedim>> &
  Cache<dim, spacedim>::get_used_vertices() const
  {
    // In the following, we will first check whether the data structure
    // in question needs to be updated (in which case we update it, and
    // reset the flag that indices that this needs to happen to zero), and
    // then return it. Make this thread-safe by using a mutex to guard
    // all of this:
    std::lock_guard<std::mutex> lock(used_vertices_mutex);

    if (update_flags & update_used_vertices)
      {
        used_vertices = GridTools::extract_used_vertices(*tria, *mapping);

        // Atomically clear the flag that indicates that this data member
        // needs to be updated:
        update_flags &= ~update_used_vertices;
      }
    return used_vertices;
  }



  template <int dim, int spacedim>
  const RTree<std::pair<Point<spacedim>, unsigned int>> &
  Cache<dim, spacedim>::get_used_vertices_rtree() const
  {
    // In the following, we will first check whether the data structure
    // in question needs to be updated (in which case we update it, and
    // reset the flag that indices that this needs to happen to zero), and
    // then return it. Make this thread-safe by using a mutex to guard
    // all of this:
    std::lock_guard<std::mutex> lock(used_vertices_rtree_mutex);

    if (update_flags & update_used_vertices_rtree)
      {
        const auto &used_vertices = get_used_vertices();
        std::vector<std::pair<Point<spacedim>, unsigned int>> vertices(
          used_vertices.size());
        unsigned int i = 0;
        for (const auto &it : used_vertices)
          vertices[i++] = std::make_pair(it.second, it.first);
        used_vertices_rtree = pack_rtree(vertices);

        // Atomically clear the flag that indicates that this data member
        // needs to be updated:
        update_flags &= ~update_used_vertices_rtree;
      }
    return used_vertices_rtree;
  }



  template <int dim, int spacedim>
  const RTree<
    std::pair<BoundingBox<spacedim>,
              typename Triangulation<dim, spacedim>::active_cell_iterator>> &
  Cache<dim, spacedim>::get_cell_bounding_boxes_rtree() const
  {
    // In the following, we will first check whether the data structure
    // in question needs to be updated (in which case we update it, and
    // reset the flag that indices that this needs to happen to zero), and
    // then return it. Make this thread-safe by using a mutex to guard
    // all of this:
    std::lock_guard<std::mutex> lock(cell_bounding_boxes_rtree_mutex);

    if (update_flags & update_cell_bounding_boxes_rtree)
      {
        std::vector<std::pair<
          BoundingBox<spacedim>,
          typename Triangulation<dim, spacedim>::active_cell_iterator>>
          boxes;
        boxes.reserve(tria->n_active_cells());
        for (const auto &cell : tria->active_cell_iterators())
          boxes.emplace_back(mapping->get_bounding_box(cell), cell);

        cell_bounding_boxes_rtree = pack_rtree(boxes);

        // Atomically clear the flag that indicates that this data member
        // needs to be updated:
        update_flags &= ~update_cell_bounding_boxes_rtree;
      }
    return cell_bounding_boxes_rtree;
  }



  template <int dim, int spacedim>
  const RTree<
    std::pair<BoundingBox<spacedim>,
              typename Triangulation<dim, spacedim>::active_cell_iterator>> &
  Cache<dim, spacedim>::get_locally_owned_cell_bounding_boxes_rtree() const
  {
    // In the following, we will first check whether the data structure
    // in question needs to be updated (in which case we update it, and
    // reset the flag that indices that this needs to happen to zero), and
    // then return it. Make this thread-safe by using a mutex to guard
    // all of this:
    std::lock_guard<std::mutex> lock(
      locally_owned_cell_bounding_boxes_rtree_mutex);

    if (update_flags & update_locally_owned_cell_bounding_boxes_rtree)
      {
        std::vector<std::pair<
          BoundingBox<spacedim>,
          typename Triangulation<dim, spacedim>::active_cell_iterator>>
          boxes;
        if (const parallel::TriangulationBase<dim, spacedim> *parallel_tria =
              dynamic_cast<const parallel::TriangulationBase<dim, spacedim> *>(
                &*tria))
          boxes.reserve(parallel_tria->n_locally_owned_active_cells());
        else
          boxes.reserve(tria->n_active_cells());
        for (const auto &cell : tria->active_cell_iterators() |
                                  IteratorFilters::LocallyOwnedCell())
          boxes.emplace_back(mapping->get_bounding_box(cell), cell);

        locally_owned_cell_bounding_boxes_rtree = pack_rtree(boxes);

        // Atomically clear the flag that indicates that this data member
        // needs to be updated:
        update_flags &= ~update_locally_owned_cell_bounding_boxes_rtree;
      }
    return locally_owned_cell_bounding_boxes_rtree;
  }



  template <int dim, int spacedim>
  const RTree<std::pair<BoundingBox<spacedim>, unsigned int>> &
  Cache<dim, spacedim>::get_covering_rtree(const unsigned int level) const
  {
    // In the following, we will first check whether the data structure
    // in question needs to be updated (in which case we update it, and
    // reset the flag that indices that this needs to happen to zero), and
    // then return it. Make this thread-safe by using a mutex to guard
    // all of this:
    std::lock_guard<std::mutex> lock(covering_rtree_mutex);

    if (update_flags & update_covering_rtree ||
        covering_rtree.find(level) == covering_rtree.end())
      {
        const auto boxes =
          extract_rtree_level(get_locally_owned_cell_bounding_boxes_rtree(),
                              level);

        if (const auto tria_mpi =
              dynamic_cast<const parallel::TriangulationBase<dim, spacedim> *>(
                &(*tria)))
          {
            covering_rtree[level] = GridTools::build_global_description_tree(
              boxes, tria_mpi->get_mpi_communicator());
          }
        else
          {
            covering_rtree[level] =
              GridTools::build_global_description_tree(boxes, MPI_COMM_SELF);
          }

        // Atomically clear the flag that indicates that this data member
        // needs to be updated:
        update_flags &= ~update_covering_rtree;
      }

    return covering_rtree[level];
  }



  template <int dim, int spacedim>
  const std::vector<std::set<unsigned int>> &
  Cache<dim, spacedim>::get_vertex_to_neighbor_subdomain() const
  {
    // In the following, we will first check whether the data structure
    // in question needs to be updated (in which case we update it, and
    // reset the flag that indices that this needs to happen to zero), and
    // then return it. Make this thread-safe by using a mutex to guard
    // all of this:
    std::lock_guard<std::mutex> lock(vertex_to_neighbor_subdomain_mutex);

    if (update_flags & update_vertex_to_neighbor_subdomain)
      {
        vertex_to_neighbor_subdomain.clear();
        vertex_to_neighbor_subdomain.resize(tria->n_vertices());
        for (const auto &cell : tria->active_cell_iterators())
          {
            if (cell->is_ghost())
              for (const unsigned int v : cell->vertex_indices())
                vertex_to_neighbor_subdomain[cell->vertex_index(v)].insert(
                  cell->subdomain_id());
          }
        // Atomically clear the flag that indicates that this data member
        // needs to be updated:
        update_flags &= ~update_vertex_to_neighbor_subdomain;
      }
    return vertex_to_neighbor_subdomain;
  }



  template <int dim, int spacedim>
  const std::map<unsigned int, std::set<types::subdomain_id>> &
  Cache<dim, spacedim>::get_vertices_with_ghost_neighbors() const
  {
    // In the following, we will first check whether the data structure
    // in question needs to be updated (in which case we update it, and
    // reset the flag that indices that this needs to happen to zero), and
    // then return it. Make this thread-safe by using a mutex to guard
    // all of this:
    std::lock_guard<std::mutex> lock(vertices_with_ghost_neighbors_mutex);

    if (update_flags & update_vertex_with_ghost_neighbors)
      {
        vertices_with_ghost_neighbors =
          GridTools::compute_vertices_with_ghost_neighbors(*tria);

        // Atomically clear the flag that indicates that this data member
        // needs to be updated:
        update_flags &= ~update_vertex_with_ghost_neighbors;
      }

    return vertices_with_ghost_neighbors;
  }

#include "grid/grid_tools_cache.inst"

} // namespace GridTools

DEAL_II_NAMESPACE_CLOSE
