// ---------------------------------------------------------------------
//
// Copyright (C) 2017 - 2022 by the deal.II authors
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

#include <deal.II/base/bounding_box.h>
#include <deal.II/base/mpi.h>

#include <deal.II/grid/filtered_iterator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/grid_tools_cache.h>

DEAL_II_NAMESPACE_OPEN

namespace GridTools
{
  template <int dim, int spacedim>
  Cache<dim, spacedim>::Cache(const Triangulation<dim, spacedim> &tria,
                              const Mapping<dim, spacedim> &      mapping)
    : update_flags(update_all)
    , tria(&tria)
    , mapping(&mapping)
  {
    tria_signal =
      tria.signals.any_change.connect([&]() { mark_for_update(update_all); });
  }

  template <int dim, int spacedim>
  Cache<dim, spacedim>::~Cache()
  {
    // Make sure that the signals that was attached to the triangulation
    // is removed here.
    if (tria_signal.connected())
      tria_signal.disconnect();
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
    if (update_flags & update_vertex_to_cell_map)
      {
        vertex_to_cells = GridTools::vertex_to_cell_map(*tria);
        update_flags    = update_flags & ~update_vertex_to_cell_map;
      }
    return vertex_to_cells;
  }



  template <int dim, int spacedim>
  const std::vector<std::vector<Tensor<1, spacedim>>> &
  Cache<dim, spacedim>::get_vertex_to_cell_centers_directions() const
  {
    if (update_flags & update_vertex_to_cell_centers_directions)
      {
        vertex_to_cell_centers = GridTools::vertex_to_cell_centers_directions(
          *tria, get_vertex_to_cell_map());
        update_flags = update_flags & ~update_vertex_to_cell_centers_directions;
      }
    return vertex_to_cell_centers;
  }



  template <int dim, int spacedim>
  const std::map<unsigned int, Point<spacedim>> &
  Cache<dim, spacedim>::get_used_vertices() const
  {
    if (update_flags & update_used_vertices)
      {
        used_vertices = GridTools::extract_used_vertices(*tria, *mapping);
        update_flags  = update_flags & ~update_used_vertices;
      }
    return used_vertices;
  }



  template <int dim, int spacedim>
  const RTree<std::pair<Point<spacedim>, unsigned int>> &
  Cache<dim, spacedim>::get_used_vertices_rtree() const
  {
    if (update_flags & update_used_vertices_rtree)
      {
        const auto &used_vertices = get_used_vertices();
        std::vector<std::pair<Point<spacedim>, unsigned int>> vertices(
          used_vertices.size());
        unsigned int i = 0;
        for (const auto &it : used_vertices)
          vertices[i++] = std::make_pair(it.second, it.first);
        used_vertices_rtree = pack_rtree(vertices);
        update_flags        = update_flags & ~update_used_vertices_rtree;
      }
    return used_vertices_rtree;
  }



  template <int dim, int spacedim>
  const RTree<
    std::pair<BoundingBox<spacedim>,
              typename Triangulation<dim, spacedim>::active_cell_iterator>> &
  Cache<dim, spacedim>::get_cell_bounding_boxes_rtree() const
  {
    if (update_flags & update_cell_bounding_boxes_rtree)
      {
        std::vector<std::pair<
          BoundingBox<spacedim>,
          typename Triangulation<dim, spacedim>::active_cell_iterator>>
                     boxes(tria->n_active_cells());
        unsigned int i = 0;
        for (const auto &cell : tria->active_cell_iterators())
          boxes[i++] = std::make_pair(mapping->get_bounding_box(cell), cell);

        cell_bounding_boxes_rtree = pack_rtree(boxes);
        update_flags = update_flags & ~update_cell_bounding_boxes_rtree;
      }
    return cell_bounding_boxes_rtree;
  }



  template <int dim, int spacedim>
  const RTree<
    std::pair<BoundingBox<spacedim>,
              typename Triangulation<dim, spacedim>::active_cell_iterator>> &
  Cache<dim, spacedim>::get_locally_owned_cell_bounding_boxes_rtree() const
  {
    if (update_flags & update_locally_owned_cell_bounding_boxes_rtree)
      {
        std::vector<std::pair<
          BoundingBox<spacedim>,
          typename Triangulation<dim, spacedim>::active_cell_iterator>>
          boxes;
        boxes.reserve(tria->n_active_cells());
        for (const auto &cell : tria->active_cell_iterators() |
                                  IteratorFilters::LocallyOwnedCell())
          boxes.emplace_back(mapping->get_bounding_box(cell), cell);

        locally_owned_cell_bounding_boxes_rtree = pack_rtree(boxes);
        update_flags =
          update_flags & ~update_locally_owned_cell_bounding_boxes_rtree;
      }
    return locally_owned_cell_bounding_boxes_rtree;
  }



  template <int dim, int spacedim>
  const RTree<std::pair<BoundingBox<spacedim>, unsigned int>> &
  Cache<dim, spacedim>::get_covering_rtree(const unsigned int level) const
  {
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
              boxes, tria_mpi->get_communicator());
          }
        else
          {
            covering_rtree[level] =
              GridTools::build_global_description_tree(boxes, MPI_COMM_SELF);
          }
        update_flags = update_flags & ~update_covering_rtree;
      }

    return covering_rtree[level];
  }

  template <int dim, int spacedim>
  const std::vector<std::set<unsigned int>> &
  Cache<dim, spacedim>::get_vertex_to_neighbor_subdomain() const
  {
    if (update_flags & update_vertex_to_neighbor_subdomain)
      {
        vertex_to_neighbor_subdomain.clear();
        vertex_to_neighbor_subdomain.resize(tria->n_vertices());
        for (const auto &cell : tria->active_cell_iterators())
          {
            if (cell->is_ghost())
              for (const unsigned int v : GeometryInfo<dim>::vertex_indices())
                vertex_to_neighbor_subdomain[cell->vertex_index(v)].insert(
                  cell->subdomain_id());
          }
        update_flags = update_flags & ~update_vertex_to_neighbor_subdomain;
      }
    return vertex_to_neighbor_subdomain;
  }

#include "grid_tools_cache.inst"

} // namespace GridTools

DEAL_II_NAMESPACE_CLOSE
