// ---------------------------------------------------------------------
//
// Copyright (C) 2017 - 2020 by the deal.II authors
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
  const RTree<std::pair<BoundingBox<spacedim>, unsigned int>> &
  Cache<dim, spacedim>::get_covering_rtree() const
  {
    if (update_flags & update_covering_rtree)
      {
        // TO DO: the locally owned portion of the domain
        // might consists of more separate pieces: a single
        // bounding box might not always be the best choice.

        // Note: we create the local variable ptr here, because gcc 4.8.5
        // fails to compile if we pass the variable directly.
        IteratorFilters::LocallyOwnedCell locally_owned_cell_predicate;
        std::function<bool(
          const typename Triangulation<dim, spacedim>::active_cell_iterator &)>
          ptr(locally_owned_cell_predicate);

        const BoundingBox<spacedim> bbox =
          GridTools::compute_bounding_box(get_triangulation(), ptr);


        std::vector<BoundingBox<spacedim>> bbox_v(1, bbox);
        if (const auto tria_mpi =
              dynamic_cast<const parallel::TriangulationBase<dim, spacedim> *>(
                &(*tria)))
          {
            covering_rtree = GridTools::build_global_description_tree(
              bbox_v, tria_mpi->get_communicator());
          }
        else
          {
            covering_rtree =
              GridTools::build_global_description_tree(bbox_v, MPI_COMM_SELF);
          }
        update_flags = update_flags & ~update_covering_rtree;
      }

    return covering_rtree;
  }

#include "grid_tools_cache.inst"

} // namespace GridTools

DEAL_II_NAMESPACE_CLOSE
