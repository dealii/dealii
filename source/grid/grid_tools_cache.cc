// ---------------------------------------------------------------------
//
// Copyright (C) 2017 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------

#include <deal.II/grid/grid_tools_cache.h>
#include <deal.II/grid/grid_tools.h>

DEAL_II_NAMESPACE_OPEN

namespace GridTools
{

  template<int dim, int spacedim>
  Cache<dim,spacedim>::Cache(
    const Triangulation<dim, spacedim> &tria,
    const TriangulationInfoCacheFlags &flags,
    const Mapping<dim, spacedim> &mapping) :
    tria(&tria),
    flags(flags),
    mapping(&mapping)
  {
    tria_signal = tria.signals.any_change.connect([&]()
    {
      update();
    });

    if (tria.n_active_cells()>0)
      update();
  }

  template<int dim, int spacedim>
  Cache<dim,spacedim>::~Cache()
  {
    // Make sure that the signals that was attached to the triangulation
    // is removed here.
    if (tria_signal.connected())
      tria_signal.disconnect();
  }



  template<int dim, int spacedim>
  void Cache<dim,spacedim>::update(bool topology_is_unchanged)
  {
    // If the triangulation is empty, just clear everything.
    if (tria->n_active_cells() == 0)
      {
        vertex_to_cells.clear();
        vertex_to_cell_centers.clear();
#ifdef DEAL_II_WITH_NANOFLANN
        vertex_kdtree.set_points(tria->get_vertices());
#endif
        return;
      }

    if (topology_is_unchanged == false)
      {
        if (cache_vertex_to_cell_map & flags)
          vertex_to_cells = GridTools::vertex_to_cell_map(*tria);

#ifdef DEAL_II_WITH_NANOFLANN
        if (cache_vertex_kdtree & flags)
          vertex_kdtree.set_points(tria->get_vertices());
#endif
      }

    if (cache_vertex_to_cell_centers_directions & flags)
      vertex_to_cell_centers =
        GridTools::vertex_to_cell_centers_directions(*tria, vertex_to_cells);

  }



  template<int dim, int spacedim>
  const std::vector<std::set<typename Triangulation<dim,spacedim>::active_cell_iterator> > &
  Cache<dim,spacedim>::get_vertex_to_cell_map() const
  {
    Assert(flags & cache_vertex_to_cell_map,
           ExcAccessToInvalidCacheObject(cache_vertex_to_cell_map, flags));
    return vertex_to_cells;
  }



  template<int dim, int spacedim>
  const std::vector<std::vector<Tensor<1,spacedim>>> &
  Cache<dim,spacedim>::get_vertex_to_cell_centers_directions() const
  {
    Assert(flags & cache_vertex_to_cell_centers_directions,
           ExcAccessToInvalidCacheObject(cache_vertex_to_cell_centers_directions, flags));
    return vertex_to_cell_centers;
  }



#ifdef DEAL_II_WITH_NANOFLANN
  template<int dim, int spacedim>
  const KDTree<spacedim> &Cache<dim,spacedim>::get_vertex_kdtree() const
  {
    Assert(flags & cache_vertex_kdtree,
           ExcAccessToInvalidCacheObject(cache_vertex_kdtree, flags));
    return vertex_kdtree;
  }
#endif

#include "grid_tools_cache.inst"

} // GridTools

DEAL_II_NAMESPACE_CLOSE

