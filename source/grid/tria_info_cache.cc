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

#include <deal.II/grid/tria_info_cache.h>
#include <deal.II/grid/grid_tools.h>

DEAL_II_NAMESPACE_OPEN

template<int dim, int spacedim>
TriangulationInfoCache<dim,spacedim>::TriangulationInfoCache(
  const Triangulation<dim, spacedim> &tria,
  const TriangulationInfoCacheFlags &flags,
  const Mapping<dim, spacedim> &mapping) :
  tria(&tria),
  flags(flags),
  mapping(&mapping)
{
  tria->signals.any_change.connect([&]()
  {
    update();
  });

  if (tria.n_active_cells()>0)
    update();
}



template<int dim, int spacedim>
void TriangulationInfoCache<dim,spacedim>::update(bool topology_is_unchanged)
{
  if (topology_is_unchanged == false)
    {
      if (cache_vertex_to_cell_map & flags)
        vertex_to_cells = GridTools::vertex_to_cell_map(*tria);
    }

  if (cache_vertex_to_cell_centers_directions & flags)
    vertex_to_cell_centers =
      GridTools::vertex_to_cell_centers_directions(*tria, vertex_to_cells);
}



template<int dim, int spacedim>
const std::vector<std::set<typename Triangulation<dim,spacedim>::active_cell_iterator> > &
TriangulationInfoCache<dim,spacedim>::get_vertex_to_cell_map() const
{
  Assert(flags & cache_vertex_to_cell_map,
         ExcAccessToInvalidCacheObject(cache_vertex_to_cell_map, flags));
  return vertex_to_cells;
}



template<int dim, int spacedim>
const std::vector<std::vector<Tensor<1,spacedim>>> &
TriangulationInfoCache<dim,spacedim>::get_vertex_to_cell_centers_directions() const
{
  Assert(flags & cache_vertex_to_cell_centers_directions,
         ExcAccessToInvalidCacheObject(cache_vertex_to_cell_centers_directions, flags));
  return vertex_to_cell_centers;
}

DEAL_II_NAMESPACE_CLOSE
