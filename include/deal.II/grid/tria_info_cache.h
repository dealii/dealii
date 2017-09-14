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

#ifndef dealii_grid_tria_info_cache_h
#define dealii_grid_tria_info_cache_h


#include <deal.II/base/config.h>
#include <deal.II/base/subscriptor.h>
#include <deal.II/base/exceptions.h>
#include <deal.II/base/point.h>

#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_info_cache_flags.h>
#include <deal.II/fe/mapping_q1.h>

#include <cmath>

DEAL_II_NAMESPACE_OPEN

/**
 * A class that caches computationally intensive information about a
 * Triangulation.
 *
 * This class attaches a signal to the Triangulation passed at construction
 * time, and rebuilds all the structures specified in the list of specified
 * TriangulationInfoCacheFlags, whenever the Triangulation is changed.
 *
 * If you try to access one of the objects that has not been cached (due to a
 * missing TriangulationInfoCacheFlags flag), an Assertion is thrown.
 *
 * Notice that this class only notices if the underlying Triangulation has
 * changed due to a Triangulation::Signals::any_change() signal being triggered.
 *
 * If the triangulation changes due to a MappingQEulerian being passed to this
 * class, or because you manually change some vertex locations, then some of
 * the structures in this class become obsolete, and you will have to call the
 * method update() manually.
 *
 * @author Luca Heltai, 2017.
 */
template <int dim, int spacedim>
class TriangulationInfoCache : public Subscriptor
{
public:
  /**
   * Constructor.
   *
   * You can specify what to cache by passing appropriate
   * TriangulationInfoCacheFlags. If you provide the optional `mapping`
   * argument, then this is used whenever a mapping is required.
   *
   * @param tria The triangulation for which to store information
   * @param flags TriangulationInfoCacheFlags that specify what to cache
   * @param mapping The mapping to use when computing cached objects
   */
  TriangulationInfoCache (const Triangulation<dim,spacedim> &tria,
                          const TriangulationInfoCacheFlags &flags = cache_nothing,
                          const Mapping<dim,spacedim> &mapping=StaticMappingQ1<dim,spacedim>::mapping);

  /**
   * Loop over all TriangulationInfoCacheFlags and rebuilds the specific data
   * structure.
   *
   * The optional parameter can be set to true to update only data that is
   * dependent on the position of the vertices. This may happen, for example,
   * if you have a MappingQEulerian class associated to the triangulation. In
   * this case, even if the Triangulation has not changed, the vector that
   * specifies the position of the vertices may actually change *independently*
   * of the triangulation, requiring you to recompute vertex dependent things.
   * In this case, the topological information is not changed, and does not
   * need to be updated.
   *
   * @param topology_is_unchanged
   */
  void update(bool topology_is_unchanged=false);


  /**
   * Return the cached vertex_to_cell_map as computed by GridTools::vertex_to_cell_map().
   */
  const std::vector< std::set<typename Triangulation<dim,spacedim>::active_cell_iterator> >
  &get_vertex_to_cell_map() const;

  /**
   * Return the cached vertex_to_cell_centers_directions as computed by
   * GridTools::vertex_to_cell_centers_directions().
   */
  const std::vector<std::vector<Tensor<1,spacedim>>>
  &get_vertex_to_cell_centers_directions() const;

  /**
   * Exception for uninitialized cache fields.
   */
  DeclException2(ExcAccessToInvalidCacheObject,
                 TriangulationInfoCacheFlags &,
                 TriangulationInfoCacheFlags &,
                 << "You tried to access a cache object that has not "
                 << "been constructed by this class. You specified "
                 << arg2 << " at construction time, but you need also "
                 << arg1 << " to access the field you are requesting now.");
private:
  /**
   * A pointer to the Triangulation.
   */
  SmartPointer<Triangulation<dim,spacedim>,
               TriangulationInfoCache<dim,spacedim>> tria;

  /**
   * Specifies what to cache.
   */
  const TriangulationInfoCacheFlags flags;

  /**
   * Mapping to use when computing on the Triangulation.
   */
  SmartPointer<Mapping<dim,spacedim>,
               TriangulationInfoCache<dim,spacedim>> mapping;


  /**
   * Store vertex to cell map information, as generated by
   * GridTools::vertex_to_cell_map()
   */
  std::vector< std::set<typename Triangulation<dim,spacedim>::active_cell_iterator> > vertex_to_cells;

  /**
   * Store vertex to cell center directions, as generated by
   * GridTools::vertex_to_cell_centers_directions.
   */
  std::vector<std::vector<Tensor<1,spacedim>>> vertex_to_cell_centers;
};

DEAL_II_NAMESPACE_CLOSE

#endif
