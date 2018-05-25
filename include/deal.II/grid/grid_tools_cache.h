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

#ifndef dealii_grid_grid_tools_cache_h
#define dealii_grid_grid_tools_cache_h


#include <deal.II/base/config.h>

#include <deal.II/base/exceptions.h>
#include <deal.II/base/point.h>
#include <deal.II/base/subscriptor.h>

#include <deal.II/fe/mapping_q1.h>

#include <deal.II/grid/grid_tools_cache_update_flags.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/numerics/kdtree.h>

#include <boost/signals2.hpp>

#include <cmath>

DEAL_II_NAMESPACE_OPEN

namespace GridTools
{
  /**
   * A class that caches computationally intensive information about a
   * Triangulation.
   *
   * This class attaches a signal to the Triangulation passed at construction
   * time to keep track about refinement changes, and allows the user to query
   * some of the data structures constructed using functions in the GridTools
   * namespace which are computed only once, and then cached inside this class
   * for faster access whenever the triangulation has not changed.
   *
   * Notice that this class only notices if the underlying Triangulation has
   * changed due to a Triangulation::Signals::any_change() signal being
   * triggered.
   *
   * If the triangulation changes for other reasons, for example because you
   * use it in conjunction with a MappingQEulerian object that sees the
   * vertices through its own transformation, or because you manually change
   * some vertex locations, then some of the structures in this class become
   * obsolete, and you will have to mark them as outdated, by calling the
   * method mark_for_update() manually.
   *
   * @author Luca Heltai, 2017.
   */
  template <int dim, int spacedim = dim>
  class Cache : public Subscriptor
  {
  public:
    /**
     * Constructor.
     *
     * If you provide the optional `mapping` argument, then this is used
     * whenever a mapping is required.
     *
     * @param tria The triangulation for which to store information
     * @param mapping The mapping to use when computing cached objects
     */
    Cache(const Triangulation<dim, spacedim> &tria,
          const Mapping<dim, spacedim> &      mapping =
            StaticMappingQ1<dim, spacedim>::mapping);

    /**
     * Destructor.
     */
    ~Cache() override;

    /**
     * Make sure that the objects marked for update are recomputed during
     * subsequent calls to the `get_*` functions defined in this class.
     *
     * Notice that no work is performed when you call this function. The actual
     * data structures are computed the next time you call the corresponding
     * `get_*` method.
     *
     * @param flags What to mark for update
     */
    void
    mark_for_update(const CacheUpdateFlags &flags = update_all);


    /**
     * Return the cached vertex_to_cell_map as computed by
     * GridTools::vertex_to_cell_map().
     */
    const std::vector<
      std::set<typename Triangulation<dim, spacedim>::active_cell_iterator>> &
    get_vertex_to_cell_map() const;

    /**
     * Return the cached vertex_to_cell_centers_directions as computed by
     * GridTools::vertex_to_cell_centers_directions().
     */
    const std::vector<std::vector<Tensor<1, spacedim>>> &
    get_vertex_to_cell_centers_directions() const;

    /**
     * Return the cached map of used vertices, as computed by
     * GridTools::extract_used_vertices().
     */
    const std::map<unsigned int, Point<spacedim>> &
    get_used_vertices() const;

    /**
     * Return a reference to the stored triangulation.
     */
    const Triangulation<dim, spacedim> &
    get_triangulation() const;

    /**
     * Return a reference to the stored mapping.
     */
    const Mapping<dim, spacedim> &
    get_mapping() const;

#ifdef DEAL_II_WITH_NANOFLANN
    /**
     * Return the cached vertex_kdtree object, constructed with the vertices of
     * the stored triangulation.
     */
    const KDTree<spacedim> &
    get_vertex_kdtree() const;
#endif

  private:
    /**
     * Keep track of what needs to be updated next.
     */
    mutable CacheUpdateFlags update_flags;

    /**
     * A pointer to the Triangulation.
     */
    SmartPointer<const Triangulation<dim, spacedim>, Cache<dim, spacedim>> tria;

    /**
     * Mapping to use when computing on the Triangulation.
     */
    SmartPointer<const Mapping<dim, spacedim>, Cache<dim, spacedim>> mapping;


    /**
     * Store vertex to cell map information, as generated by
     * GridTools::vertex_to_cell_map()
     */
    mutable std::vector<
      std::set<typename Triangulation<dim, spacedim>::active_cell_iterator>>
      vertex_to_cells;

    /**
     * Store vertex to cell center directions, as generated by
     * GridTools::vertex_to_cell_centers_directions().
     */
    mutable std::vector<std::vector<Tensor<1, spacedim>>>
      vertex_to_cell_centers;

#ifdef DEAL_II_WITH_NANOFLANN
    /**
     * A KDTree object constructed with the vertices of the triangulation.
     */
    mutable KDTree<spacedim> vertex_kdtree;
#endif

    /**
     * Store the used vertices of the Triangulation, as generated by
     * GridTools::extract_used_vertices().
     */
    mutable std::map<unsigned int, Point<spacedim>> used_vertices;

    /**
     * Storage for the status of the triangulation signal.
     */
    boost::signals2::connection tria_signal;
  };



  // Inline functions
  template <int dim, int spacedim>
  inline const Triangulation<dim, spacedim> &
  Cache<dim, spacedim>::get_triangulation() const
  {
    return *tria;
  }



  template <int dim, int spacedim>
  inline const Mapping<dim, spacedim> &
  Cache<dim, spacedim>::get_mapping() const
  {
    return *mapping;
  }
} // namespace GridTools



DEAL_II_NAMESPACE_CLOSE

#endif
