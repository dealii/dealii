// ---------------------------------------------------------------------
//
// Copyright (C) 2017 - 2019 by the deal.II authors
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

#include <deal.II/numerics/rtree.h>

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
            (ReferenceCells::get_hypercube<dim>()
               .template get_default_linear_mapping<dim, spacedim>()));

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
     * Return the cached RTree object for the vertices, constructed using the
     * used vertices of the triangulation.
     */
    const RTree<std::pair<Point<spacedim>, unsigned int>> &
    get_used_vertices_rtree() const;

    /**
     * Return the cached RTree object of the cell bounding boxes, constructed
     * using the active cell iterators of the stored triangulation. For
     * parallel::distributed::Triangulation objects, this function will return
     * also the bounding boxes of ghost and artificial cells.
     */
    const RTree<
      std::pair<BoundingBox<spacedim>,
                typename Triangulation<dim, spacedim>::active_cell_iterator>> &
    get_cell_bounding_boxes_rtree() const;

    /**
     * Return the cached RTree object of bounding boxes containing locally owned
     * active cells, constructed using the active cell iterators of the stored
     * triangulation.
     *
     * In contrast to the previous function, this function builds the RTree
     * using only the locally owned cells, i.e., not including ghost or
     * artificial cells. The two functions return the same result in serial
     * computations.
     */
    const RTree<
      std::pair<BoundingBox<spacedim>,
                typename Triangulation<dim, spacedim>::active_cell_iterator>> &
    get_locally_owned_cell_bounding_boxes_rtree() const;


    /**
     * Returns the vector of set of integer containing the subdomain id
     * to which each vertex is connected to. This feature is used extensively
     * in the particle_handler to detect on which processors ghost particles
     * must be built.
     */
    const std::vector<std::set<unsigned int>> &
    get_vertex_to_neighbor_subdomain() const;

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


    /**
     * This function returns an object that allows identifying
     * which process(es) in a parallel computation may own the
     * cell that surrounds a given point. The elements of this
     * object -- an Rtree -- are pairs of bounding boxes denoting
     * areas that cover all or parts of the local portion of a
     * parallel triangulation, and an unsigned int representing
     * the process or subdomain that owns the cells falling within the given
     * bounding box.
     *
     * Given a point on a parallel::TriangulationBase, this tree
     * allows to identify one, or few candidate processes, for
     * which the point lies on a locally owned cell.
     *
     * Constructing or updating the rtree requires a call to
     * GridTools::build_global_description_tree(), which exchanges
     * bounding boxes between all processes using
     * Utilities::MPI::all_gather(), a collective operation.
     * Therefore this function must be called by all processes
     * at the same time.
     *
     * The local bounding boxes are constructed by extracting the
     * specified @p level from the rtree object returned by the
     * get_locally_owned_cell_bounding_boxes_rtree(). Notice that @p level
     * in this context does not refer to the level of the triangulation, but
     * refer to the level of the RTree object (see, e.g.,
     * https://en.wikipedia.org/wiki/R-tree).
     */
    const RTree<std::pair<BoundingBox<spacedim>, unsigned int>> &
    get_covering_rtree(const unsigned int level = 0) const;

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

    /**
     * A collection of rtree objects covering the whole mesh.
     *
     * Each entry of the map is constructed from the function
     * extract_rtree_level() applied to the rtree returned by the function
     * get_locally_owned_cell_bounding_boxes_rtree(), with the specified level
     * in input.
     */
    mutable std::map<unsigned int,
                     RTree<std::pair<BoundingBox<spacedim>, unsigned int>>>
      covering_rtree;

    /**
     * Store the used vertices of the Triangulation, as generated by
     * GridTools::extract_used_vertices().
     */
    mutable std::map<unsigned int, Point<spacedim>> used_vertices;

    /**
     * Store an RTree object, containing the used vertices of the triangulation.
     */
    mutable RTree<std::pair<Point<spacedim>, unsigned int>> used_vertices_rtree;

    /**
     * Store an RTree object, containing the bounding boxes of the cells of the
     * triangulation.
     */
    mutable RTree<
      std::pair<BoundingBox<spacedim>,
                typename Triangulation<dim, spacedim>::active_cell_iterator>>
      cell_bounding_boxes_rtree;

    /**
     * Store an RTree object, containing the bounding boxes of the locally owned
     * cells of the triangulation.
     */
    mutable RTree<
      std::pair<BoundingBox<spacedim>,
                typename Triangulation<dim, spacedim>::active_cell_iterator>>
      locally_owned_cell_bounding_boxes_rtree;


    /**
     * Store an std::vector of std::set of integer containing the id of all
     * subdomain to which a vertex is connected to.
     */
    mutable std::vector<std::set<unsigned int>> vertex_to_neighbor_subdomain;

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
