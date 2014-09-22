// ---------------------------------------------------------------------
//
// Copyright (C) 2001 - 2014 by the deal.II authors
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

#ifndef __deal2__grid_tools_H
#define __deal2__grid_tools_H


#include <deal.II/base/config.h>
#include <deal.II/fe/mapping.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/fe/mapping_q1.h>

#include <bitset>
#include <list>

DEAL_II_NAMESPACE_OPEN


template <int, int> class DoFHandler;
template <int, int> class Mapping;
namespace hp
{
  template <int, int> class DoFHandler;
  template <int, int> class MappingCollection;
}

class SparsityPattern;


/**
 * This namespace is a collection of algorithms working on triangulations,
 * such as shifting or rotating triangulations, but also finding a
 * cell that contains a given point. See the descriptions of the
 * individual functions for more information.
 *
 * @ingroup grid
 */
namespace GridTools
{
  /**
   *  @name Information about meshes and cells
   */
  /*@{*/

  /**
   * Return the diameter of a
   * triangulation. The diameter is
   * computed using only the
   * vertices, i.e. if the diameter
   * should be larger than the
   * maximal distance between
   * boundary vertices due to a
   * higher order mapping, then
   * this function will not catch
   * this.
   */
  template <int dim, int spacedim>
  double diameter (const Triangulation<dim, spacedim> &tria);

  /**
   * Compute the volume (i.e. the dim-dimensional measure) of the
   * triangulation. We compute the measure using the integral
   * $\int 1 \; dx$. The integral approximated is approximated
   * via quadrature for which we need the mapping argument.
   *
   * This function also works for objects of type
   * parallel::distributed::Triangulation, in which case the
   * function is a collective operation.
   */
  template <int dim, int spacedim>
  double volume (const Triangulation<dim,spacedim> &tria,
                 const Mapping<dim,spacedim> &mapping = (StaticMappingQ1<dim,spacedim>::mapping));

  /**
   * Return the diamater of the smallest
   * active cell of a triangulation. See
   * step-24 for an example
   * of use of this function.
   */
  template <int dim, int spacedim>
  double
  minimal_cell_diameter (const Triangulation<dim, spacedim> &triangulation);

  /**
   * Return the diamater of the largest
   * active cell of a triangulation.
   */
  template <int dim, int spacedim>
  double
  maximal_cell_diameter (const Triangulation<dim, spacedim> &triangulation);

  /**
   * Given a list of vertices (typically
   * obtained using
   * Triangulation::get_vertices) as the
   * first, and a list of vertex indices
   * that characterize a single cell as the
   * second argument, return the measure
   * (area, volume) of this cell. If this
   * is a real cell, then you can get the
   * same result using
   * <code>cell-@>measure()</code>, but
   * this function also works for cells
   * that do not exist except that you make
   * it up by naming its vertices from the
   * list.
   */
  template <int dim>
  double cell_measure (const std::vector<Point<dim> > &all_vertices,
                       const unsigned int (&vertex_indices)[GeometryInfo<dim>::vertices_per_cell]);

  /*@}*/
  /**
   *  @name Functions supporting the creation of meshes
   */
  /*@{*/

  /**
   * Remove vertices that are not
   * referenced by any of the
   * cells. This function is called
   * by all <tt>GridIn::read_*</tt>
   * functions to eliminate
   * vertices that are listed in
   * the input files but are not
   * used by the cells in the input
   * file. While these vertices
   * should not be in the input
   * from the beginning, they
   * sometimes are, most often when
   * some cells have been removed
   * by hand without wanting to
   * update the vertex lists, as
   * they might be lengthy.
   *
   * This function is called by all
   * <tt>GridIn::read_*</tt>
   * functions as the triangulation
   * class requires them to be
   * called with used vertices
   * only. This is so, since the
   * vertices are copied verbatim
   * by that class, so we have to
   * eliminate unused vertices
   * beforehand.
   *
   * Not implemented for the
   * codimension one case.
   */
  template <int dim, int spacedim>
  void delete_unused_vertices (std::vector<Point<spacedim> >    &vertices,
                               std::vector<CellData<dim> > &cells,
                               SubCellData                 &subcelldata);

  /**
   * Remove vertices that are duplicated,
   * due to the input of a structured grid,
   * for example. If these vertices are not
   * removed, the faces bounded by these
   * vertices become part of the boundary,
   * even if they are in the interior of
   * the mesh.
   *
   * This function is called by some
   * <tt>GridIn::read_*</tt> functions. Only
   * the vertices with indices in @p
   * considered_vertices are tested for
   * equality. This speeds up the algorithm,
   * which is quadratic and thus quite slow
   * to begin with. However, if you wish to
   * consider all vertices, simply pass an
   * empty vector.
   *
   * Two vertices are considered equal if
   * their difference in each coordinate
   * direction is less than @p tol.
   */
  template <int dim, int spacedim>
  void delete_duplicated_vertices (std::vector<Point<spacedim> >    &all_vertices,
                                   std::vector<CellData<dim> > &cells,
                                   SubCellData                 &subcelldata,
                                   std::vector<unsigned int>   &considered_vertices,
                                   const double                 tol=1e-12);

  /*@}*/
  /**
   *  @name Rotating, stretching and otherwise transforming meshes
   */
  /*@{*/

  /**
   * Transform the vertices of the given
   * triangulation by applying the
   * function object provided as first argument to all its vertices.
   *
   * The transformation given as
   * argument is used to transform
   * each vertex. Its respective
   * type has to offer a
   * function-like syntax, i.e. the
   * predicate is either an object
   * of a type that has an
   * <tt>operator()</tt>, or it is a
   * pointer to the function. In
   * either case, argument and
   * return value have to be of
   * type <tt>Point@<spacedim@></tt>.
   *
   * @note If you are using a parallel::distributed::Triangulation you will have
   * hanging nodes in your local Triangulation even if your "global" mesh has
   * no hanging nodes. This will cause issues with wrong positioning of hanging
   * nodes in ghost cells if you call the current functions: The vertices of
   * all locally owned cells will be correct,
   * but the vertices of some ghost cells may not. This means that
   * computations like KellyErrorEstimator may give wrong answers. A safe approach is
   * to use this function prior to any refinement in parallel, if that is possible, but
   * not after you refine the mesh.
   *
   * This function is used in the
   * "Possibilities for extensions" section
   * of step-38. It is also used in step-49 and step-53.
   */
  template <int dim, typename Transformation, int spacedim>
  void transform (const Transformation        &transformation,
                  Triangulation<dim,spacedim> &triangulation);

  /**
   * Shift each vertex of the
   * triangulation by the given
   * shift vector. This function
   * uses the transform()
   * function above, so the
   * requirements on the
   * triangulation stated there
   * hold for this function as
   * well.
   */
  template <int dim, int spacedim>
  void shift (const Point<spacedim>   &shift_vector,
              Triangulation<dim,spacedim> &triangulation);


  /**
   * Rotate all vertices of the
   * given two-dimensional
   * triangulation in
   * counter-clockwise sense around
   * the origin of the coordinate
   * system by the given angle
   * (given in radians, rather than
   * degrees). This function uses
   * the transform() function
   * above, so the requirements on
   * the triangulation stated there
   * hold for this function as
   * well.
   */
  void rotate (const double      angle,
               Triangulation<2> &triangulation);

  /**
   * Transform the given triangulation smoothly to a different domain where,
   * typically, each of the vertices at the boundary of the triangulation is mapped to
   * the corresponding points in the @p new_points map.
   *
   * The way this function works is that it solves a Laplace equation for each
   * of the dim components of a displacement field that maps the current
   * domain into one described by @p new_points . The @p new_points array
   * therefore represents the boundary values of this displacement field.
   * The function then evaluates this displacement field at each vertex in
   * the interior and uses it to place the mapped vertex where the
   * displacement field locates it. Because the solution of the Laplace
   * equation is smooth, this guarantees a smooth mapping from the old
   * domain to the new one.
   *
   * @param[in] new_points The locations where a subset of the existing vertices
   * are to be placed. Typically, this would be a map from the vertex indices
   * of all nodes on the boundary to their new locations, thus completely
   * specifying the geometry of the mapped domain. However, it may also include
   * interior points if necessary and it does not need to include all boundary
   * vertices (although you then lose control over the exact shape of the mapped
   * domain).
   *
   * @param[in,out] tria The Triangulation object. This object is changed in-place,
   * i.e., the previous locations of vertices are overwritten.
   *
   * @param[in] coefficient An optional coefficient for the Laplace problem.
   * Larger values make cells less prone to deformation (effectively increasing their stiffness).
   * The coefficient is evaluated in the coordinate system of the old,
   * undeformed configuration of the triangulation as input, i.e., before
   * the transformation is applied.
   * Should this function be provided, sensible results can only be expected if
   * all coefficients are positive.
   *
   * @note This function is not currently implemented for the 1d case.
   */
  template <int dim>
  void laplace_transform (const std::map<unsigned int,Point<dim> > &new_points,
                          Triangulation<dim> &tria,
                          const Function<dim,double> *coefficient = 0);

  /**
   * Scale the entire triangulation
   * by the given factor. To
   * preserve the orientation of
   * the triangulation, the factor
   * must be positive.
   *
   * This function uses the
   * transform() function
   * above, so the requirements on
   * the triangulation stated there
   * hold for this function as
   * well.
   */
  template <int dim, int spacedim>
  void scale (const double        scaling_factor,
              Triangulation<dim, spacedim> &triangulation);

  /**
   * Distort the given triangulation by randomly
   * moving around all the vertices
   * of the grid.  The direction of
   * movement of each vertex is random, while the
   * length of the shift vector has
   * a value of @p factor times
   * the minimal length of the
   * active edges adjacent to this
   * vertex. Note that @p factor
   * should obviously be well below
   * <tt>0.5</tt>.
   *
   * If @p keep_boundary is set to
   * @p true (which is the
   * default), then boundary
   * vertices are not moved.
   */
  template <int dim, int spacedim>
  void distort_random (const double factor,
                       Triangulation<dim, spacedim> &triangulation,
                       const bool   keep_boundary=true);

  /*@}*/
  /**
   *  @name Finding cells and vertices of a triangulation
   */
  /*@{*/

  /**
   * Find and return the number of
   * the used vertex in a given
   * mesh that is located closest
   * to a given point.
   *
   * @param container A variable of a type that satisfies the
   *   requirements of a mesh container (see @ref GlossMeshAsAContainer).
   * @param p The point for which we want to find the closest vertex.
   * @return The index of the closest vertex found.
   *
   * @author Ralf B. Schulz, 2006
   */
  template <int dim, template <int, int> class Container, int spacedim>
  unsigned int
  find_closest_vertex (const Container<dim, spacedim> &container,
                       const Point<spacedim>     &p);

  /**
   * Find and return a vector of
   * iterators to active cells that
   * surround a given vertex with index @p vertex_index.
   *
   * For locally refined grids, the
   * vertex itself might not be a vertex
   * of all adjacent cells that are returned.
   * However, it will
   * always be either a vertex of a cell or be
   * a hanging node located on a face or an
   * edge of it.
   *
   * @param container A variable of a type that satisfies the
   *   requirements of a mesh container (see @ref GlossMeshAsAContainer).
   * @param vertex_index The index of the vertex for which we try to
   *   find adjacent cells.
   * @return A vector of cells that lie adjacent to the given vertex.
   *
   * @note If the point requested does not lie in any of the cells of
   * the mesh given, then this function throws an exception of type
   * GridTools::ExcPointNotFound. You can catch this exception and
   * decide what to do in that case.
   *
   * @note It isn't entirely clear at this time whether the function
   * does the right thing with anisotropically refined meshes. It needs
   * to be checked for this case.
   */
  template<int dim, template <int, int> class Container, int spacedim>
  std::vector<typename Container<dim,spacedim>::active_cell_iterator>
  find_cells_adjacent_to_vertex (const Container<dim,spacedim> &container,
                                 const unsigned int    vertex_index);


  /**
   * Find and return an iterator to the active cell that surrounds a
   * given point.
   *
   * This is solely a wrapper function for the function of same name
   * below.  A Q1 mapping is used for the boundary, and the iterator
   * to the cell in which the point resides is returned.
   *
   * It is recommended to use the other version of this function, as
   * it simultaneously delivers the local coordinate of the given
   * point without additional computational cost.
   *
   * @param container A variable of a type that satisfies the
   *   requirements of a mesh container (see @ref GlossMeshAsAContainer).
   * @param p The point for which we want to find the surrounding cell.
   * @return An iterator into the mesh container that points to the
   *   surrounding cell.
   *
   * @note If the point requested does not lie in any of the cells of
   * the mesh given, then this function throws an exception of type
   * GridTools::ExcPointNotFound. You can catch this exception and
   * decide what to do in that case.
   *
   * @note When applied to a triangulation or DoF handler object based
   * on a parallel::distributed::Triangulation object, the cell
   * returned may in fact be a ghost or artificial cell (see
   * @ref GlossArtificialCell and @ref GlossGhostCell). If so,
   * many of the operations one may want to do on this cell
   * (e.g., evaluating the solution) may not be possible and you
   * will have to decide what to do in that case.
   */
  template <int dim, template <int,int> class Container, int spacedim>
  typename Container<dim,spacedim>::active_cell_iterator
  find_active_cell_around_point (const Container<dim,spacedim>  &container,
                                 const Point<spacedim> &p);

  /**
   * Find and return an iterator to the active cell that surrounds a
   * given point @p p.
   *
   * The algorithm used in this function proceeds by first looking for
   * vertex located closest to the given point, see
   * find_closest_vertex(). Secondly, all adjacent cells to this point
   * are found in the mesh, see find_cells_adjacent_to_vertex().
   * Lastly, for each of these cells, it is tested whether the point
   * is inside. This check is performed using arbitrary boundary
   * mappings.  Still, it is possible that due to roundoff errors, the
   * point cannot be located exactly inside the unit cell. In this
   * case, even points at a very small distance outside the unit cell
   * are allowed.
   *
   * If a point lies on the boundary of two or more cells, then the
   * algorithm tries to identify the cell that is of highest
   * refinement level.
   *
   * @param mapping The mapping used to determine whether the given
   *   point is inside a given cell.
   * @param container A variable of a type that satisfies the
   *   requirements of a mesh container (see @ref GlossMeshAsAContainer).
   * @param p The point for which we want to find the surrounding cell.
   * @return An pair of an iterator into the mesh container that points to the
   *   surrounding cell, and of the coordinates of that point inside the cell
   *   in the reference coordinates of that cell. This local
   *   position might be located slightly outside an actual unit cell,
   *   due to numerical roundoff.  Therefore, the point returned by this
   *   function should be projected onto the unit cell, using
   *   GeometryInfo::project_to_unit_cell().  This is not automatically
   *   performed by the algorithm.
   *
   * @note If the point requested does not lie in any of the cells of
   * the mesh given, then this function throws an exception of type
   * GridTools::ExcPointNotFound. You can catch this exception and
   * decide what to do in that case.
   *
   * @note When applied to a triangulation or DoF handler object based
   * on a parallel::distributed::Triangulation object, the cell
   * returned may in fact be a ghost or artificial cell (see
   * @ref GlossArtificialCell and @ref GlossGhostCell). If so,
   * many of the operations one may want to do on this cell
   * (e.g., evaluating the solution) may not be possible and you
   * will have to decide what to do in that case.
   */
  template <int dim, template<int, int> class Container, int spacedim>
  std::pair<typename Container<dim,spacedim>::active_cell_iterator, Point<dim> >
  find_active_cell_around_point (const Mapping<dim,spacedim>   &mapping,
                                 const Container<dim,spacedim> &container,
                                 const Point<spacedim>     &p);

  /**
   * A version of the previous function where we use that mapping on a
   * given cell that corresponds to the active finite element index of
   * that cell. This is obviously only useful for hp problems, since
   * the active finite element index for all other DoF handlers is
   * always zero.
   *
   * @note If the point requested does not lie in any of the cells of
   * the mesh given, then this function throws an exception of type
   * GridTools::ExcPointNotFound. You can catch this exception and
   * decide what to do in that case.
   *
   * @note When applied to a triangulation or DoF handler object based
   * on a parallel::distributed::Triangulation object, the cell
   * returned may in fact be a ghost or artificial cell (see
   * @ref GlossArtificialCell and @ref GlossGhostCell). If so,
   * many of the operations one may want to do on this cell
   * (e.g., evaluating the solution) may not be possible and you
   * will have to decide what to do in that case.
   */
  template <int dim, int spacedim>
  std::pair<typename hp::DoFHandler<dim,spacedim>::active_cell_iterator, Point<dim> >
  find_active_cell_around_point (const hp::MappingCollection<dim,spacedim>   &mapping,
                                 const hp::DoFHandler<dim,spacedim> &container,
                                 const Point<spacedim>     &p);

  /**
   * Return a list of all descendants of
   * the given cell that are active. For
   * example, if the current cell is once
   * refined but none of its children are
   * any further refined, then the returned
   * list will contain all its children.
   *
   * If the current cell is already active,
   * then the returned list is empty
   * (because the cell has no children that
   * may be active).
   *
   * @tparam Container A type that satisfies the
   *   requirements of a mesh container (see @ref GlossMeshAsAContainer).
   * @param cell An iterator pointing to a cell of the mesh container.
   * @return A list of active descendants of the given cell
   *
   * @note Since in C++ the type of the Container
   * template argument can
   * not be deduced from a function call,
   * you will have to specify it after the
   * function name, as for example in
   * @code
   *   GridTools::get_active_child_cells<DoFHandler<dim> > (cell)
   * @endcode
   */
  template <class Container>
  std::vector<typename Container::active_cell_iterator>
  get_active_child_cells (const typename Container::cell_iterator &cell);

  /**
   * Extract the active cells around a given
   * cell @p cell and return them in the
   * vector @p active_neighbors.
   *
   * @tparam Container A type that satisfies the
   *   requirements of a mesh container (see @ref GlossMeshAsAContainer).
   * @param cell[in] An iterator pointing to a cell of the mesh container.
   * @param active_neighbors[out] A list of active descendants of the given cell
   */
  template <class Container>
  void
  get_active_neighbors (const typename Container::active_cell_iterator        &cell,
                        std::vector<typename Container::active_cell_iterator> &active_neighbors);

  /*@}*/
  /**
   *  @name Partitions and subdomains of triangulations
   */
  /*@{*/

  /**
   * Produce a sparsity pattern in which
   * nonzero entries indicate that two
   * cells are connected via a common
   * face. The diagonal entries of the
   * sparsity pattern are also set.
   *
   * The rows and columns refer to the
   * cells as they are traversed in their
   * natural order using cell iterators.
   */
  template <int dim, int spacedim>
  void
  get_face_connectivity_of_cells (const Triangulation<dim, spacedim> &triangulation,
                                  SparsityPattern                    &connectivity);

  /**
   * Use the METIS partitioner to generate
   * a partitioning of the active cells
   * making up the entire domain. After
   * calling this function, the subdomain
   * ids of all active cells will have
   * values between zero and
   * @p n_partitions-1. You can access the
   * subdomain id of a cell by using
   * <tt>cell-@>subdomain_id()</tt>.
   *
   * This function will generate an error
   * if METIS is not installed unless
   * @p n_partitions is one. I.e., you can
   * write a program so that it runs in the
   * single-processor single-partition case
   * without METIS installed, and only
   * requires METIS when multiple
   * partitions are required.
   */
  template <int dim, int spacedim>
  void
  partition_triangulation (const unsigned int  n_partitions,
                           Triangulation<dim, spacedim> &triangulation);

  /**
   * This function does the same as the
   * previous one, i.e. it partitions a
   * triangulation using METIS into a
   * number of subdomains identified by the
   * <code>cell-@>subdomain_id()</code>
   * flag.
   *
   * The difference to the previous
   * function is the second argument, a
   * sparsity pattern that represents the
   * connectivity pattern between cells.
   *
   * While the function above builds it
   * directly from the triangulation by
   * considering which cells neighbor each
   * other, this function can take a more
   * refined connectivity graph. The
   * sparsity pattern needs to be of size
   * $N\times N$, where $N$ is the number
   * of active cells in the
   * triangulation. If the sparsity pattern
   * contains an entry at position $(i,j)$,
   * then this means that cells $i$ and $j$
   * (in the order in which they are
   * traversed by active cell iterators)
   * are to be considered connected; METIS
   * will then try to partition the domain
   * in such a way that (i) the subdomains
   * are of roughly equal size, and (ii) a
   * minimal number of connections are
   * broken.
   *
   * This function is mainly useful in
   * cases where connections between cells
   * exist that are not present in the
   * triangulation alone (otherwise the
   * previous function would be the simpler
   * one to use). Such connections may
   * include that certain parts of the
   * boundary of a domain are coupled
   * through symmetric boundary conditions
   * or integrals (e.g. friction contact
   * between the two sides of a crack in
   * the domain), or if a numerical scheme
   * is used that not only connects
   * immediate neighbors but a larger
   * neighborhood of cells (e.g. when
   * solving integral equations).
   *
   * In addition, this function may be
   * useful in cases where the default
   * sparsity pattern is not entirely
   * sufficient. This can happen because
   * the default is to just consider face
   * neighbors, not neighboring cells that
   * are connected by edges or
   * vertices. While the latter couple when
   * using continuous finite elements, they
   * are typically still closely connected
   * in the neighborship graph, and METIS
   * will not usually cut important
   * connections in this case. However, if
   * there are vertices in the mesh where
   * many cells (many more than the common
   * 4 or 6 in 2d and 3d, respectively)
   * come together, then there will be a
   * significant number of cells that are
   * connected across a vertex, but several
   * degrees removed in the connectivity
   * graph built only using face
   * neighbors. In a case like this, METIS
   * may sometimes make bad decisions and
   * you may want to build your own
   * connectivity graph.
   */
  template <int dim, int spacedim>
  void
  partition_triangulation (const unsigned int     n_partitions,
                           const SparsityPattern &cell_connection_graph,
                           Triangulation<dim,spacedim>    &triangulation);

  /**
   * For each active cell, return in the
   * output array to which subdomain (as
   * given by the <tt>cell->subdomain_id()</tt>
   * function) it belongs. The output array
   * is supposed to have the right size
   * already when calling this function.
   *
   * This function returns the association
   * of each cell with one subdomain. If
   * you are looking for the association of
   * each @em DoF with a subdomain, use the
   * <tt>DoFTools::get_subdomain_association</tt>
   * function.
   */
  template <int dim, int spacedim>
  void
  get_subdomain_association (const Triangulation<dim, spacedim>  &triangulation,
                             std::vector<types::subdomain_id> &subdomain);

  /**
   * Count how many cells are uniquely
   * associated with the given @p subdomain
   * index.
   *
   * This function may return zero
   * if there are no cells with the
   * given @p subdomain index. This
   * can happen, for example, if
   * you try to partition a coarse
   * mesh into more partitions (one
   * for each processor) than there
   * are cells in the mesh.
   *
   * This function returns the number of
   * cells associated with one
   * subdomain. If you are looking for the
   * association of @em DoFs with this
   * subdomain, use the
   * <tt>DoFTools::count_dofs_with_subdomain_association</tt>
   * function.
   */
  template <int dim, int spacedim>
  unsigned int
  count_cells_with_subdomain_association (const Triangulation<dim, spacedim> &triangulation,
                                          const types::subdomain_id         subdomain);

  /*@}*/
  /**
   *  @name Comparing different meshes
   */
  /*@{*/

  /**
   * Given two mesh containers
   * (i.e. objects of type
   * Triangulation, DoFHandler,
   * hp::DoFHandler, or
   * MGDoFHandler) that are based
   * on the same coarse mesh, this
   * function figures out a set of
   * cells that are matched between
   * the two meshes and where at
   * most one of the meshes is more
   * refined on this cell. In other
   * words, it finds the smallest
   * cells that are common to both
   * meshes, and that together
   * completely cover the domain.
   *
   * This function is useful, for
   * example, in time-dependent or
   * nonlinear application, where
   * one has to integrate a
   * solution defined on one mesh
   * (e.g., the one from the
   * previous time step or
   * nonlinear iteration) against
   * the shape functions of another
   * mesh (the next time step, the
   * next nonlinear iteration). If,
   * for example, the new mesh is
   * finer, then one has to obtain
   * the solution on the coarse
   * mesh (mesh_1) and interpolate
   * it to the children of the
   * corresponding cell of
   * mesh_2. Conversely, if the new
   * mesh is coarser, one has to
   * express the coarse cell shape
   * function by a linear
   * combination of fine cell shape
   * functions. In either case, one
   * needs to loop over the finest
   * cells that are common to both
   * triangulations. This function
   * returns a list of pairs of
   * matching iterators to cells in
   * the two meshes that can be
   * used to this end.
   *
   * Note that the list of these
   * iterators is not necessarily
   * ordered, and does also not
   * necessarily coincide with the
   * order in which cells are
   * traversed in one, or both, of
   * the meshes given as arguments.
   *
   * @tparam Container A type that satisfies the
   *   requirements of a mesh container (see @ref GlossMeshAsAContainer).
   */
  template <typename Container>
  std::list<std::pair<typename Container::cell_iterator,
      typename Container::cell_iterator> >
      get_finest_common_cells (const Container &mesh_1,
                               const Container &mesh_2);

  /**
   * Return true if the two
   * triangulations are based on
   * the same coarse mesh. This is
   * determined by checking whether
   * they have the same number of
   * cells on the coarsest level,
   * and then checking that they
   * have the same vertices.
   *
   * The two meshes may have
   * different refinement histories
   * beyond the coarse mesh.
   */
  template <int dim, int spacedim>
  bool
  have_same_coarse_mesh (const Triangulation<dim, spacedim> &mesh_1,
                         const Triangulation<dim, spacedim> &mesh_2);

  /**
   * The same function as above,
   * but working on arguments of
   * type DoFHandler,
   * hp::DoFHandler, or
   * MGDoFHandler. This function is
   * provided to allow calling
   * have_same_coarse_mesh for all
   * types of containers
   * representing triangulations or
   * the classes built on
   * triangulations.
   *
   * @tparam Container A type that satisfies the
   *   requirements of a mesh container (see @ref GlossMeshAsAContainer).
   */
  template <typename Container>
  bool
  have_same_coarse_mesh (const Container &mesh_1,
                         const Container &mesh_2);

  /**
   * @deprecated Use GridGenerator::create_union_triangulation().
   */
  template <int dim, int spacedim>
  void
  create_union_triangulation (const Triangulation<dim, spacedim> &triangulation_1,
                              const Triangulation<dim, spacedim> &triangulation_2,
                              Triangulation<dim, spacedim>       &result)  DEAL_II_DEPRECATED;

  /*@}*/
  /**
   *  @name Dealing with distorted cells
   */
  /*@{*/

  /**
   * Given a triangulation and a
   * list of cells whose children
   * have become distorted as a
   * result of mesh refinement, try
   * to fix these cells up by
   * moving the center node around.
   *
   * The function returns a list of
   * cells with distorted children
   * that couldn't be fixed up for
   * whatever reason. The returned
   * list is therefore a subset of
   * the input argument.
   *
   * For a definition of the
   * concept of distorted cells,
   * see the
   * @ref GlossDistorted "glossary entry".
   * The first argument passed to the
   * current function is typically
   * the exception thrown by the
   * Triangulation::execute_coarsening_and_refinement
   * function.
   */
  template <int dim, int spacedim>
  typename Triangulation<dim,spacedim>::DistortedCellList
  fix_up_distorted_child_cells (const typename Triangulation<dim,spacedim>::DistortedCellList &distorted_cells,
                                Triangulation<dim,spacedim> &triangulation);

  /*@}*/
  /**
   *  @name Extracting and creating patches of cells surrounding a single cell
   */
  /*@{*/


  /**
   * This function returns a list of all the active neighbor cells of the
   * given, active cell.  Here, a neighbor is defined as one having at least
   * part of a face in common with the given cell, but not edge (in 3d) or
   * vertex neighbors (in 2d and 3d).
   *
   * The first element of the returned list is the cell provided as
   * argument. The remaining ones are neighbors: The function loops over all
   * faces of that given cell and checks if that face is not on the boundary of
   * the domain. Then, if the neighbor cell does not have any children (i.e., it
   * is either at the same refinement level as the current cell, or coarser)
   * then this neighbor cell is added to the list of cells. Otherwise, if the
   * neighbor cell is refined and therefore has children, then this function
   * loops over all subfaces of current face adds the neighbors behind these
   * sub-faces to the list to be returned.
   *
   * @tparam Container A type that satisfies the
   *   requirements of a mesh container (see @ref GlossMeshAsAContainer).
   *   In C++, the compiler can not determine the type of <code>Container</code>
   *   from the function call. You need to specify it as an explicit template
   *   argument following the function name.
   * @param cell[in] An iterator pointing to a cell of the mesh container.
   * @return A list of active cells that form the patch around the given cell
   *
   * @note Patches are often used in defining error estimators that require the
   * solution of a local problem on the patch surrounding each of the cells of
   * the mesh. This also requires manipulating the degrees of freedom associated
   * with the cells of a patch. To this end, there are further functions working
   * on patches in namespace DoFTools.
   *
   * @note In the context of a parallel distributed computation, it only makes
   * sense to call this function on locally owned cells. This is because the
   * neighbors of locally owned cells are either locally owned themselves, or
   * ghost cells. For both, we know that these are in fact the real cells of
   * the complete, parallel triangulation. We can also query the degrees of
   * freedom on these.
   *
   * @author Arezou Ghesmati, Wolfgang Bangerth, 2014
   */
  template <class Container>
  std::vector<typename Container::active_cell_iterator>
  get_patch_around_cell(const typename Container::active_cell_iterator &cell);


  /*@}*/
  /**
   *  @name Lower-dimensional meshes for parts of higher-dimensional meshes
   */
  /*@{*/


#ifdef _MSC_VER
  // Microsoft's VC++ has a bug where it doesn't want to recognize that
  // an implementation (definition) of the extract_boundary_mesh function
  // matches a declaration. This can apparently only be avoided by
  // doing some contortion with the return type using the following
  // intermediate type. This is only used when using MS VC++ and uses
  // the direct way of doing it otherwise
  template <template <int,int> class Container, int dim, int spacedim>
  struct ExtractBoundaryMesh
  {
    typedef
    std::map<typename Container<dim-1,spacedim>::cell_iterator,
        typename Container<dim,spacedim>::face_iterator>
        return_type;
  };
#endif

  /**
   * @deprecated Use GridGenerator::extract_boundary_mesh() instead.
   */
  template <template <int,int> class Container, int dim, int spacedim>
#ifndef _MSC_VER
  std::map<typename Container<dim-1,spacedim>::cell_iterator,
      typename Container<dim,spacedim>::face_iterator>
#else
  typename ExtractBoundaryMesh<Container,dim,spacedim>::return_type
#endif
      extract_boundary_mesh (const Container<dim,spacedim> &volume_mesh,
                             Container<dim-1,spacedim>     &surface_mesh,
                             const std::set<types::boundary_id> &boundary_ids
                             = std::set<types::boundary_id>()) DEAL_II_DEPRECATED;

  /*@}*/
  /**
   *  @name Dealing with periodic domains
   */
  /*@{*/

  /**
   * Data type that provides all the information that is needed
   * to create periodicity constraints and a periodic p4est forest
   * with respect to two periodic cell faces
   */
  template<typename CellIterator>
  struct PeriodicFacePair
  {
    CellIterator cell[2];
    unsigned int face_idx[2];
    std::bitset<3> orientation;
  };

  /**
   * An orthogonal equality test for faces.
   *
   * @p face1 and @p face2 are considered equal, if a one to one matching
   * between its vertices can be achieved via an orthogonal equality
   * relation: Two vertices <tt>v_1</tt> and <tt>v_2</tt> are considered
   * equal, if <code> (v_1 + offset) - v_2</code> is parallel to the unit
   * vector in @p direction.
   *
   * If the matching was successful, the _relative_ orientation of @p face1
   * with respect to @p face2 is returned in the bitset @p orientation,
   * where
   * @code
   * orientation[0] -> face_orientation
   * orientation[1] -> face_flip
   * orientation[2] -> face_rotation
   * @endcode
   *
   * In 2D <tt>face_orientation</tt> is always <tt>true</tt>,
   * <tt>face_rotation</tt> is always <tt>false</tt>, and face_flip has the
   * meaning of <tt>line_flip</tt>. More precisely in 3d:
   *
   * <tt>face_orientation</tt>:
   * <tt>true</tt> if @p face1 and @p face2 have the same orientation.
   * Otherwise, the vertex indices of @p face1 match the vertex indices of
   * @p face2 in the following manner:
   *
   * @code
   * face1:           face2:
   *
   * 1 - 3            2 - 3
   * |   |    <-->    |   |
   * 0 - 2            0 - 1
   * @endcode
   *
   * <tt>face_flip</tt>:
   * <tt>true</tt> if the matched vertices are rotated by 180 degrees:
   *
   * @code
   * face1:           face2:
   *
   * 1 - 0            2 - 3
   * |   |    <-->    |   |
   * 3 - 2            0 - 1
   * @endcode
   *
   * <tt>face_rotation</tt>: <tt>true</tt> if the matched vertices are
   * rotated by 90 degrees counterclockwise:
   *
   * @code
   * face1:           face2:
   *
   * 0 - 2            2 - 3
   * |   |    <-->    |   |
   * 1 - 3            0 - 1
   * @endcode
   *
   * and any combination of that...
   * More information on the topic can be found in the
   * @ref GlossFaceOrientation "glossary" article.
   *
   * @author Matthias Maier, 2012
   */
  template<typename FaceIterator>
  bool
  orthogonal_equality (std::bitset<3>     &orientation,
                       const FaceIterator &face1,
                       const FaceIterator &face2,
                       const int          direction,
                       const dealii::Tensor<1,FaceIterator::AccessorType::space_dimension> &offset);


  /**
   * Same function as above, but doesn't return the actual orientation
   */
  template<typename FaceIterator>
  bool
  orthogonal_equality (const FaceIterator &face1,
                       const FaceIterator &face2,
                       const int          direction,
                       const dealii::Tensor<1,FaceIterator::AccessorType::space_dimension> &offset);


  /**
   * This function will collect periodic face pairs on the
   * coarsest mesh level of the given @p container (a Triangulation or
   * DoFHandler) and add them to the vector @p matched_pairs leaving the
   * original contents intact.
   *
   * Define a 'first' boundary as all boundary faces having boundary_id
   * @p b_id1 and a 'second' boundary consisting of all faces belonging
   * to @p b_id2.
   *
   * This function tries to match all faces belonging to the first
   * boundary with faces belonging to the second boundary with the help
   * of orthogonal_equality().
   *
   * The bitset that is returned inside of PeriodicFacePair encodes the
   * _relative_ orientation of the first face with respect to the second
   * face, see the documentation of orthogonal_equality for further details.
   *
   * The @p direction refers to the space direction in which periodicity
   * is enforced.
   *
   * The @p offset is a vector tangential to the faces that is added to the
   * location of vertices of the 'first' boundary when attempting to match
   * them to the corresponding vertices of the 'second' boundary. This can
   * be used to implement conditions such as $u(0,y)=u(1,y+1)$.
   *
   * @tparam Container A type that satisfies the
   *   requirements of a mesh container (see @ref GlossMeshAsAContainer).
   *
   * @note The created std::vector can be used in
   * DoFTools::make_periodicity_constraints and in
   * parallel::distributed::Triangulation::add_periodicity to enforce
   * periodicity algebraically.
   *
   * @note Because elements will be added to @p matched_pairs (and existing
   * entries will be preserved), it is possible to call this function several
   * times with different boundary ids to generate a vector with all periodic
   * pairs.
   *
   * @author Daniel Arndt, Matthias Maier, 2013
   */
  template <typename Container>
  void
  collect_periodic_faces
  (const Container &container,
   const types::boundary_id b_id1,
   const types::boundary_id b_id2,
   const int                direction,
   std::vector<PeriodicFacePair<typename Container::cell_iterator> > &matched_pairs,
   const dealii::Tensor<1,Container::space_dimension> &offset = dealii::Tensor<1,Container::space_dimension>());


  /**
   * This compatibility version of collect_periodic_face_pairs only works
   * on grids with cells in @ref GlossFaceOrientation "standard orientation".
   *
   * Instead of defining a 'first' and 'second' boundary with the help of
   * two boundary_indicators this function defines a 'left' boundary as all
   * faces with local face index <code>2*dimension</code> and boundary
   * indicator @p b_id and, similarly, a 'right' boundary consisting of all
   * face with local face index <code>2*dimension+1</code> and boundary
   * indicator @p b_id.
   *
   * This function will collect periodic face pairs on the coarsest mesh level
   * and add them to @p matched_pairs leaving the original contents intact.
   *
   * @tparam Container A type that satisfies the
   *   requirements of a mesh container (see @ref GlossMeshAsAContainer).
   *
   * @note This version of collect_periodic_face_pairs  will not work on
   * meshes with cells not in @ref GlossFaceOrientation
   * "standard orientation".
   *
   * @author Daniel Arndt, Matthias Maier, 2013
   */
  template <typename Container>
  void
  collect_periodic_faces
  (const Container         &container,
   const types::boundary_id b_id,
   const int                direction,
   std::vector<PeriodicFacePair<typename Container::cell_iterator> > &matched_pairs,
   const dealii::Tensor<1,Container::space_dimension> &offset = dealii::Tensor<1,Container::space_dimension>());

  /*@}*/
  /**
   *  @name Exceptions
   */
  /*@{*/


  /**
   * Exception
   */
  DeclException1 (ExcInvalidNumberOfPartitions,
                  int,
                  << "The number of partitions you gave is " << arg1
                  << ", but must be greater than zero.");
  /**
   * Exception
   */
  DeclException1 (ExcNonExistentSubdomain,
                  int,
                  << "The subdomain id " << arg1
                  << " has no cells associated with it.");
  /**
   * Exception
   */
  DeclException0 (ExcTriangulationHasBeenRefined);

  /**
   * Exception
   */
  DeclException1 (ExcScalingFactorNotPositive,
                  double,
                  << "The scaling factor must be positive, but is " << arg1);
  /**
   * Exception
   */
  template <int N>
  DeclException1 (ExcPointNotFoundInCoarseGrid,
                  Point<N>,
                  << "The point <" << arg1
                  << "> could not be found inside any of the "
                  << "coarse grid cells.");
  /**
   * Exception
   */
  template <int N>
  DeclException1 (ExcPointNotFound,
                  Point<N>,
                  << "The point <" << arg1
                  << "> could not be found inside any of the "
                  << "subcells of a coarse grid cell.");

  /**
   * Exception
   */
  DeclException1 (ExcVertexNotUsed,
                  unsigned int,
                  << "The given vertex " << arg1
                  << " is not used in the given triangulation");


  /*@}*/

} /*namespace GridTools*/



/* ----------------- Template function --------------- */

#ifndef DOXYGEN

namespace GridTools
{
  template <int dim, typename Predicate, int spacedim>
  void transform (const Predicate    &predicate,
                  Triangulation<dim, spacedim> &triangulation)
  {
    std::vector<bool> treated_vertices (triangulation.n_vertices(),
                                        false);

    // loop over all active cells, and
    // transform those vertices that
    // have not yet been touched. note
    // that we get to all vertices in
    // the triangulation by only
    // visiting the active cells.
    typename Triangulation<dim, spacedim>::active_cell_iterator
    cell = triangulation.begin_active (),
    endc = triangulation.end ();
    for (; cell!=endc; ++cell)
      for (unsigned int v=0; v<GeometryInfo<dim>::vertices_per_cell; ++v)
        if (treated_vertices[cell->vertex_index(v)] == false)
          {
            // transform this vertex
            cell->vertex(v) = predicate(cell->vertex(v));
            // and mark it as treated
            treated_vertices[cell->vertex_index(v)] = true;
          };


    // now fix any vertices on hanging nodes so that we don't create any holes
    if (dim==2)
      {
        typename Triangulation<dim,spacedim>::active_cell_iterator
        cell = triangulation.begin_active(),
        endc = triangulation.end();
        for (; cell!=endc; ++cell)
          for (unsigned int face=0; face<GeometryInfo<dim>::faces_per_cell; ++face)
            if (cell->face(face)->has_children() &&
                !cell->face(face)->at_boundary())
              {
                // this line has children
                cell->face(face)->child(0)->vertex(1)
                  = (cell->face(face)->vertex(0) +
                     cell->face(face)->vertex(1)) / 2;
              }
      }
    else if (dim==3)
      {
        typename Triangulation<dim,spacedim>::active_cell_iterator
        cell = triangulation.begin_active(),
        endc = triangulation.end();
        for (; cell!=endc; ++cell)
          for (unsigned int face=0; face<GeometryInfo<dim>::faces_per_cell; ++face)
            if (cell->face(face)->has_children() &&
                !cell->face(face)->at_boundary())
              {
                // this face has hanging nodes
                cell->face(face)->child(0)->vertex(1)
                  = (cell->face(face)->vertex(0) + cell->face(face)->vertex(1)) / 2.0;
                cell->face(face)->child(0)->vertex(2)
                  = (cell->face(face)->vertex(0) + cell->face(face)->vertex(2)) / 2.0;
                cell->face(face)->child(1)->vertex(3)
                  = (cell->face(face)->vertex(1) + cell->face(face)->vertex(3)) / 2.0;
                cell->face(face)->child(2)->vertex(3)
                  = (cell->face(face)->vertex(2) + cell->face(face)->vertex(3)) / 2.0;

                // center of the face
                cell->face(face)->child(0)->vertex(3)
                  = (cell->face(face)->vertex(0) + cell->face(face)->vertex(1)
                     + cell->face(face)->vertex(2) + cell->face(face)->vertex(3)) / 4.0;
              }
      }
  }



  template <class Container>
  std::vector<typename Container::active_cell_iterator>
  get_active_child_cells (const typename Container::cell_iterator &cell)
  {
    std::vector<typename Container::active_cell_iterator> child_cells;

    if (cell->has_children())
      {
        for (unsigned int child=0;
             child<cell->n_children(); ++child)
          if (cell->child (child)->has_children())
            {
              const std::vector<typename Container::active_cell_iterator>
              children = get_active_child_cells<Container> (cell->child(child));
              child_cells.insert (child_cells.end(),
                                  children.begin(), children.end());
            }
          else
            child_cells.push_back (cell->child(child));
      }

    return child_cells;
  }



  template <class Container>
  void
  get_active_neighbors(const typename Container::active_cell_iterator        &cell,
                       std::vector<typename Container::active_cell_iterator> &active_neighbors)
  {
    active_neighbors.clear ();
    for (unsigned int n=0; n<GeometryInfo<Container::dimension>::faces_per_cell; ++n)
      if (! cell->at_boundary(n))
        {
          if (Container::dimension == 1)
            {
              // check children of neighbor. note
              // that in 1d children of the neighbor
              // may be further refined. In 1d the
              // case is simple since we know what
              // children bound to the present cell
              typename Container::cell_iterator
              neighbor_child = cell->neighbor(n);
              if (!neighbor_child->active())
                {
                  while (neighbor_child->has_children())
                    neighbor_child = neighbor_child->child (n==0 ? 1 : 0);

                  Assert (neighbor_child->neighbor(n==0 ? 1 : 0)==cell,
                          ExcInternalError());
                }
              active_neighbors.push_back (neighbor_child);
            }
          else
            {
              if (cell->face(n)->has_children())
                // this neighbor has children. find
                // out which border to the present
                // cell
                for (unsigned int c=0; c<cell->face(n)->number_of_children(); ++c)
                  active_neighbors.push_back (cell->neighbor_child_on_subface(n,c));
              else
                {
                  // the neighbor must be active
                  // himself
                  Assert(cell->neighbor(n)->active(), ExcInternalError());
                  active_neighbors.push_back(cell->neighbor(n));
                }
            }
        }
  }




// declaration of explicit specializations

  template <>
  double
  cell_measure<3>(const std::vector<Point<3> > &all_vertices,
                  const unsigned int (&vertex_indices) [GeometryInfo<3>::vertices_per_cell]);

  template <>
  double
  cell_measure<2>(const std::vector<Point<2> > &all_vertices,
                  const unsigned int (&vertex_indices) [GeometryInfo<2>::vertices_per_cell]);
}

#endif

DEAL_II_NAMESPACE_CLOSE

/*----------------------------   grid_tools.h     ---------------------------*/
/* end of #ifndef __deal2__grid_tools_H */
#endif
/*----------------------------   grid_tools.h     ---------------------------*/
