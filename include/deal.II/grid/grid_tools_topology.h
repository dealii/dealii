// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2023 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#ifndef dealii_grid_tools_topology_h
#define dealii_grid_tools_topology_h

#include <deal.II/base/config.h>

#include <deal.II/base/point.h>
#include <deal.II/base/template_constraints.h>

#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_description.h>

#include <deal.II/lac/dynamic_sparsity_pattern.h>

#include <list>
#include <map>
#include <set>
#include <tuple>
#include <utility>
#include <vector>

DEAL_II_NAMESPACE_OPEN

namespace GridTools
{
  /**
   * @name Querying or modifying topological information
   */
  /** @{ */
  /**
   * Return the arrays that define the coarse mesh of a Triangulation. This
   * function is the inverse of Triangulation::create_triangulation() in the
   * sense that if one called this function on a triangulation, then that
   * triangulation could be recreated by some kind of refinement from the
   * results of this function.
   *
   * The return value is a tuple with the vector of vertices, the vector of
   * cells, and a SubCellData structure. The latter contains additional
   * information about faces and lines. These three objects are exactly
   * the arguments to Triangulation::create_triangulation().
   *
   * This function is useful in cases where one needs to deconstruct a
   * Triangulation or manipulate the numbering of the vertices in some way: an
   * example is GridGenerator::merge_triangulations().
   */
  template <int dim, int spacedim>
  std::
    tuple<std::vector<Point<spacedim>>, std::vector<CellData<dim>>, SubCellData>
    get_coarse_mesh_description(const Triangulation<dim, spacedim> &tria);

  /**
   * Remove vertices that are not referenced by any of the cells. This
   * function is called by all <tt>GridIn::read_*</tt> functions to eliminate
   * vertices that are listed in the input files but are not used by the cells
   * in the input file. While these vertices should not be in the input from
   * the beginning, they sometimes are, most often when some cells have been
   * removed by hand without wanting to update the vertex lists, as they might
   * be lengthy.
   *
   * This function is called by all <tt>GridIn::read_*</tt> functions as the
   * triangulation class requires them to be called with used vertices only.
   * This is so, since the vertices are copied verbatim by that class, so we
   * have to eliminate unused vertices beforehand.
   *
   * Not implemented for the codimension one case.
   */
  template <int dim, int spacedim>
  void
  delete_unused_vertices(std::vector<Point<spacedim>> &vertices,
                         std::vector<CellData<dim>>   &cells,
                         SubCellData                  &subcelldata);

  /**
   * Remove vertices that are duplicated, due to the input of a structured
   * grid, for example. If these vertices are not removed, the faces bounded
   * by these vertices become part of the boundary, even if they are in the
   * interior of the mesh.
   *
   * This function is called by some <tt>GridIn::read_*</tt> functions. Only
   * the vertices with indices in @p considered_vertices are tested for
   * equality. This speeds up the algorithm, which is, for worst-case hyper
   * cube geometries $O(N^{3/2})$ in 2d and $O(N^{5/3})$ in 3d: quite slow.
   * However, if you wish to consider all vertices, simply pass an empty
   * vector. In that case, the function fills @p considered_vertices with all
   * vertices.
   *
   * Two vertices are considered equal if their difference in each coordinate
   * direction is less than @p tol. This implies that nothing happens if
   * the tolerance is set to zero.
   */
  template <int dim, int spacedim>
  void
  delete_duplicated_vertices(std::vector<Point<spacedim>> &all_vertices,
                             std::vector<CellData<dim>>   &cells,
                             SubCellData                  &subcelldata,
                             std::vector<unsigned int>    &considered_vertices,
                             const double                  tol = 1e-12);

  /**
   * Remove vertices that are duplicated.
   *
   * Two vertices are considered equal if their difference in each coordinate
   * direction is less than @p tol. This implies that nothing happens if
   * the tolerance is set to zero.
   */
  template <int dim>
  void
  delete_duplicated_vertices(std::vector<Point<dim>> &vertices,
                             const double             tol = 1e-12);

  /**
   * Grids generated by grid generators may have an orientation of cells which
   * is the inverse of the orientation required by deal.II.
   *
   * In 2d and 3d this function checks whether all cells have negative or
   * positive measure/volume. In the former case, all cells are inverted. It
   * does nothing in 1d.
   *
   * The inversion of cells might also work when only a subset of all cells
   * have negative volume. However, grids consisting of a mixture of negative
   * and positively oriented cells are very likely to be broken. Therefore, an
   * exception is thrown, in case cells are not uniformly oriented.
   *
   * @note This function should be called before GridTools::consistently_order_cells().
   *
   * @param all_vertices The vertices of the mesh.
   * @param cells The array of CellData objects that describe the mesh's topology.
   */
  template <int dim, int spacedim>
  void
  invert_all_negative_measure_cells(
    const std::vector<Point<spacedim>> &all_vertices,
    std::vector<CellData<dim>>         &cells);

  /**
   * Check the given cells and inverts any cell that is considered to have
   * negative measure/volume in the orientation required by deal.II.
   *
   * This function is identical to invert_all_negative_measure_cells() except it
   * does not throw an error if only some of the cells are inverted.  Instead,
   * this function returns how many cells were inverted.  Additionally, it will
   * always throw an exception outside of codimension 0.
   */
  template <int dim, int spacedim>
  std::size_t
  invert_cells_with_negative_measure(
    const std::vector<Point<spacedim>> &all_vertices,
    std::vector<CellData<dim>>         &cells);

  /**
   * Given a vector of CellData objects describing a mesh, reorder their
   * vertices so that all lines are consistently oriented.
   *
   * The expectations on orientation and a discussion of this function are
   * available in the
   * @ref reordering "reordering topic".
   *
   * @param cells The array of CellData objects that describe the mesh's topology.
   * @ingroup reordering
   */
  template <int dim>
  void
  consistently_order_cells(std::vector<CellData<dim>> &cells);

  /**
   * Return a std::map with all vertices of faces located in the boundary
   *
   * @param[in] tria The Triangulation object.
   */
  template <int dim, int spacedim>
  std::map<unsigned int, Point<spacedim>>
  get_all_vertices_at_boundary(const Triangulation<dim, spacedim> &tria);

  /**
   * Return a std::vector that contains a map for each closed boundary.
   * Each closed boundary is described by a std:vector of pairs
   * `vertex index , Point<spacedim>`.
   *
   * The vertices are in counter-clockwise ordering for the outer boundary.
   * The vertices are in clockwise ordering for inner boundaries (holes).
   *
   * It is generally not guaranteed that the first entry of the vector is the
   * outer boundary. However, when cell 0 is located on the outer boundary this
   * is the case.
   *
   * The mapping argument enables the use of e.g. MappingQEulerian.
   *
   * @param[in] tria The Triangulation object.
   * @param[in] mapping The mapping used to map the vertices of the cell
   */
  template <int dim, int spacedim>
  std::vector<std::vector<std::pair<unsigned int, Point<spacedim>>>>
  extract_ordered_boundary_vertices(
    const Triangulation<dim, spacedim> &tria,
    const Mapping<dim, spacedim>       &mapping =
      (ReferenceCells::get_hypercube<dim>()
#ifndef _MSC_VER
         .template get_default_linear_mapping<dim, spacedim>()
#else
         .ReferenceCell::get_default_linear_mapping<dim, spacedim>()
#endif
         ));

  /**
   * Remove hanging nodes from a grid. If the @p isotropic parameter is set
   * to @p false (default) this function detects cells with hanging nodes and
   * refines the neighbours in the direction that removes hanging nodes.
   * If the @p isotropic parameter is set
   * to @p true, the neighbours refinement is made in each directions.
   * In order to remove all hanging nodes this procedure has to be repeated:
   * this could require a large number of iterations.
   * To avoid this a max number (@p max_iterations) of iteration is provided.
   *
   * Consider the following grid:
   * @image html remove_hanging_nodes-hanging.png
   *
   * @p isotropic == @p false would return:
   * @image html remove_hanging_nodes-aniso.png
   *
   * @p isotropic == @p true would return:
   * @image html remove_hanging_nodes-isotro.png
   *
   * @param[in,out] tria Triangulation to refine.
   *
   * @param[in] isotropic If true refine cells in each directions, otherwise
   * (default value) refine the cell in the direction that removes hanging node.
   *
   * @param[in] max_iterations At each step only closest cells to hanging nodes
   * are refined. The code may require a lot of iterations to remove all
   * hanging nodes. @p max_iterations is the maximum number of iteration
   * allowed. If @p max_iterations == numbers::invalid_unsigned_int this
   * function continues refining until there are no hanging nodes.
   *
   * @note In the case of parallel codes, this function should be combined
   * with GridGenerator::flatten_triangulation.
   */
  template <int dim, int spacedim>
  void
  remove_hanging_nodes(Triangulation<dim, spacedim> &tria,
                       const bool                    isotropic      = false,
                       const unsigned int            max_iterations = 100);

  /**
   * Refine a mesh anisotropically such that the resulting mesh is composed by
   * cells with maximum ratio between dimensions less than @p max_ratio.
   * This procedure requires an algorithm that may not terminate. Consequently,
   * it is possible to set a maximum number of iterations through the
   * @p max_iterations parameter.
   *
   * Starting from a cell like this:
   * @image html remove_anisotropy-coarse.png
   *
   * This function would return:
   * @image html remove_anisotropy-refined.png
   *
   * @param[in,out] tria Triangulation to refine.
   *
   * @param[in] max_ratio Maximum value allowed among the ratio between
   * the dimensions of each cell.
   *
   * @param[in] max_iterations Maximum number of iterations allowed.
   *
   * @note In the case of parallel codes, this function should be combined
   * with GridGenerator::flatten_triangulation and
   * GridTools::remove_hanging_nodes.
   */
  template <int dim, int spacedim>
  void
  remove_anisotropy(Triangulation<dim, spacedim> &tria,
                    const double                  max_ratio      = 1.6180339887,
                    const unsigned int            max_iterations = 5);

  /**
   * Return a map `vertex index -> Point<spacedim>` containing the used
   * vertices of the given `container`. The key of the returned map (i.e.,
   * the first element of the pair above) is the global index in the
   * triangulation, whereas the value of each pair is the physical
   * location of the corresponding vertex. The used vertices are obtained by
   * looping over all cells,
   * and querying for each cell where its vertices are through the (optional)
   * `mapping` argument.
   *
   * In serial Triangulation objects and parallel::shared::Triangulation
   * objects, the size of the returned map
   * equals Triangulation::n_used_vertices() (not Triangulation::n_vertices()).
   * Note that in parallel::distributed::Triangulation objects, only vertices in
   * locally owned cells and ghost cells are returned, as for all other vertices
   * their real location might not be known (e.g. for distributed computations
   * using MappingQEulerian).
   *
   * If you use the default `mapping`, the returned map satisfies the following
   * equality:
   *
   * @code
   * const auto used_vertices = extract_used_vertices(tria);
   * auto all_vertices = tria.get_vertices();
   *
   * for(const auto &id_and_v : used_vertices)
   *   all_vertices[id_and_v.first] == id_and_v.second; // true
   * @endcode
   *
   * Notice that the above is not satisfied for mappings that change the
   * location of vertices, like MappingQEulerian.
   *
   * @ref ConceptMeshType "MeshType concept".
   * @param container The container to extract vertices from.
   * @param mapping The mapping to use to compute the points locations.
   */
  template <int dim, int spacedim>
  std::map<unsigned int, Point<spacedim>>
  extract_used_vertices(
    const Triangulation<dim, spacedim> &container,
    const Mapping<dim, spacedim>       &mapping =
      (ReferenceCells::get_hypercube<dim>()
#ifndef _MSC_VER
         .template get_default_linear_mapping<dim, spacedim>()
#else
         .ReferenceCell::get_default_linear_mapping<dim, spacedim>()
#endif
         ));

  /**
   * Return the adjacent cells of all the vertices. If a vertex is also a
   * hanging node, the associated coarse cell is also returned. The vertices
   * are ordered by the vertex index. This is the number returned by the
   * function <code>cell-@>vertex_index()</code>. Notice that only the indices
   * marked in the array returned by
   * Triangulation<dim,spacedim>::get_used_vertices() are used.
   */
  template <int dim, int spacedim>
  std::vector<
    std::set<typename Triangulation<dim, spacedim>::active_cell_iterator>>
  vertex_to_cell_map(const Triangulation<dim, spacedim> &triangulation);

  /**
   * Produce a sparsity pattern in which nonzero entries indicate that two
   * cells are connected via a common face. The diagonal entries of the
   * sparsity pattern are also set.
   *
   * The rows and columns refer to the cells as they are traversed in their
   * natural order using cell iterators.
   */
  template <int dim, int spacedim>
  void
  get_face_connectivity_of_cells(
    const Triangulation<dim, spacedim> &triangulation,
    DynamicSparsityPattern             &connectivity);

  /**
   * Produce a sparsity pattern in which nonzero entries indicate that two
   * cells are connected via a common vertex. The diagonal entries of the
   * sparsity pattern are also set.
   *
   * The rows and columns refer to the cells as they are traversed in their
   * natural order using cell iterators.
   */
  template <int dim, int spacedim>
  void
  get_vertex_connectivity_of_cells(
    const Triangulation<dim, spacedim> &triangulation,
    DynamicSparsityPattern             &connectivity);

  /**
   * Produce a sparsity pattern for a given level mesh in which nonzero entries
   * indicate that two cells are connected via a common vertex. The diagonal
   * entries of the sparsity pattern are also set.
   *
   * The rows and columns refer to the cells as they are traversed in their
   * natural order using cell iterators.
   */
  template <int dim, int spacedim>
  void
  get_vertex_connectivity_of_cells_on_level(
    const Triangulation<dim, spacedim> &triangulation,
    const unsigned int                  level,
    DynamicSparsityPattern             &connectivity);

  /** @} */

  /**
   * @name Comparing different meshes
   */
  /** @{ */

  /**
   * Given two meshes (i.e. objects of type Triangulation or DoFHandler) that
   * are based on the same coarse mesh, this function figures out a set of cells
   * that are matched between the two meshes and where at most one of the meshes
   * is more refined on this cell. In other words, it finds the smallest cells
   * that are common to both meshes, and that together completely cover the
   * domain.
   *
   * This function is useful, for example, in time-dependent or nonlinear
   * application, where one has to integrate a solution defined on one mesh
   * (e.g., the one from the previous time step or nonlinear iteration)
   * against the shape functions of another mesh (the next time step, the next
   * nonlinear iteration). If, for example, the new mesh is finer, then one
   * has to obtain the solution on the coarse mesh (mesh_1) and interpolate it
   * to the children of the corresponding cell of mesh_2. Conversely, if the
   * new mesh is coarser, one has to express the coarse cell shape function by
   * a linear combination of fine cell shape functions. In either case, one
   * needs to loop over the finest cells that are common to both
   * triangulations. This function returns a list of pairs of matching
   * iterators to cells in the two meshes that can be used to this end.
   *
   * Note that the list of these iterators is not necessarily ordered, and
   * does also not necessarily coincide with the order in which cells are
   * traversed in one, or both, of the meshes given as arguments.
   *
   * @tparam MeshType A type that satisfies the requirements of the
   * @ref ConceptMeshType "MeshType concept".
   *
   * @note This function can only be used with
   * parallel::distributed::Triangulation when both meshes use the same
   * Triangulation since, with a distributed Triangulation, not all cells are
   * stored locally, so the resulting list may not cover the entire domain.
   *
   * @dealiiConceptRequires{concepts::is_triangulation_or_dof_handler<MeshType>}
   */
  template <typename MeshType>
  DEAL_II_CXX20_REQUIRES(concepts::is_triangulation_or_dof_handler<MeshType>)
  std::list<std::pair<
    typename MeshType::cell_iterator,
    typename MeshType::cell_iterator>> get_finest_common_cells(const MeshType
                                                                 &mesh_1,
                                                               const MeshType
                                                                 &mesh_2);

  /**
   * Return true if the two triangulations are based on the same coarse mesh.
   * This is determined by checking whether they have the same number of cells
   * on the coarsest level, and then checking that they have the same
   * vertices.
   *
   * The two meshes may have different refinement histories beyond the coarse
   * mesh.
   */
  template <int dim, int spacedim>
  bool
  have_same_coarse_mesh(const Triangulation<dim, spacedim> &mesh_1,
                        const Triangulation<dim, spacedim> &mesh_2);

  /**
   * The same function as above, but working on arguments of type DoFHandler.
   * This function is provided to allow calling have_same_coarse_mesh for all
   * types of containers representing triangulations or the classes built on
   * triangulations.
   *
   * @tparam MeshType A type that satisfies the requirements of the
   * @ref ConceptMeshType "MeshType concept".
   *
   * @dealiiConceptRequires{concepts::is_triangulation_or_dof_handler<MeshType>}
   */
  template <typename MeshType>
  DEAL_II_CXX20_REQUIRES(concepts::is_triangulation_or_dof_handler<MeshType>)
  bool have_same_coarse_mesh(const MeshType &mesh_1, const MeshType &mesh_2);

  /** @} */

  /**
   * @name Exceptions
   */
  /** @{ */

  /**
   * An exception that is thrown whenever the edges of a mesh are not
   * orientable.
   *
   * @ingroup Exceptions
   */
  DeclExceptionMsg(ExcMeshNotOrientable,
                   "The edges of the mesh are not consistently orientable.");

  /** @} */
} // namespace GridTools

DEAL_II_NAMESPACE_CLOSE

#endif
