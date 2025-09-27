// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 1999 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#ifndef dealii_intergrid_map_h
#define dealii_intergrid_map_h

#include <deal.II/base/config.h>

#include <deal.II/base/observer_pointer.h>

#include <deal.II/dofs/dof_accessor.h>

#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

DEAL_II_NAMESPACE_OPEN


/**
 * This class provides a map between two grids which are derived from the same
 * coarse grid. For each cell iterator of the source map, it provides the
 * respective cell iterator on the destination map, through its <tt>operator
 * []</tt>.
 *
 * Usually, the two grids will be refined differently. Then, the value
 * returned for an iterator on the source grid will be either:
 * <ul>
 * <li> The same cell on the destination grid, if it exists there;
 * <li> The most refined cell of the destination grid from which the pendant
 * of the source cell could be obtained by refinement. This cell is always
 * active and has a refinement level less than that of the source cell.
 * </ul>
 * Keys for this map are all cells on the source grid, whether active or not.
 *
 * For example, consider these two one-dimensional grids:
 * @verbatim
 * Grid 1:
 *   x--x--x-----x-----------x
 *    1  2    3        4
 *
 * Grid 2:
 *   x-----x-----x-----x-----x
 *      1     2     3     4
 * @endverbatim
 * (Cell numbers are only given as an example and will not correspond to real
 * cell iterator's indices.) The mapping from grid 1 to grid 2 will then be as
 * follows:
 * @verbatim
 *    Cell on grid 1         Cell on grid 2
 *          1  ------------------>  1
 *          2  ------------------>  1
 *          3  ------------------>  2
 *          4  ------------------>  parent cell of cells 3 and 4
 *                                  (a non-active cell, not shown here)
 * @endverbatim
 * Besides the mappings shown here, the non-active cells on grid 1 are also
 * valid keys. For example, the mapping for the parent cell of cells 1 and 2
 * on the first grid will point to cell 1 on the second grid.
 *
 * @tparam MeshType This class may be used with any class that satisfies the
 * @ref ConceptMeshType "MeshType concept".
 * The extension to other classes offering iterator functions and some minor
 * additional requirements is simple.
 *
 * Note that this class could in principle be based on the C++
 * <tt>std::map<Key,Value></tt> data type. Instead, it uses another data
 * format which is more effective both in terms of computing time for access
 * as well as with regard to memory consumption.
 *
 *
 * <h3>Usage</h3>
 *
 * In practice, use of this class is as follows:
 * @code
 *   // have two grids, which are derived from the same coarse grid
 *   Triangulation<dim> tria1, tria2;
 *   DoFHandler<dim> dof_handler_1 (tria1), dof_handler_2 (tria2);
 *   ...
 *   // do something with these objects, e.g. refine the triangulations
 *   // differently, distribute degrees of freedom, etc
 *   ...
 *   // create the mapping
 *   InterGridMap<DoFHandler<dim> > grid_1_to_2_map;
 *   grid_1_to_2_map.make_mapping (dof_handler_1,
 *                                 dof_handler_2);
 *   ...
 *   typename DoFHandler<dim>::cell_iterator cell = dof_handler_1.begin(),
 *                                           endc = dof_handler_1.end();
 *   for (; cell!=endc; ++cell)
 *     // now do something with the cell of dof_handler_2 corresponding to
 *     // cell (which is one of dof_handler_1's cells)
 *     f (grid_1_to_2_map[cell]);
 * @endcode
 *
 * Note that the template parameters to this class have to be given as
 * <tt>InterGridMap<DoFHandler<2> ></tt>, which here is DoFHandler (and could
 * equally well be Triangulation or PersistentTriangulation).
 *
 * @ingroup grid
 *
 * @dealiiConceptRequires{concepts::is_triangulation_or_dof_handler<MeshType>}
 */
template <typename MeshType>
DEAL_II_CXX20_REQUIRES(concepts::is_triangulation_or_dof_handler<MeshType>)
class InterGridMap : public EnableObserverPointer
{
public:
  /**
   * Typedef to the iterator type of the grid class under consideration.
   */
  using cell_iterator = typename MeshType::cell_iterator;

  /**
   * Constructor setting the class name arguments in the ObserverPointer
   * members.
   */
  InterGridMap();

  /**
   * Create the mapping between the two grids.
   */
  void
  make_mapping(const MeshType &source_grid, const MeshType &destination_grid);

  /**
   * Access operator: give a cell on the source grid and receive the
   * respective cell on the other grid, or if that does not exist, the most
   * refined cell of which the source cell would be created if it were further
   * refined.
   */
  cell_iterator
  operator[](const cell_iterator &source_cell) const;

  /**
   * Delete all data of this class.
   */
  void
  clear();

  /**
   * Return a reference to the source grid.
   */
  const MeshType &
  get_source_grid() const;

  /**
   * Return a reference to the destination grid.
   */
  const MeshType &
  get_destination_grid() const;

  /**
   * Determine an estimate for the memory consumption (in bytes) of this
   * object.
   */
  std::size_t
  memory_consumption() const;

  /**
   * Exception
   */
  DeclException1(ExcInvalidKey,
                 cell_iterator,
                 << "The iterator " << arg1 << " is not valid as key for "
                 << "this map.");
  /**
   * Exception
   */
  DeclException0(ExcIncompatibleGrids);

private:
  /**
   * The actual data. Hold one iterator for each cell on each level.
   */
  std::vector<std::vector<cell_iterator>> mapping;

  /**
   * Store a pointer to the source grid.
   */
  ObserverPointer<const MeshType, InterGridMap<MeshType>> source_grid;

  /**
   * Likewise for the destination grid.
   */
  ObserverPointer<const MeshType, InterGridMap<MeshType>> destination_grid;

  /**
   * Set the mapping for the pair of cells given. These shall match in level
   * of refinement and all other properties.
   */
  void
  set_mapping(const cell_iterator &src_cell, const cell_iterator &dst_cell);

  /**
   * Set the value of the key @p src_cell to @p dst_cell. Do so as well for
   * all the children and their children of @p src_cell. This function is used
   * for cells which are more refined on @p src_grid than on @p dst_grid; then
   * all values of the hierarchy of cells and their children point to one cell
   * on the @p dst_grid.
   */
  void
  set_entries_to_cell(const cell_iterator &src_cell,
                      const cell_iterator &dst_cell);
};


DEAL_II_NAMESPACE_CLOSE

#endif
