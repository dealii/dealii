// ---------------------------------------------------------------------
//
// Copyright (C) 2000 - 2020 by the deal.II authors
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


#include <deal.II/base/timer.h>
#include <deal.II/base/utilities.h>

#include <deal.II/grid/grid_reordering.h>
#include <deal.II/grid/grid_tools.h>

#include <algorithm>
#include <fstream>
#include <functional>
#include <iostream>
#include <set>

DEAL_II_NAMESPACE_OPEN


// anonymous namespace for internal helper functions
namespace
{
  /**
   * A set of functions that
   * reorder the data from the
   * "current" to the "classic"
   * format of vertex numbering of
   * cells and faces. These functions
   * do the reordering of their
   * arguments in-place.
   */
  void reorder_new_to_old_style(std::vector<CellData<1>> &)
  {}


  void reorder_new_to_old_style(std::vector<CellData<2>> &cells)
  {
    for (auto &cell : cells)
      std::swap(cell.vertices[2], cell.vertices[3]);
  }


  void reorder_new_to_old_style(std::vector<CellData<3>> &cells)
  {
    unsigned int tmp[GeometryInfo<3>::vertices_per_cell];
    for (auto &cell : cells)
      {
        for (const unsigned int i : GeometryInfo<3>::vertex_indices())
          tmp[i] = cell.vertices[i];
        for (const unsigned int i : GeometryInfo<3>::vertex_indices())
          cell.vertices[i] = tmp[GeometryInfo<3>::ucd_to_deal[i]];
      }
  }


  /**
   * And now also in the opposite direction.
   */
  void reorder_old_to_new_style(std::vector<CellData<1>> &)
  {}


  void reorder_old_to_new_style(std::vector<CellData<2>> &cells)
  {
    // just invert the permutation:
    reorder_new_to_old_style(cells);
  }


  void reorder_old_to_new_style(std::vector<CellData<3>> &cells)
  {
    // undo the ordering above
    unsigned int tmp[GeometryInfo<3>::vertices_per_cell];
    for (auto &cell : cells)
      {
        for (const unsigned int i : GeometryInfo<3>::vertex_indices())
          tmp[i] = cell.vertices[i];
        for (const unsigned int i : GeometryInfo<3>::vertex_indices())
          cell.vertices[GeometryInfo<3>::ucd_to_deal[i]] = tmp[i];
      }
  }
} // namespace



template <int dim, int spacedim>
void
GridReordering<dim, spacedim>::reorder_cells(std::vector<CellData<dim>> &cells,
                                             const bool use_new_style_ordering)
{
  Assert(cells.size() != 0,
         ExcMessage("List of elements to orient must have at least one cell"));

  // there is nothing for us to do in 1d
  if (dim == 1)
    return;

  // if necessary, convert to new-style format
  if (use_new_style_ordering == false)
    reorder_old_to_new_style(cells);

  GridTools::consistently_order_cells(cells);

  // and convert back if necessary
  if (use_new_style_ordering == false)
    reorder_new_to_old_style(cells);
}



template <>
void
GridReordering<1>::invert_all_cells_of_negative_grid(
  const std::vector<Point<1>> &,
  std::vector<CellData<1>> &,
  const bool)
{
  // nothing to be done in 1d
}



template <>
void
GridReordering<1, 2>::invert_all_cells_of_negative_grid(
  const std::vector<Point<2>> &,
  std::vector<CellData<1>> &,
  const bool)
{
  // nothing to be done in 1d
}



template <>
void
GridReordering<1, 3>::invert_all_cells_of_negative_grid(
  const std::vector<Point<3>> &,
  std::vector<CellData<1>> &,
  const bool)
{
  // nothing to be done in 1d
}


template <>
void
GridReordering<2>::invert_all_cells_of_negative_grid(
  const std::vector<Point<2>> &all_vertices,
  std::vector<CellData<2>> &   cells,
  const bool                   use_new_style_ordering)
{
  if (!use_new_style_ordering)
    reorder_old_to_new_style(cells);

  GridTools::invert_all_negative_measure_cells(all_vertices, cells);

  if (!use_new_style_ordering)
    reorder_new_to_old_style(cells);
}



template <>
void
GridReordering<2, 3>::invert_all_cells_of_negative_grid(
  const std::vector<Point<3>> &,
  std::vector<CellData<2>> &,
  const bool)
{
  Assert(false, ExcNotImplemented());
}



template <>
void
GridReordering<3>::invert_all_cells_of_negative_grid(
  const std::vector<Point<3>> &all_vertices,
  std::vector<CellData<3>> &   cells,
  const bool                   use_new_style_ordering)
{
  if (!use_new_style_ordering)
    reorder_old_to_new_style(cells);

  GridTools::invert_all_negative_measure_cells(all_vertices, cells);

  if (!use_new_style_ordering)
    reorder_new_to_old_style(cells);
}



/* ------------------------ explicit instantiations ------------------- */
template class GridReordering<1, 1>;
template class GridReordering<1, 2>;
template class GridReordering<1, 3>;
template class GridReordering<2, 2>;
template class GridReordering<2, 3>;
template class GridReordering<3, 3>;

DEAL_II_NAMESPACE_CLOSE
