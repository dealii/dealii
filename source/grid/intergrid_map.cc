// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 1999 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#include <deal.II/base/memory_consumption.h>
#include <deal.II/base/observer_pointer.h>

#include <deal.II/distributed/fully_distributed_tria.h>
#include <deal.II/distributed/shared_tria.h>
#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe.h>

#include <deal.II/grid/intergrid_map.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>


DEAL_II_NAMESPACE_OPEN


template <typename MeshType>
DEAL_II_CXX20_REQUIRES(concepts::is_triangulation_or_dof_handler<MeshType>)
InterGridMap<MeshType>::InterGridMap()
  : source_grid(nullptr, typeid(*this).name())
  , destination_grid(nullptr, typeid(*this).name())
{}



template <typename MeshType>
DEAL_II_CXX20_REQUIRES(concepts::is_triangulation_or_dof_handler<MeshType>)
void InterGridMap<MeshType>::make_mapping(const MeshType &source_grid,
                                          const MeshType &destination_grid)
{
  // first delete all contents
  clear();

  // next store pointers to grids
  this->source_grid      = &source_grid;
  this->destination_grid = &destination_grid;

  // then set up the meshes from
  // scratch and fill them with end-iterators
  const unsigned int n_levels = source_grid.get_triangulation().n_levels();
  mapping.resize(n_levels);
  for (unsigned int level = 0; level < n_levels; ++level)
    {
      // first find out about the highest
      // index used on this level. We could
      // in principle ask the triangulation
      // about this, but we would have to
      // know the underlying data structure
      // for this and we would like to
      // avoid such knowledge here
      unsigned int  n_cells = 0;
      cell_iterator cell    = source_grid.begin(level),
                    endc    = source_grid.end(level);
      for (; cell != endc; ++cell)
        if (static_cast<unsigned int>(cell->index()) > n_cells)
          n_cells = cell->index();

      // note: n_cells is now the largest
      // zero-based index, but we need the
      // number of cells, which is one larger
      mapping[level].resize(n_cells + 1, destination_grid.end());
    }

  // now make up the mapping
  // loop over all cells and set the user
  // pointers as well as the contents of
  // the two arrays. note that the function
  // takes a *reference* to the int and
  // this may change it
  cell_iterator src_cell = source_grid.begin(0),
                dst_cell = destination_grid.begin(0), endc = source_grid.end(0);
  for (; src_cell != endc; ++src_cell, ++dst_cell)
    set_mapping(src_cell, dst_cell);

  // little assertion that the two grids
  // are indeed related:
  Assert(dst_cell == destination_grid.end(0), ExcIncompatibleGrids());
}



template <typename MeshType>
DEAL_II_CXX20_REQUIRES(concepts::is_triangulation_or_dof_handler<MeshType>)
void InterGridMap<MeshType>::set_mapping(const cell_iterator &src_cell,
                                         const cell_iterator &dst_cell)
{
  // first set the map for this cell
  mapping[src_cell->level()][src_cell->index()] = dst_cell;

  // if both cells have children, we may
  // recurse further into the hierarchy
  if (src_cell->has_children() && dst_cell->has_children())
    {
      Assert(src_cell->n_children() ==
               GeometryInfo<MeshType::dimension>::max_children_per_cell,
             ExcNotImplemented());
      Assert(dst_cell->n_children() ==
               GeometryInfo<MeshType::dimension>::max_children_per_cell,
             ExcNotImplemented());
      Assert(src_cell->refinement_case() == dst_cell->refinement_case(),
             ExcNotImplemented());
      for (unsigned int c = 0;
           c < GeometryInfo<MeshType::dimension>::max_children_per_cell;
           ++c)
        set_mapping(src_cell->child(c), dst_cell->child(c));
    }
  else if (src_cell->has_children() && !dst_cell->has_children())
    // src grid is more refined here.
    // set entries for all children
    // of this cell to the one
    // dst_cell
    for (unsigned int c = 0; c < src_cell->n_children(); ++c)
      set_entries_to_cell(src_cell->child(c), dst_cell);
  // else (no cell is refined or
  // dst_cell is refined): no pointers
  // to be set
}



template <typename MeshType>
DEAL_II_CXX20_REQUIRES(concepts::is_triangulation_or_dof_handler<MeshType>)
void InterGridMap<MeshType>::set_entries_to_cell(const cell_iterator &src_cell,
                                                 const cell_iterator &dst_cell)
{
  // first set the map for this cell
  mapping[src_cell->level()][src_cell->index()] = dst_cell;

  // then do so for the children as well
  // if there are any
  if (src_cell->has_children())
    for (unsigned int c = 0; c < src_cell->n_children(); ++c)
      set_entries_to_cell(src_cell->child(c), dst_cell);
}


template <typename MeshType>
DEAL_II_CXX20_REQUIRES(concepts::is_triangulation_or_dof_handler<MeshType>)
typename InterGridMap<MeshType>::cell_iterator
  InterGridMap<MeshType>::operator[](const cell_iterator &source_cell) const
{
  Assert(source_cell.state() == IteratorState::valid,
         ExcInvalidKey(source_cell));
  Assert(source_cell->level() <= static_cast<int>(mapping.size()),
         ExcInvalidKey(source_cell));
  Assert(source_cell->index() <=
           static_cast<int>(mapping[source_cell->level()].size()),
         ExcInvalidKey(source_cell));

  return mapping[source_cell->level()][source_cell->index()];
}



template <typename MeshType>
DEAL_II_CXX20_REQUIRES(concepts::is_triangulation_or_dof_handler<MeshType>)
void InterGridMap<MeshType>::clear()
{
  mapping.clear();
  source_grid      = nullptr;
  destination_grid = nullptr;
}



template <typename MeshType>
DEAL_II_CXX20_REQUIRES(concepts::is_triangulation_or_dof_handler<MeshType>)
const MeshType &InterGridMap<MeshType>::get_source_grid() const
{
  return *source_grid;
}



template <typename MeshType>
DEAL_II_CXX20_REQUIRES(concepts::is_triangulation_or_dof_handler<MeshType>)
const MeshType &InterGridMap<MeshType>::get_destination_grid() const
{
  return *destination_grid;
}



template <typename MeshType>
DEAL_II_CXX20_REQUIRES(concepts::is_triangulation_or_dof_handler<MeshType>)
std::size_t InterGridMap<MeshType>::memory_consumption() const
{
  return (MemoryConsumption::memory_consumption(mapping) +
          MemoryConsumption::memory_consumption(source_grid) +
          MemoryConsumption::memory_consumption(destination_grid));
}



// explicit instantiations
#include "grid/intergrid_map.inst"

DEAL_II_NAMESPACE_CLOSE
