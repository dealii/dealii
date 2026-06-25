// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
// Copyright (C) 2021 - 2026 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Detailed license information governing the source code and contributions
// can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
//
// -----------------------------------------------------------------------------

#ifndef dealii_cell_id_translator_h
#define dealii_cell_id_translator_h

#include <deal.II/base/exceptions.h>
#include <deal.II/base/utilities.h>

#include <deal.II/grid/cell_id.h>
#include <deal.II/grid/tria_accessor.h>

#include <cstdint>
#include <limits>
#include <vector>

DEAL_II_NAMESPACE_OPEN

namespace internal
{
  /**
   * A class that helps to translate between CellIDs and globally unique
   * indices. The resulting index space might be non-contiguous in the case
   * of locally refined meshes.
   *
   * The following code snippet shows how to set up an IndexSet for active
   * cells, using this class:
   * @code
   * // set up translator
   * CellIDTranslator<dim> translator(tria.n_global_coarse_cells(),
   *                                  tria.n_global_levels());
   *
   * // initialize the index set with the correct size
   * IndexSet index_set(translator.size());
   *
   * // convert IDs to indices and add them to the index set
   * for(const auto & cell : tria.active_cell_iterators())
   *   index_set.add_index(translator.translate(cell));
   * @endcode
   *
   * @note The definition of the indices are not relevant.
   */
  template <int dim>
  class CellIDTranslator
  {
  public:
    /**
     * Constructor taking the number of global coarse cells and number of
     * global levels of a triangulation.
     */
    CellIDTranslator(const types::global_cell_index n_coarse_cells,
                     const types::global_cell_index n_global_levels);

    /**
     * Return maximum number of cells, i.e., in the case the triangulation
     * were globally refined `n_global_levels - 1` times.
     */
    types::global_cell_index
    size() const;

    /**
     * Convert a @p cell of type TriaAccessor, CellAccessor, DoFAccessor, or
     * DoFCellAccessor to an index.
     */
    template <typename Accessor>
    types::global_cell_index
    translate(const TriaIterator<Accessor> &cell) const;

    /**
     * Convert a @p cell of type TriaAccessor, CellAccessor, DoFAccessor, or
     * DoFCellAccessor to an index unique on a level.
     */
    template <typename Accessor>
    static types::global_cell_index
    translate_level(const TriaIterator<Accessor> &cell);

    /**
     * Convert the @p i-th child of @p to an index.
     */
    template <typename Accessor>
    types::global_cell_index
    translate(const TriaIterator<Accessor>  &cell,
              const types::global_cell_index i) const;

    /**
     * Convert an index to a CellId.
     */
    CellId
    to_cell_id(const types::global_cell_index id) const;

  private:
    /**
     * Convert a CellId to a index on that cell's level.
     */
    static types::global_cell_index
    to_level_cell_index(const CellId &cell_id);

    /**
     * Number of global coarse cells.
     */
    const types::global_cell_index n_coarse_cells;

    /**
     * Number of global levels.
     */
    const types::global_cell_index n_global_levels;

    /**
     * Offset of each level in the index space.
     */
    std::vector<types::global_cell_index> tree_sizes;
  };



  template <int dim>
  types::global_cell_index
  CellIDTranslator<dim>::size() const
  {
    return n_coarse_cells *
           (Utilities::pow<types::global_cell_index>(
              ReferenceCells::max_n_children<dim>(), n_global_levels) -
            1);
  }



  template <int dim>
  template <typename Accessor>
  types::global_cell_index
  CellIDTranslator<dim>::translate(const TriaIterator<Accessor> &cell) const
  {
    static_assert(dim == Accessor::dimension &&
                    dim == Accessor::structure_dimension,
                  "The information can only be queried for cells.");

    return translate_level(cell) + tree_sizes[cell->level()];
  }



  template <int dim>
  template <typename Accessor>
  types::global_cell_index
  CellIDTranslator<dim>::translate_level(const TriaIterator<Accessor> &cell)
  {
    static_assert(dim == Accessor::dimension &&
                    dim == Accessor::structure_dimension,
                  "The information can only be queried for cells.");

    return to_level_cell_index(
      CellAccessor<Accessor::dimension, Accessor::space_dimension>(*cell).id());
  }



  template <int dim>
  template <typename Accessor>
  types::global_cell_index
  CellIDTranslator<dim>::translate(const TriaIterator<Accessor>  &cell,
                                   const types::global_cell_index i) const
  {
    static_assert(dim == Accessor::dimension &&
                    dim == Accessor::structure_dimension,
                  "The information can only be queried for cells.");

    return (translate(cell) - tree_sizes[cell->level()]) *
             ReferenceCells::max_n_children<dim>() +
           i + tree_sizes[cell->level() + 1];
  }
} // namespace internal

DEAL_II_NAMESPACE_CLOSE

#endif
