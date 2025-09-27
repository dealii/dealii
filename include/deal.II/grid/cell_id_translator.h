// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2021 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#ifndef dealii_cell_id_translator_h
#define dealii_cell_id_translator_h

#include <deal.II/base/exceptions.h>
#include <deal.II/base/geometry_info.h>

#include <deal.II/grid/cell_id.h>
#include <deal.II/grid/tria_accessor.h>

#include <cstdint>
#include <limits>

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
     * Convert a binary representation to a index.
     */
    static types::global_cell_index
    convert_cell_id_binary_type_to_level_coarse_cell_id(
      const typename CellId::binary_type &binary_representation);

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
  CellIDTranslator<dim>::CellIDTranslator(
    const types::global_cell_index n_coarse_cells,
    const types::global_cell_index n_global_levels)
    : n_coarse_cells(n_coarse_cells)
    , n_global_levels(n_global_levels)
  {
    // The class stores indices as types::global_cell_index variables,
    // but when configuring deal.II with default flags, this is a 32-bit
    // data type and it is possible with highly (locally) refined meshes
    // that we exceed the maximal 32-bit numbers even with relatively
    // modest numbers of cells. Check for this by first calculating
    // the maximal index we will get in 64-bit arithmetic and testing
    // that it is representable in 32-bit arithmetic:
    std::uint64_t max_cell_index = 0;

    for (unsigned int i = 0; i < n_global_levels; ++i)
      max_cell_index +=
        Utilities::pow<std::uint64_t>(GeometryInfo<dim>::max_children_per_cell,
                                      i) *
        n_coarse_cells;

    max_cell_index -= 1;

    AssertThrow(
      max_cell_index <= std::numeric_limits<types::global_cell_index>::max(),
      ExcMessage(
        "You have exceeded the maximal number of possible indices this function "
        "can handle. The current setup (n_coarse_cells=" +
        std::to_string(n_coarse_cells) +
        ", n_global_levels=" + std::to_string(n_global_levels) + ") requires " +
        std::to_string(max_cell_index + 1) +
        " indices but the current deal.II configuration only supports " +
        std::to_string(std::numeric_limits<types::global_cell_index>::max()) +
        " indices. You may want to consider to build deal.II with 64bit "
        "indices (-D DEAL_II_WITH_64BIT_INDICES=\"ON\") to increase the limit "
        "of indices."));

    // Now do the whole computation again, but for real:
    tree_sizes.reserve(n_global_levels + 1);
    tree_sizes.push_back(0);
    for (unsigned int i = 0; i < n_global_levels; ++i)
      tree_sizes.push_back(tree_sizes.back() +
                           Utilities::pow<types::global_cell_index>(
                             GeometryInfo<dim>::max_children_per_cell, i) *
                             n_coarse_cells);
  }



  template <int dim>
  types::global_cell_index
  CellIDTranslator<dim>::size() const
  {
    return n_coarse_cells *
           (Utilities::pow<types::global_cell_index>(
              GeometryInfo<dim>::max_children_per_cell, n_global_levels) -
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

    return convert_cell_id_binary_type_to_level_coarse_cell_id(
      CellAccessor<Accessor::dimension, Accessor::space_dimension>(*cell)
        .id()
        .template to_binary<dim>());
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
             GeometryInfo<dim>::max_children_per_cell +
           i + tree_sizes[cell->level() + 1];
  }



  template <int dim>
  CellId
  CellIDTranslator<dim>::to_cell_id(const types::global_cell_index id) const
  {
    std::vector<std::uint8_t> child_indices;

    types::global_cell_index id_temp = id;

    types::global_cell_index level = 0;

    for (; level < n_global_levels; ++level)
      if (id < tree_sizes[level])
        break;
    level -= 1;

    id_temp -= tree_sizes[level];

    for (types::global_cell_index l = 0; l < level; ++l)
      {
        child_indices.push_back(id_temp %
                                GeometryInfo<dim>::max_children_per_cell);
        id_temp /= GeometryInfo<dim>::max_children_per_cell;
      }

    std::reverse(child_indices.begin(), child_indices.end());

    return {id_temp, child_indices};
  }



  template <int dim>
  types::global_cell_index
  CellIDTranslator<dim>::convert_cell_id_binary_type_to_level_coarse_cell_id(
    const typename CellId::binary_type &binary_representation)
  {
    // exploiting the structure of CellId::binary_type
    // see also the documentation of CellId

    // actual coarse-grid id
    const unsigned int coarse_cell_id  = binary_representation[0];
    const unsigned int n_child_indices = binary_representation[1] >> 2;

    const unsigned int children_per_value =
      sizeof(CellId::binary_type::value_type) * 8 / dim;
    unsigned int child_level  = 0;
    unsigned int binary_entry = 2;

    // compute new coarse-grid id: c_{i+1} = c_{i}*2^dim + q on path to cell
    types::global_cell_index level_coarse_cell_id = coarse_cell_id;
    while (child_level < n_child_indices)
      {
        Assert(binary_entry < binary_representation.size(), ExcInternalError());

        for (unsigned int j = 0; j < children_per_value; ++j)
          {
            unsigned int cell_index =
              (((binary_representation[binary_entry] >> (j * dim))) &
               (GeometryInfo<dim>::max_children_per_cell - 1));
            level_coarse_cell_id =
              level_coarse_cell_id * GeometryInfo<dim>::max_children_per_cell +
              cell_index;
            ++child_level;
            if (child_level == n_child_indices)
              break;
          }
        ++binary_entry;
      }

    return level_coarse_cell_id;
  }

} // namespace internal

DEAL_II_NAMESPACE_CLOSE

#endif
