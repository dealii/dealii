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

#include <deal.II/base/exceptions.h>
#include <deal.II/base/std_cxx26/inplace_vector.h>
#include <deal.II/base/utilities.h>

#include <deal.II/grid/cell_id_translator.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>

#include <algorithm>
#include <cstdint>
#include <limits>

DEAL_II_NAMESPACE_OPEN

namespace internal
{
  template <int dim>
  template <int spacedim>
  CellIDTranslator<dim>::CellIDTranslator(
    const Triangulation<dim, spacedim> &tria)
    : n_coarse_cells(tria.n_global_coarse_cells())
    , n_global_levels(tria.n_global_levels())
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
        Utilities::pow<std::uint64_t>(ReferenceCells::max_n_children<dim>(),
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
                             ReferenceCells::max_n_children<dim>(), i) *
                             n_coarse_cells);
  }



  template <int dim>
  CellId
  CellIDTranslator<dim>::to_cell_id(const types::global_cell_index id) const
  {
    std_cxx26::inplace_vector<std::uint8_t, numbers::max_n_levels - 1>
      child_indices;

    std::uint8_t level = 0;
    for (; level < n_global_levels; ++level)
      if (id < tree_sizes[level])
        break;
    Assert(level > 0, ExcInternalError());
    level -= 1;

    types::coarse_cell_id id_temp = id - tree_sizes[level];
    for (std::uint8_t l = 0; l < level; ++l)
      {
        child_indices.push_back(id_temp %
                                ReferenceCells::max_n_children<dim>());
        id_temp /= ReferenceCells::max_n_children<dim>();
      }

    std::reverse(child_indices.begin(), child_indices.end());

    return CellId(id_temp,
                  static_cast<unsigned int>(child_indices.size()),
                  child_indices.data());
  }



  template <int dim>
  types::global_cell_index
  CellIDTranslator<dim>::to_level_cell_index(const CellId &cell_id)
  {
    // compute level id: c_{i+1} = c_{i}*(max_n_children) + q on path to cell
    auto level_cell_id = cell_id.get_coarse_cell_id();
    for (const auto &child_index : cell_id.get_child_indices())
      level_cell_id =
        level_cell_id * ReferenceCells::max_n_children<dim>() + child_index;

    return level_cell_id;
  }

} // namespace internal


// explicit instantiations
#include "grid/cell_id_translator.inst"

DEAL_II_NAMESPACE_CLOSE
