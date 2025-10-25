// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



#include <deal.II/matrix_free/patch_storage.h>

DEAL_II_NAMESPACE_OPEN

namespace internal
{

  void
  orient_patch2D(
    std::vector<typename Triangulation<2>::cell_iterator> &patch_cells,
    const types::global_vertex_index                      &vertex)
  {
    const static unsigned int lookup_rotations[4] = {0, 3, 1, 2};
    const unsigned int        n_rotations =
      lookup_rotations[internal::compute_vertex_index<2>(patch_cells[0],
                                                         vertex)];
    for (unsigned int i = 0; i < n_rotations; ++i)
      rotate_patch<2>(patch_cells);
  }

  template <>
  std::vector<typename Triangulation<2>::cell_iterator>
  order_patch<2>(
    const std::set<typename Triangulation<2>::cell_iterator> &patch_cells,
    const types::global_vertex_index                         &vertex_index)
  {
    const constexpr int                 dim = 2;
    const static constexpr unsigned int n_vertices =
      GeometryInfo<dim>::vertices_per_cell;


    // in case of non-standard patch just copy cell into vector
    if (patch_cells.size() != n_vertices)
      return std::vector<typename Triangulation<2>::cell_iterator>(
        patch_cells.begin(), patch_cells.end());


    std::set<typename Triangulation<dim>::cell_iterator> patch_copy =
      patch_cells;
    // we have a regular patch, thus  we need to order the cells:
    std::vector<typename Triangulation<dim>::cell_iterator> ordered_cells;
    ordered_cells.reserve(n_vertices);
    ordered_cells.push_back(*patch_copy.begin());
    patch_copy.erase(patch_copy.begin());


    // Look-up table. If face i is on the left of the cell
    // right2upper[i] tells which face is on the upper side of cell
    const constexpr std::array<int, n_vertices> right2upper = {{3, 2, 0, 1}};

    for (unsigned int i = 0; i < n_vertices; ++i)
      {
        if (patch_copy.find(ordered_cells[0]->neighbor(i)) !=
              patch_copy.end() &&
            patch_copy.find(ordered_cells[0]->neighbor(right2upper[i])) !=
              patch_copy.end())
          {
            ordered_cells.push_back(ordered_cells[0]->neighbor(right2upper[i]));
            patch_copy.erase(ordered_cells[0]->neighbor(right2upper[i]));

            ordered_cells.push_back(ordered_cells[0]->neighbor(i));
            patch_copy.erase(ordered_cells[0]->neighbor(i));
          }
      }

    AssertDimension(patch_copy.size(), 1);

    ordered_cells.push_back(*patch_copy.begin());

    orient_patch2D(ordered_cells, vertex_index);
    return ordered_cells;
  }


  // counter-clockwise rotation of the given patch
  template <>
  void
  rotate_patch<2>(
    std::vector<typename Triangulation<2>::cell_iterator> &patch_cells)
  {
    static const std::vector<size_t> patch_rotation = {2, 0, 3, 1};

    internal::reorder(patch_cells, patch_rotation);
  }


  namespace GaussSeidel
  {
    unsigned int
    min_category(const unsigned int a, const unsigned int b)
    {
      if ((a == 1 && b == 2) && (b == 1 && a == 2))
        return 0;

      return std::min(a, b);
    }
  } // namespace GaussSeidel

} // namespace internal

DEAL_II_NAMESPACE_CLOSE
