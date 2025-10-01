// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2020 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

// check basic functionality of PatchStorage: initialization on a mesh with a
// single patch

#include <deal.II/base/numbers.h>
#include <deal.II/base/point.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_renumbering.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>

#include <deal.II/matrix_free/fe_evaluation.h>

#include <iostream>
#include <limits>
#include <map>
#include <numeric>
#include <set>
#include <vector>

#include "../tests.h"

template <int dim>
void
generate(unsigned int fe_degre)
{
  using IndexType         = unsigned int;
  const IndexType n_cells = 1 << dim;


// Choose output target at compile time. Define USE_DEALLOG to write to deallog,
// otherwise output goes to std::cout.
#ifndef USE_STDCOUT
  auto &output = deallog;
#else
  auto &output = std::cout;
#endif

  output << "=========================================" << std::endl;
  output << "===========START GENERATION==============" << std::endl;
  output << "=========================================" << std::endl;

  output << std::endl;
  output << std::endl;


  Triangulation<dim> triangulation;
  GridGenerator::hyper_cube(triangulation, 0.0, 1.0);
  triangulation.refine_global(1);

  FE_Q<dim>       fe(fe_degre);
  DoFHandler<dim> dof_handler(triangulation);


  dof_handler.distribute_dofs(fe);

  output << "Number of DoFs: " << dof_handler.n_dofs() << std::endl;
  DoFRenumbering::lexicographic(dof_handler);


  std::array<std::vector<types::global_dof_index>, n_cells> cell2patch;
  auto                                                      dof_numbering =
    FEEvaluation<dim, -1, -1, 1, double>(fe,
                                         QGauss<1>(fe_degre + 1),
                                         update_default)
      .get_internal_dof_numbering();


  output << "Internal DoF numbering: ";
  for (const auto &index : dof_numbering)
    output << index << " ";
  output << std::endl;


  const auto boundary_dofs_set = DoFTools::extract_boundary_dofs(dof_handler);
  IndexType  n_boundary_dofs   = 0;

  std::vector<bool> selected_dofs(dof_handler.n_dofs(), false);
  for (const auto &dof_index : boundary_dofs_set)
    {
      AssertIndexRange(dof_index, selected_dofs.size());
      selected_dofs[dof_index] = true;
      ++n_boundary_dofs;
    }
  DoFRenumbering::sort_selected_dofs_back(dof_handler, selected_dofs);


  const IndexType n_interior_dofs = dof_handler.n_dofs() - n_boundary_dofs;

  output << "Number of interior DoFs: " << n_interior_dofs << std::endl;
  output << "Number of boundary DoFs: " << n_boundary_dofs << std::endl;


  auto parent_cell = dof_handler.begin(0);
  for (IndexType c = 0; c < GeometryInfo<dim>::max_children_per_cell; ++c)
    {
      typename DoFHandler<dim>::active_cell_iterator cell =
        parent_cell->child(c);
      std::vector<types::global_dof_index> local_dof_indices(
        fe.n_dofs_per_cell());
      cell->get_dof_indices(local_dof_indices);

      cell2patch[c].resize(fe.n_dofs_per_cell());
      for (IndexType i = 0; i < fe.n_dofs_per_cell(); ++i)
        cell2patch[c][i] = local_dof_indices[dof_numbering[i]];

      output << "DoF indices for cell " << c << " after renumbering: ";
      for (const auto &index : cell2patch[c])
        output << index << " ";
      output << std::endl;
    }


  const std::map<types::global_dof_index, Point<dim>> dof_location_map =
    DoFTools::map_dofs_to_support_points(MappingQ1<dim>(), dof_handler);

  std::ofstream dof_location_file("dof_locations.gnuplot");
  DoFTools::write_gnuplot_dof_support_point_info(dof_location_file,
                                                 dof_location_map);


  // Data structures for lookup tables
  std::array<std::vector<IndexType>, n_cells> skipped_dofs;

  std::array<IndexType, n_cells> n_skipped_dofs;
  // Populate skipped_dofs based on boundary DoFs
  for (IndexType c = 0; c < n_cells; ++c)
    {
      for (IndexType i = 0; i < cell2patch[c].size(); ++i)
        {
          const types::global_dof_index dof_index = cell2patch[c][i];
          if (dof_index >= n_interior_dofs)
            skipped_dofs[c].push_back(i); // Store the local index i
        }
      // Optional: Sort skipped_dofs for consistency if needed
      std::sort(skipped_dofs[c].begin(), skipped_dofs[c].end());
      n_skipped_dofs[c] = skipped_dofs[c].size();
    }

  // {cell_dof, patch_dof, overplap}
  std::array<std::vector<std::tuple<IndexType, IndexType, bool>>, n_cells>
                                 cell_to_patch;
  std::array<IndexType, n_cells> n_cell_to_patch;
  std::array<IndexType, n_cells> n_non_overlapping;

  // Keep track of patch DoFs that have already been assigned to a cell
  // to mark overlaps. The index corresponds to the patch_dof index (which
  // is the global interior dof index).
  std::vector<bool> patch_dof_seen(n_interior_dofs, false);

  for (IndexType c = 0; c < n_cells; ++c)
    {
      for (IndexType i = 0; i < fe.n_dofs_per_cell(); ++i)
        {
          const types::global_dof_index global_dof = cell2patch[c][i];
          // Only consider interior DoFs for the patch mapping
          if (global_dof < n_interior_dofs)
            {
              // The patch DoF index is the global DoF index itself for interior
              // DoFs after renumbering.
              const IndexType patch_dof = global_dof;

              // Check if this patch_dof has been seen before.
              // The first time it's seen, is_overlap is false.
              // Subsequent times, is_overlap is true.
              const bool is_overlap = patch_dof_seen[patch_dof];

              // Add the mapping: {local_cell_dof_index, patch_dof_index,
              // overlap_flag}
              cell_to_patch[c].emplace_back(i, patch_dof, is_overlap);

              // Mark this patch_dof as seen.
              patch_dof_seen[patch_dof] = true;
            }
          n_cell_to_patch[c] = cell_to_patch[c].size();
        }

      // Sort the tuples: non-overlapping first, then by cell_dof index.
      std::sort(cell_to_patch[c].begin(),
                cell_to_patch[c].end(),
                [](const auto &a, const auto &b) {
                  const bool overlap_a = std::get<2>(a);
                  const bool overlap_b = std::get<2>(b);

                  // Primary sort key: non-overlapping first (!overlap_a means
                  // true if a is non-overlapping)
                  if (overlap_a != overlap_b)
                    {
                      return !overlap_a;
                    }

                  // Secondary sort key: local cell DoF index
                  return std::get<0>(a) < std::get<0>(b);
                });
    }
  // Count the number of non-overlapping DoFs per cell
  for (IndexType c = 0; c < n_cells; ++c)
    {
      n_non_overlapping[c] = 0;
      for (const auto &mapping : cell_to_patch[c])
        {
          if (!std::get<2>(mapping)) // If not an overlap
            {
              n_non_overlapping[c]++;
            }
        }
    }

  // Print skipped DoFs for verification (optional)
  for (IndexType c = 0; c < n_cells; ++c)
    {
      output << "Skipped DoF local indices for cell " << c << ": ";
      for (const auto &index : skipped_dofs[c])
        output << index << " ";
      output << std::endl;
    }

  // Print cell_to_patch for verification
  for (IndexType c = 0; c < n_cells; ++c)
    {
      output << "Cell " << c
             << " to Patch mapping (cell_dof_idx, patch_dof_idx): ";
      for (const auto &pair : cell_to_patch[c])
        {
          output << "(" << std::get<0>(pair) << ", " << std::get<1>(pair)
                 << " , " << std::get<2>(pair) << ") ";
        }
      output << std::endl;
      output << "Number of non-overlapping DoFs : " << n_non_overlapping[c]
             << std::endl;
    }


  output << "\n================================\n" << std::endl;

  const auto index_string = "Index";
  const auto pair_string  = "Pair";


  // Print header:
  output << "template <>\n"
            "struct CellPatchLookup<"
         << dim << ", " << fe_degre << ">\n {\n";

  output << "using " << index_string << " = unsigned int;" << std::endl;
  output << "  using " << pair_string << " = std::pair<" << index_string << ", "
         << index_string << ">;\n\n";

  output << " static constexpr  " << index_string << "  n_cells = " << n_cells
         << " ;" << std::endl;

  output << " static constexpr  " << index_string
         << "  n_patch_dofs = " << n_interior_dofs << " ;" << std::endl;


  IndexType max_skipped_dofs = 0;
  if (!n_skipped_dofs.empty())
    max_skipped_dofs =
      *std::max_element(n_skipped_dofs.begin(), n_skipped_dofs.end());

  output << "static constexpr std::array<" << index_string << ", " << n_cells
         << "> n_skipped_dofs = {{";
  for (IndexType c = 0; c < n_cells; ++c)
    {
      output << n_skipped_dofs[c];
      if (c < n_cells - 1)
        output << ", ";
    }
  output << "}}; \n" << std::endl;


  output << "static constexpr std::array<std::array< " << index_string << " , "
         << max_skipped_dofs << ">, " << n_cells << "> skipped_dofs = {{"
         << std::endl;
  for (IndexType c = 0; c < n_cells; ++c)
    {
      output << "{{";
      for (size_t i = 0; i < max_skipped_dofs; ++i)
        {
          if (i >= n_skipped_dofs[c])
            output << numbers::invalid_unsigned_int;
          else
            output << skipped_dofs[c][i];
          if (i < max_skipped_dofs - 1)
            output << ", ";
        }
      output << "}}";
      if (c < n_cells - 1)
        output << ",";
      output << std::endl;
    }
  output << "}};" << std::endl;

  output << "\n " << std::endl;


  IndexType max_cell_to_patch =
    *std::max_element(n_cell_to_patch.begin(), n_cell_to_patch.end());

  output << "static constexpr std::array< " << index_string << ", " << n_cells
         << "> n_cell_to_patch = {{";
  for (IndexType c = 0; c < n_cells; ++c)
    {
      output << n_cell_to_patch[c];
      if (c < n_cells - 1)
        output << ", ";
    }
  output << "}}; \n" << std::endl;

  output << "static constexpr std::array<std::array<" << pair_string << ", "
         << max_cell_to_patch << ">, " << n_cells << "> cell_to_patch = {{"
         << std::endl;
  for (IndexType c = 0; c < n_cells; ++c)
    {
      output << "{{";
      for (size_t i = 0; i < max_cell_to_patch; ++i)
        {
          if (i < cell_to_patch[c].size())
            output << pair_string << "{" << std::get<0>(cell_to_patch[c][i])
                   << ", " << std::get<1>(cell_to_patch[c][i]) << "}";
          else
            output << "{" << numbers::invalid_unsigned_int << ", "
                   << numbers::invalid_unsigned_int << "}";
          if (i < max_cell_to_patch - 1)
            output << ", ";
        }
      output << "}}";
      if (c < n_cells - 1)
        output << ",";
      output << std::endl;
    }
  output << "}};" << std::endl;

  output << std::endl;

  output << "static constexpr std::array<" << index_string << ", " << n_cells
         << "> n_non_overlapping = {{";

  for (IndexType c = 0; c < n_cells; ++c)
    {
      output << n_non_overlapping[c];
      if (c < n_cells - 1)
        output << ", ";
    }
  output << "}};\n" << std::endl;

  output << "};" << std::endl;
  output << std::endl;
  output << std::endl;
  output << "=========================================" << std::endl;
  output << "=============END GENERATION==============" << std::endl;
  output << "=========================================" << std::endl;
}

int
main(int argc, char *argv[])
{
  const int fe_degre = argc > 1 ? std::stoi(argv[1]) : 2;
  initlog();


  generate<2>(fe_degre);

  generate<3>(fe_degre);
  return 0;
}