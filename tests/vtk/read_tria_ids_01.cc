// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
// Copyright (C) 2026 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Detailed license information governing the source code and contributions
// can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
//
// -----------------------------------------------------------------------------

// Test reading VTK files with material, boundary, and manifold ids through
// the public read_tria() interface.

#include <deal.II/grid/tria.h>

#include <deal.II/vtk/utilities.h>

#include <map>

#include "../tests.h"

using namespace dealii;

namespace
{
  template <int dim, int spacedim>
  void
  test_file(const std::string &filename, const std::string &manifold_field)
  {
    const auto full_name = std::string(SOURCE_DIR) + "/" + filename;

    Triangulation<dim, spacedim> triangulation;
    VTKWrappers::read_tria(full_name,
                           triangulation,
                           true,
                           0.0,
                           "MaterialID",
                           "Boundary ID",
                           manifold_field);

    deallog << "Reading file: " << filename << std::endl;
    deallog << "  Triangulation: " << triangulation.n_active_cells()
            << " cells, " << triangulation.n_vertices() << " vertices"
            << std::endl;

    std::map<int, unsigned int> material_id_counts;
    std::map<int, unsigned int> cell_manifold_id_counts;
    std::map<int, unsigned int> boundary_id_counts;
    std::map<int, unsigned int> boundary_manifold_id_counts;

    for (const auto &cell : triangulation.active_cell_iterators())
      {
        ++material_id_counts[static_cast<int>(cell->material_id())];
        ++cell_manifold_id_counts[static_cast<int>(cell->manifold_id())];

        for (unsigned int f = 0; f < cell->n_faces(); ++f)
          if (cell->face(f)->at_boundary())
            {
              ++boundary_id_counts[static_cast<int>(
                cell->face(f)->boundary_id())];
              ++boundary_manifold_id_counts[static_cast<int>(
                cell->face(f)->manifold_id())];
            }
      }

    deallog << "  Material ids:" << std::endl;
    for (const auto &[id, count] : material_id_counts)
      deallog << "    id " << id << ": " << count << std::endl;

    deallog << "  Cell manifold ids:" << std::endl;
    for (const auto &[id, count] : cell_manifold_id_counts)
      deallog << "    id " << id << ": " << count << std::endl;

    deallog << "  Boundary ids:" << std::endl;
    for (const auto &[id, count] : boundary_id_counts)
      deallog << "    id " << id << ": " << count << std::endl;

    deallog << "  Boundary manifold ids:" << std::endl;
    for (const auto &[id, count] : boundary_manifold_id_counts)
      deallog << "    id " << id << ": " << count << std::endl;
  }
} // namespace



int
main()
{
  initlog();

  test_file<2, 2>("data/square_with_ids.vtk", "ManifoldID");
  test_file<2, 2>("data/square_with_ids.vtu", "ManifoldID");
  test_file<3, 3>("data/cube_with_ids.vtk", "");
  test_file<3, 3>("data/cube_with_ids.vtu", "");

  return 0;
}
