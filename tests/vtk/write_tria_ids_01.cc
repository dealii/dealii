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

// Test writing VTK files with material, boundary, and manifold ids through
// dealii_triangulation_to_unstructured_grid() + write_vtk(), then reading them
// back through read_tria().

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>

#include <deal.II/vtk/utilities.h>

#include <map>

#include "../tests.h"

using namespace dealii;

namespace
{
  template <int dim, int spacedim>
  void
  test_file(const std::string &filename)
  {
    Triangulation<dim, spacedim> tria_in;
    GridGenerator::subdivided_hyper_cube(tria_in, 1);

    auto cell = tria_in.begin_active();
    if constexpr (dim == 2)
      cell->set_material_id(7);
    else if constexpr (dim == 3)
      cell->set_material_id(11);

    for (unsigned int f = 0; f < cell->n_faces(); ++f)
      if (cell->face(f)->at_boundary())
        {
          if constexpr (dim == 2)
            cell->face(f)->set_boundary_id(13);
          else if constexpr (dim == 3)
            cell->face(f)->set_manifold_id(5);
        }

    VTKWrappers::write_vtk(
      filename, tria_in, "material", "boundary", "manifold");

    Triangulation<dim, spacedim> tria_out;
    VTKWrappers::read_tria(
      filename, tria_out, true, 0.0, "material", "boundary", "manifold");

    deallog << "Reading file: " << filename << std::endl;
    deallog << "  Triangulation: " << tria_out.n_active_cells() << " cells, "
            << tria_out.n_vertices() << " vertices" << std::endl;

    std::map<int, unsigned int> material_id_counts;
    std::map<int, unsigned int> cell_manifold_id_counts;
    std::map<int, unsigned int> boundary_id_counts;
    std::map<int, unsigned int> boundary_manifold_id_counts;

    for (const auto &out_cell : tria_out.active_cell_iterators())
      {
        ++material_id_counts[static_cast<int>(out_cell->material_id())];
        ++cell_manifold_id_counts[static_cast<int>(out_cell->manifold_id())];

        for (unsigned int f = 0; f < out_cell->n_faces(); ++f)
          if (out_cell->face(f)->at_boundary())
            {
              ++boundary_id_counts[static_cast<int>(
                out_cell->face(f)->boundary_id())];
              ++boundary_manifold_id_counts[static_cast<int>(
                out_cell->face(f)->manifold_id())];
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

    std::remove(filename.c_str());
  }
} // namespace



int
main()
{
  initlog();

  test_file<2, 2>("write_tria_ids_01_2d.vtk");
  test_file<2, 2>("write_tria_ids_01_2d.vtu");
  test_file<3, 3>("write_tria_ids_01_3d.vtk");
  test_file<3, 3>("write_tria_ids_01_3d.vtu");

  return 0;
}
