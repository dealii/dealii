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

// Test reading VTK files with boundary-id cell data and mapping the
// "Boundary ID" field to deal.II boundary ids.

#include <deal.II/grid/tria.h>

#include <deal.II/vtk/utilities.h>

#include <vtkCellData.h>
#include <vtkDataArray.h>

#include <map>

#include "../tests.h"

using namespace dealii;

namespace
{
  std::string
  find_boundary_field_name(const vtkUnstructuredGrid &grid)
  {
    auto        *grid_ptr  = const_cast<vtkUnstructuredGrid *>(&grid);
    vtkCellData *cell_data = grid_ptr->GetCellData();
    if (cell_data == nullptr)
      return "";

    for (int i = 0; i < cell_data->GetNumberOfArrays(); ++i)
      {
        vtkDataArray *array = cell_data->GetArray(i);
        if (array == nullptr || array->GetName() == nullptr)
          continue;

        const std::string name = array->GetName();
        if (name == "Boundary ID" || name == "Boundary%20ID")
          return name;
      }

    return "";
  }



  template <int dim, int spacedim>
  void
  test_file(const std::string &filename)
  {
    const auto full_name = std::string(SOURCE_DIR) + "/" + filename;
    auto grid = VTKWrappers::internal::load_vtk_file(full_name, false, 0.0);

    deallog << "Reading file: " << filename << std::endl;
    deallog << "  Unstructured grid: " << grid->GetNumberOfPoints()
            << " points, " << grid->GetNumberOfCells() << " cells" << std::endl;

    vtkCellData *cell_data = grid->GetCellData();
    deallog << "  Cell arrays:";
    if (cell_data != nullptr)
      for (int i = 0; i < cell_data->GetNumberOfArrays(); ++i)
        {
          vtkDataArray *array = cell_data->GetArray(i);
          if (array != nullptr && array->GetName() != nullptr)
            deallog << " " << array->GetName();
        }
    deallog << std::endl;

    const std::string boundary_field = find_boundary_field_name(*grid);
    deallog << "  Boundary field: " << boundary_field << std::endl;

    Triangulation<dim, spacedim> triangulation;
    VTKWrappers::unstructured_grid_to_dealii_triangulation(*grid,
                                                           triangulation,
                                                           "",
                                                           boundary_field);

    deallog << "  Triangulation: " << triangulation.n_active_cells()
            << " cells, " << triangulation.n_vertices() << " vertices"
            << std::endl;

    std::map<int, unsigned int> boundary_id_counts;
    unsigned int                n_boundary_faces = 0;
    for (const auto &cell : triangulation.active_cell_iterators())
      for (unsigned int f = 0; f < cell->n_faces(); ++f)
        if (cell->face(f)->at_boundary())
          {
            ++n_boundary_faces;
            ++boundary_id_counts[static_cast<int>(
              cell->face(f)->boundary_id())];
          }

    deallog << "  Boundary faces: " << n_boundary_faces << std::endl;
    for (const auto &[boundary_id, count] : boundary_id_counts)
      deallog << "    id " << boundary_id << ": " << count << std::endl;
  }
} // namespace



int
main()
{
  initlog();

  test_file<2, 2>("data/square_with_ids.vtk");
  test_file<2, 2>("data/square_with_ids.vtu");
  test_file<3, 3>("data/cube_with_ids.vtk");
  test_file<3, 3>("data/cube_with_ids.vtu");

  return 0;
}
