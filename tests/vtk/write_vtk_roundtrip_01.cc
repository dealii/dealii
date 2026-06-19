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

// Test write_vtk and load_vtk_file round-trip: build a VTK dataset,
// write it via write_vtk, read it back via load_vtk_file, convert both
// to unstructured grids, and compare cell/point counts.

#include <deal.II/grid/tria.h>

#include <deal.II/vtk/utilities.h>

#include <vtkCellData.h>
#include <vtkDoubleArray.h>
#include <vtkIntArray.h>
#include <vtkPointData.h>
#include <vtkPolyData.h>
#include <vtkSmartPointer.h>
#include <vtkStructuredGrid.h>
#include <vtkTriangle.h>
#include <vtkUnstructuredGrid.h>

#include <string>

#include "../tests.h"


void
test_polydata_roundtrip(const std::string &ext)
{
  // Build a polydata with 3 triangles
  auto polydata = vtkSmartPointer<vtkPolyData>::New();

  auto points = vtkSmartPointer<vtkPoints>::New();
  points->InsertNextPoint(0.0, 0.0, 0.0);
  points->InsertNextPoint(1.0, 0.0, 0.0);
  points->InsertNextPoint(0.5, 1.0, 0.0);
  points->InsertNextPoint(1.5, 0.5, 0.0);
  points->InsertNextPoint(2.0, 0.0, 0.0);
  polydata->SetPoints(points);
  polydata->Allocate(3);

  auto tri0 = vtkSmartPointer<vtkTriangle>::New();
  tri0->GetPointIds()->SetId(0, 0);
  tri0->GetPointIds()->SetId(1, 1);
  tri0->GetPointIds()->SetId(2, 2);
  polydata->InsertNextCell(tri0->GetCellType(), tri0->GetPointIds());

  auto tri1 = vtkSmartPointer<vtkTriangle>::New();
  tri1->GetPointIds()->SetId(0, 1);
  tri1->GetPointIds()->SetId(1, 2);
  tri1->GetPointIds()->SetId(2, 3);
  polydata->InsertNextCell(tri1->GetCellType(), tri1->GetPointIds());

  auto tri2 = vtkSmartPointer<vtkTriangle>::New();
  tri2->GetPointIds()->SetId(0, 1);
  tri2->GetPointIds()->SetId(1, 3);
  tri2->GetPointIds()->SetId(2, 4);
  polydata->InsertNextCell(tri2->GetCellType(), tri2->GetPointIds());

  // Add point data
  auto scalars = vtkSmartPointer<vtkDoubleArray>::New();
  scalars->SetName("temperature");
  scalars->SetNumberOfComponents(1);
  scalars->InsertValue(0, 100.0);
  scalars->InsertValue(1, 200.0);
  scalars->InsertValue(2, 150.0);
  scalars->InsertValue(3, 175.0);
  scalars->InsertValue(4, 125.0);
  polydata->GetPointData()->AddArray(scalars);

  const std::string path = "write_vtk_roundtrip_01_polydata" + ext;
  deallog << "Writing polydata" << std::endl;

  VTKWrappers::internal::write_vtk(path, polydata);

  // Read back
  auto loaded = VTKWrappers::internal::load_vtk_file(path, false);
  deallog << "  Loaded: " << loaded->GetNumberOfPoints() << " points, "
          << loaded->GetNumberOfCells() << " cells" << std::endl;

  // Convert original to unstructured grid
  auto original_grid =
    VTKWrappers::internal::convert_to_unstructured_grid(polydata);
  deallog << "  Original -> unstructured: "
          << original_grid->GetNumberOfPoints() << " points, "
          << original_grid->GetNumberOfCells() << " cells" << std::endl;

  // Compare
  deallog << "  Point count match: "
          << (loaded->GetNumberOfPoints() == original_grid->GetNumberOfPoints())
          << std::endl;
  deallog << "  Cell count match: "
          << (loaded->GetNumberOfCells() == original_grid->GetNumberOfCells())
          << std::endl;
  deallog << "  Converted point array 'temperature': "
          << original_grid->GetPointData()->HasArray("temperature")
          << std::endl;

  // Convert to deal.II and compare
  Triangulation<2, 2> tria_loaded;
  VTKWrappers::unstructured_grid_to_dealii_triangulation(*loaded, tria_loaded);

  Triangulation<2, 2> tria_original;
  VTKWrappers::unstructured_grid_to_dealii_triangulation(*original_grid,
                                                         tria_original);

  deallog << "  Triangulation cell count match: "
          << (tria_loaded.n_active_cells() == tria_original.n_active_cells())
          << std::endl;
  deallog << "  Triangulation vertex count match: "
          << (tria_loaded.n_vertices() == tria_original.n_vertices())
          << std::endl;

  deallog << "  Round-trip OK: " << tria_loaded.n_active_cells() << " cells, "
          << tria_loaded.n_vertices() << " vertices" << std::endl;

  std::remove(path.c_str());
}


void
test_structured_grid_roundtrip(const std::string &ext)
{
  // Build a 2x2 structured grid (4 hex cells)
  auto sg = vtkSmartPointer<vtkStructuredGrid>::New();

  auto points = vtkSmartPointer<vtkPoints>::New();
  for (int j = 0; j < 3; ++j)
    for (int i = 0; i < 3; ++i)
      points->InsertNextPoint(static_cast<double>(i),
                              static_cast<double>(j),
                              0.0);
  sg->SetPoints(points);

  sg->SetDimensions(3, 3, 1);

  // Add cell data
  auto cell_ids = vtkSmartPointer<vtkIntArray>::New();
  cell_ids->SetName("region");
  cell_ids->SetNumberOfComponents(1);
  for (int i = 0; i < 4; ++i)
    cell_ids->InsertValue(i, i % 2);
  sg->GetCellData()->AddArray(cell_ids);

  const std::string path = "write_vtk_roundtrip_01_structured_grid" + ext;
  deallog << "Writing StructuredGrid" << std::endl;

  VTKWrappers::internal::write_vtk(path, sg);

  auto loaded = VTKWrappers::internal::load_vtk_file(path, false);
  deallog << "  Loaded: " << loaded->GetNumberOfPoints() << " points, "
          << loaded->GetNumberOfCells() << " cells" << std::endl;

  auto original_grid = VTKWrappers::internal::convert_to_unstructured_grid(sg);
  deallog << "  Original -> unstructured: "
          << original_grid->GetNumberOfPoints() << " points, "
          << original_grid->GetNumberOfCells() << " cells" << std::endl;

  deallog << "  Point count match: "
          << (loaded->GetNumberOfPoints() == original_grid->GetNumberOfPoints())
          << std::endl;
  deallog << "  Cell count match: "
          << (loaded->GetNumberOfCells() == original_grid->GetNumberOfCells())
          << std::endl;
  deallog << "  Converted cell array 'region': "
          << original_grid->GetCellData()->HasArray("region") << std::endl;

  Triangulation<2, 2> tria_loaded;
  VTKWrappers::unstructured_grid_to_dealii_triangulation(*loaded, tria_loaded);
  Triangulation<2, 2> tria_original;
  VTKWrappers::unstructured_grid_to_dealii_triangulation(*original_grid,
                                                         tria_original);
  deallog << "  Triangulation cell count match: "
          << (tria_loaded.n_active_cells() == tria_original.n_active_cells())
          << std::endl;
  deallog << "  Triangulation vertex count match: "
          << (tria_loaded.n_vertices() == tria_original.n_vertices())
          << std::endl;

  deallog << "  Round-trip OK: " << tria_loaded.n_active_cells() << " cells, "
          << tria_loaded.n_vertices() << " vertices" << std::endl;

  std::remove(path.c_str());
}


void
test_unstructured_grid_roundtrip(const std::string &ext)
{
  // Build an unstructured grid directly
  auto ugrid = vtkSmartPointer<vtkUnstructuredGrid>::New();

  auto points = vtkSmartPointer<vtkPoints>::New();
  points->InsertNextPoint(0.0, 0.0, 0.0);
  points->InsertNextPoint(1.0, 0.0, 0.0);
  points->InsertNextPoint(1.0, 1.0, 0.0);
  points->InsertNextPoint(0.0, 1.0, 0.0);
  points->InsertNextPoint(0.5, 0.5, 1.0);
  ugrid->SetPoints(points);

  // Insert a tetrahedron
  vtkIdType tet_ids[4] = {0, 1, 2, 4};
  ugrid->InsertNextCell(VTK_TETRA, 4, tet_ids);

  // Insert a pyramid
  vtkIdType pyr_ids[5] = {0, 1, 2, 3, 4};
  ugrid->InsertNextCell(VTK_PYRAMID, 5, pyr_ids);

  const std::string path = "write_vtk_roundtrip_01_unstructured_grid" + ext;
  deallog << "Writing UnstructuredGrid" << std::endl;

  VTKWrappers::internal::write_vtk(path, ugrid);

  auto loaded = VTKWrappers::internal::load_vtk_file(path, false);
  deallog << "  Loaded: " << loaded->GetNumberOfPoints() << " points, "
          << loaded->GetNumberOfCells() << " cells" << std::endl;

  auto original_grid =
    VTKWrappers::internal::convert_to_unstructured_grid(ugrid);
  deallog << "  Original -> unstructured: "
          << original_grid->GetNumberOfPoints() << " points, "
          << original_grid->GetNumberOfCells() << " cells" << std::endl;

  deallog << "  Point count match: "
          << (loaded->GetNumberOfPoints() == original_grid->GetNumberOfPoints())
          << std::endl;
  deallog << "  Cell count match: "
          << (loaded->GetNumberOfCells() == original_grid->GetNumberOfCells())
          << std::endl;

  Triangulation<3, 3> tria_loaded;
  VTKWrappers::unstructured_grid_to_dealii_triangulation(*loaded, tria_loaded);
  Triangulation<3, 3> tria_original;
  VTKWrappers::unstructured_grid_to_dealii_triangulation(*original_grid,
                                                         tria_original);
  deallog << "  Triangulation cell count match: "
          << (tria_loaded.n_active_cells() == tria_original.n_active_cells())
          << std::endl;
  deallog << "  Triangulation vertex count match: "
          << (tria_loaded.n_vertices() == tria_original.n_vertices())
          << std::endl;

  deallog << "  Round-trip OK: " << tria_loaded.n_active_cells() << " cells, "
          << tria_loaded.n_vertices() << " vertices" << std::endl;

  std::remove(path.c_str());
}


int
main()
{
  initlog();

  deallog << "=== Legacy .vtk ===" << std::endl;
  test_polydata_roundtrip(".vtk");
  test_structured_grid_roundtrip(".vtk");
  test_unstructured_grid_roundtrip(".vtk");

  return 0;
}
