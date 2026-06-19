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

// Test loading various VTK dataset types (polydata, structured grid,
// rectilinear grid, structured points) from files, converting them to
// deal.II triangulations, and logging cell/face/vertex counts.

#include <deal.II/base/point.h>

#include <deal.II/grid/tria.h>

#include <deal.II/vtk/utilities.h>

#include <vtkPolyData.h>
#include <vtkRectilinearGrid.h>
#include <vtkSmartPointer.h>
#include <vtkStructuredGrid.h>
#include <vtkStructuredPoints.h>

#include <array>
#include <cstdio>
#include <string>

#include "../tests.h"


void
test_polydata()
{
  // Build a polydata with 2 triangles and 1 line
  auto polydata = vtkSmartPointer<vtkPolyData>::New();

  auto points = vtkSmartPointer<vtkPoints>::New();
  points->InsertNextPoint(0.0, 0.0, 0.0);
  points->InsertNextPoint(1.0, 0.0, 0.0);
  points->InsertNextPoint(0.5, 1.0, 0.0);
  points->InsertNextPoint(1.5, 0.0, 0.0);
  polydata->SetPoints(points);

  auto tri = vtkSmartPointer<vtkTriangle>::New();
  tri->GetPointIds()->SetId(0, 0);
  tri->GetPointIds()->SetId(1, 1);
  tri->GetPointIds()->SetId(2, 2);
  polydata->InsertNextCell(tri->GetCellType(), tri->GetPointIds());

  auto tri2 = vtkSmartPointer<vtkTriangle>::New();
  tri2->GetPointIds()->SetId(0, 1);
  tri2->GetPointIds()->SetId(1, 2);
  tri2->GetPointIds()->SetId(2, 3);
  polydata->InsertNextCell(tri2->GetCellType(), tri2->GetPointIds());

  auto line = vtkSmartPointer<vtkLine>::New();
  line->GetPointIds()->SetId(0, 0);
  line->GetPointIds()->SetId(1, 3);
  polydata->InsertNextCell(line->GetCellType(), line->GetPointIds());

  const std::string path = "load_vtk_dataset_roundtrip_01_polydata.vtk";
  {
    auto writer = vtkSmartPointer<vtkUnstructuredGridWriter>::New();
    writer->SetFileName(path.c_str());
    auto grid = VTKWrappers::internal::convert_to_unstructured_grid(polydata);
    writer->SetInputData(grid);
    writer->Write();
  }

  deallog << "Polydata written" << std::endl;

  auto loaded = VTKWrappers::internal::load_vtk_file(path, false);
  deallog << "Polydata loaded: " << loaded->GetNumberOfPoints() << " points, "
          << loaded->GetNumberOfCells() << " cells" << std::endl;

  Triangulation<2, 2> tria;
  VTKWrappers::unstructured_grid_to_dealii_triangulation(*loaded, tria);
  deallog << "Polydata -> deal.II: " << tria.n_active_cells() << " cells, "
          << tria.n_vertices() << " vertices" << std::endl;

  std::remove(path.c_str());
}


void
test_structured_grid()
{
  // Build a 3x2 structured grid (3 hex cells)
  auto sg = vtkSmartPointer<vtkStructuredGrid>::New();

  auto points = vtkSmartPointer<vtkPoints>::New();
  for (int j = 0; j < 2; ++j)
    for (int i = 0; i < 3; ++i)
      points->InsertNextPoint(static_cast<double>(i),
                              static_cast<double>(j),
                              0.0);
  sg->SetPoints(points);

  sg->SetDimensions(3, 2, 1);
  sg->SetNumberOfCells(3, 1, 1); // 3x1 = 3 hex cells

  const std::string path = "load_vtk_dataset_roundtrip_01_structured_grid.vtk";
  {
    auto writer = vtkSmartPointer<vtkUnstructuredGridWriter>::New();
    writer->SetFileName(path.c_str());
    auto grid = VTKWrappers::internal::convert_to_unstructured_grid(sg);
    writer->SetInputData(grid);
    writer->Write();
  }

  deallog << "StructuredGrid written" << std::endl;

  auto loaded = VTKWrappers::internal::load_vtk_file(path, false);
  deallog << "StructuredGrid loaded: " << loaded->GetNumberOfPoints()
          << " points, " << loaded->GetNumberOfCells() << " cells" << std::endl;

  Triangulation<2, 2> tria;
  VTKWrappers::unstructured_grid_to_dealii_triangulation(*loaded, tria);
  deallog << "StructuredGrid -> deal.II: " << tria.n_active_cells()
          << " cells, " << tria.n_vertices() << " vertices" << std::endl;

  std::remove(path.c_str());
}


void
test_rectilinear_grid()
{
  // Build a 2x2x2 rectilinear grid
  auto rg = vtkSmartPointer<vtkRectilinearGrid>::New();

  auto x = vtkSmartPointer<vtkDoubleArray>::New();
  x->SetName("x");
  x->InsertValue(0, 0.0);
  x->InsertValue(1, 1.0);
  rg->SetXCoordinates(x);

  auto y = vtkSmartPointer<vtkDoubleArray>::New();
  y->SetName("y");
  y->InsertValue(0, 0.0);
  y->InsertValue(1, 1.0);
  rg->SetYCoordinates(y);

  auto z = vtkSmartPointer<vtkDoubleArray>::New();
  z->SetName("z");
  z->InsertValue(0, 0.0);
  z->InsertValue(1, 1.0);
  rg->SetZCoordinates(z);

  rg->SetDimensions(2, 2, 2); // 1x1x1 = 1 hex cell

  const std::string path = "load_vtk_dataset_roundtrip_01_rectilinear_grid.vtk";
  {
    auto writer = vtkSmartPointer<vtkUnstructuredGridWriter>::New();
    writer->SetFileName(path.c_str());
    auto grid = VTKWrappers::internal::convert_to_unstructured_grid(rg);
    writer->SetInputData(grid);
    writer->Write();
  }

  deallog << "RectilinearGrid written" << std::endl;

  auto loaded = VTKWrappers::internal::load_vtk_file(path, false);
  deallog << "RectilinearGrid loaded: " << loaded->GetNumberOfPoints()
          << " points, " << loaded->GetNumberOfCells() << " cells" << std::endl;

  Triangulation<3, 3> tria;
  VTKWrappers::unstructured_grid_to_dealii_triangulation(*loaded, tria);
  deallog << "RectilinearGrid -> deal.II: " << tria.n_active_cells()
          << " cells, " << tria.n_vertices() << " vertices" << std::endl;

  std::remove(path.c_str());
}


void
test_structured_points()
{
  // Build a 3x3 structured points (2x2 = 4 hex cells)
  auto sp = vtkSmartPointer<vtkStructuredPoints>::New();
  sp->SetDimensions(3, 3, 1); // 2x2x1 = 4 hex cells
  sp->SetSpacing(1.0, 1.0, 1.0);
  sp->SetOrigin(0.0, 0.0, 0.0);

  const std::string path =
    "load_vtk_dataset_roundtrip_01_structured_points.vtk";
  {
    auto writer = vtkSmartPointer<vtkUnstructuredGridWriter>::New();
    writer->SetFileName(path.c_str());
    auto grid = VTKWrappers::internal::convert_to_unstructured_grid(sp);
    writer->SetInputData(grid);
    writer->Write();
  }

  deallog << "StructuredPoints written" << std::endl;

  auto loaded = VTKWrappers::internal::load_vtk_file(path, false);
  deallog << "StructuredPoints loaded: " << loaded->GetNumberOfPoints()
          << " points, " << loaded->GetNumberOfCells() << " cells" << std::endl;

  Triangulation<2, 2> tria;
  VTKWrappers::unstructured_grid_to_dealii_triangulation(*loaded, tria);
  deallog << "StructuredPoints -> deal.II: " << tria.n_active_cells()
          << " cells, " << tria.n_vertices() << " vertices" << std::endl;

  std::remove(path.c_str());
}


int
main()
{
  initlog();

  test_polydata();
  test_structured_grid();
  test_rectilinear_grid();
  test_structured_points();

  return 0;
}
