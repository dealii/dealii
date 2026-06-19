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

// Test conversion of VTK cell data arrays to material, boundary, and manifold
// ids while converting a vtkUnstructuredGrid to a deal.II Triangulation.

#include <deal.II/grid/tria.h>

#include <deal.II/vtk/utilities.h>

#include <vtkCellData.h>
#include <vtkCellType.h>
#include <vtkIntArray.h>
#include <vtkPoints.h>
#include <vtkSmartPointer.h>
#include <vtkUnstructuredGrid.h>

#include <algorithm>
#include <array>

#include "../tests.h"


vtkSmartPointer<vtkIntArray>
make_cell_array(const std::string &name, const std::vector<int> &values)
{
  auto array = vtkSmartPointer<vtkIntArray>::New();
  array->SetName(name.c_str());
  array->SetNumberOfComponents(1);
  array->SetNumberOfTuples(values.size());
  for (vtkIdType i = 0; i < static_cast<vtkIdType>(values.size()); ++i)
    array->SetTypedComponent(i, 0, values[i]);

  return array;
}


template <int dim>
bool
has_vertices(const TriaIterator<TriaAccessor<dim - 1, dim, dim>> &face,
             const std::vector<unsigned int>                     &vertices)
{
  std::vector<unsigned int> face_vertices(face->n_vertices());
  for (unsigned int v = 0; v < face->n_vertices(); ++v)
    face_vertices[v] = face->vertex_index(v);

  std::sort(face_vertices.begin(), face_vertices.end());

  std::vector<unsigned int> sorted_vertices = vertices;
  std::sort(sorted_vertices.begin(), sorted_vertices.end());

  return face_vertices == sorted_vertices;
}


vtkSmartPointer<vtkPoints>
make_points(const std::vector<std::array<double, 3>> &point_locations)
{
  auto points = vtkSmartPointer<vtkPoints>::New();
  points->SetNumberOfPoints(point_locations.size());
  for (vtkIdType i = 0; i < static_cast<vtkIdType>(point_locations.size()); ++i)
    points->SetPoint(i, point_locations[i].data());

  return points;
}


void
test_2d_cell(const std::string                        &name,
             const int                                 vtk_cell_type,
             const std::vector<std::array<double, 3>> &point_locations,
             const std::vector<vtkIdType>             &cell_vertices)
{
  const std::vector<vtkIdType> line_vertices = {0, 1};

  auto grid = vtkSmartPointer<vtkUnstructuredGrid>::New();
  grid->SetPoints(make_points(point_locations));

  grid->InsertNextCell(vtk_cell_type,
                       cell_vertices.size(),
                       cell_vertices.data());

  grid->InsertNextCell(VTK_LINE, line_vertices.size(), line_vertices.data());

  grid->GetCellData()->AddArray(make_cell_array("material", {7, 99}));
  grid->GetCellData()->AddArray(make_cell_array("boundary", {9, 42}));
  grid->GetCellData()->AddArray(make_cell_array("manifold", {17, 142}));

  Triangulation<2> triangulation;
  VTKWrappers::unstructured_grid_to_dealii_triangulation(
    *grid, triangulation, "material", "boundary", "manifold");

  const auto cell = triangulation.begin_active();
  deallog << "2d " << name << " cells: " << triangulation.n_active_cells()
          << std::endl;
  deallog << "2d material id: " << cell->material_id() << std::endl;
  deallog << "2d cell manifold id: " << cell->manifold_id() << std::endl;

  for (unsigned int f = 0; f < cell->n_faces(); ++f)
    if (has_vertices<2>(cell->face(f),
                        std::vector<unsigned int>(line_vertices.begin(),
                                                  line_vertices.end())))
      deallog << "2d boundary line ids: " << cell->face(f)->boundary_id() << ' '
              << cell->face(f)->manifold_id() << std::endl;

  Triangulation<2> default_triangulation;
  VTKWrappers::unstructured_grid_to_dealii_triangulation(*grid,
                                                         default_triangulation);
  deallog << "2d " << name
          << " default cells: " << default_triangulation.n_active_cells()
          << std::endl;
  deallog << "2d default material id: "
          << default_triangulation.begin_active()->material_id() << std::endl;
}


void
test_2d()
{
  test_2d_cell("quad",
               VTK_QUAD,
               {{{0.0, 0.0, 0.0}},
                {{1.0, 0.0, 0.0}},
                {{1.0, 1.0, 0.0}},
                {{0.0, 1.0, 0.0}}},
               {0, 1, 2, 3});

  test_2d_cell("triangle",
               VTK_TRIANGLE,
               {{{0.0, 0.0, 0.0}}, {{1.0, 0.0, 0.0}}, {{0.0, 1.0, 0.0}}},
               {0, 1, 2});
}


void
test_3d_cell(const std::string                        &name,
             const int                                 vtk_cell_type,
             const std::vector<std::array<double, 3>> &point_locations,
             const std::vector<vtkIdType>             &cell_vertices,
             const int                                 vtk_face_type,
             const std::vector<vtkIdType>             &face_vertices)
{
  auto grid = vtkSmartPointer<vtkUnstructuredGrid>::New();
  grid->SetPoints(make_points(point_locations));

  grid->InsertNextCell(vtk_cell_type,
                       cell_vertices.size(),
                       cell_vertices.data());

  grid->InsertNextCell(vtk_face_type,
                       face_vertices.size(),
                       face_vertices.data());

  grid->GetCellData()->AddArray(make_cell_array("material", {11, 88}));
  grid->GetCellData()->AddArray(make_cell_array("boundary", {5, 55}));
  grid->GetCellData()->AddArray(make_cell_array("manifold", {21, 155}));

  Triangulation<3> triangulation;
  VTKWrappers::unstructured_grid_to_dealii_triangulation(
    *grid, triangulation, "material", "boundary", "manifold");

  const auto cell = triangulation.begin_active();
  deallog << "3d " << name << " cells: " << triangulation.n_active_cells()
          << std::endl;
  deallog << "3d material id: " << cell->material_id() << std::endl;
  deallog << "3d cell manifold id: " << cell->manifold_id() << std::endl;

  for (unsigned int f = 0; f < cell->n_faces(); ++f)
    if (has_vertices<3>(cell->face(f),
                        std::vector<unsigned int>(face_vertices.begin(),
                                                  face_vertices.end())))
      deallog << "3d boundary face ids: " << cell->face(f)->boundary_id() << ' '
              << cell->face(f)->manifold_id() << std::endl;

  Triangulation<3> default_triangulation;
  VTKWrappers::unstructured_grid_to_dealii_triangulation(*grid,
                                                         default_triangulation);
  deallog << "3d " << name
          << " default cells: " << default_triangulation.n_active_cells()
          << std::endl;
  deallog << "3d default material id: "
          << default_triangulation.begin_active()->material_id() << std::endl;
}


void
test_3d()
{
  test_3d_cell("hex",
               VTK_HEXAHEDRON,
               {{{0.0, 0.0, 0.0}},
                {{1.0, 0.0, 0.0}},
                {{1.0, 1.0, 0.0}},
                {{0.0, 1.0, 0.0}},
                {{0.0, 0.0, 1.0}},
                {{1.0, 0.0, 1.0}},
                {{1.0, 1.0, 1.0}},
                {{0.0, 1.0, 1.0}}},
               {0, 1, 2, 3, 4, 5, 6, 7},
               VTK_QUAD,
               {0, 1, 2, 3});

  test_3d_cell("tet",
               VTK_TETRA,
               {{{0.0, 0.0, 0.0}},
                {{1.0, 0.0, 0.0}},
                {{0.0, 1.0, 0.0}},
                {{0.0, 0.0, 1.0}}},
               {0, 1, 2, 3},
               VTK_TRIANGLE,
               {0, 1, 2});

  test_3d_cell("pyramid",
               VTK_PYRAMID,
               {{{0.0, 0.0, 0.0}},
                {{1.0, 0.0, 0.0}},
                {{1.0, 1.0, 0.0}},
                {{0.0, 1.0, 0.0}},
                {{0.5, 0.5, 1.0}}},
               {0, 1, 2, 3, 4},
               VTK_QUAD,
               {0, 1, 2, 3});

  test_3d_cell("wedge",
               VTK_WEDGE,
               {{{0.0, 0.0, 0.0}},
                {{1.0, 0.0, 0.0}},
                {{0.0, 1.0, 0.0}},
                {{0.0, 0.0, 1.0}},
                {{1.0, 0.0, 1.0}},
                {{0.0, 1.0, 1.0}}},
               {0, 1, 2, 3, 4, 5},
               VTK_TRIANGLE,
               {0, 1, 2});
}


int
main()
{
  initlog();

  test_2d();
  test_3d();
}
