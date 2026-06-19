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


// Test: export material, boundary and manifold ids from a deal.II
// Triangulation to a VTK UnstructuredGrid. Verifies that cell-data arrays
// and appended subcells (boundary faces) contain the expected id values.
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>

#include <deal.II/vtk/utilities.h>

#include <vtkCellData.h>
#include <vtkIntArray.h>
#include <vtkPoints.h>
#include <vtkSmartPointer.h>
#include <vtkUnstructuredGrid.h>

#include "../tests.h"

using namespace dealii;

template <int dim, int spacedim>
void
test()
{
  Triangulation<dim, spacedim> tria;
  GridGenerator::subdivided_hyper_cube(tria, 1);

  // set material id on the single cell
  auto cell = tria.begin_active();
  if constexpr (dim == 2)
    cell->set_material_id(7);
  else if constexpr (dim == 3)
    cell->set_material_id(11);

  if constexpr (dim == 2)
    for (unsigned int f = 0; f < cell->n_faces(); ++f)
      if (cell->face(f)->at_boundary())
        cell->face(f)->set_boundary_id(13);

  if constexpr (dim == 3)
    for (unsigned int f = 0; f < cell->n_faces(); ++f)
      if (cell->face(f)->at_boundary())
        cell->face(f)->set_manifold_id(5);

  auto grid = VTKWrappers::dealii_triangulation_to_unstructured_grid(
    tria, "material", "boundary", "manifold");

  auto material =
    vtkIntArray::SafeDownCast(grid->GetCellData()->GetArray("material"));
  auto boundary =
    vtkIntArray::SafeDownCast(grid->GetCellData()->GetArray("boundary"));
  auto manifold =
    vtkIntArray::SafeDownCast(grid->GetCellData()->GetArray("manifold"));

  const bool has_material = (material != nullptr);
  const bool has_boundary = (boundary != nullptr);
  const bool has_manifold = (manifold != nullptr);

  deallog << dim << "d arrays present: material=" << has_material
          << " boundary=" << has_boundary << " manifold=" << has_manifold
          << std::endl;

  if (has_material)
    deallog << dim << "d material id (cell 0): " << material->GetValue(0)
            << std::endl;
  if (has_boundary)
    deallog << dim << "d boundary id (cell 0): " << boundary->GetValue(0)
            << std::endl;
  if (has_manifold)
    deallog << dim << "d manifold id (cell 0): " << manifold->GetValue(0)
            << std::endl;

  // Determine number of tuples from an available array
  vtkIdType n_tuples = 0;
  if (has_material)
    n_tuples = material->GetNumberOfTuples();
  else if (has_boundary)
    n_tuples = boundary->GetNumberOfTuples();
  else if (has_manifold)
    n_tuples = manifold->GetNumberOfTuples();

  deallog << dim << "d n_tuples: " << n_tuples << std::endl;

  bool found = false;
  if constexpr (dim == 2)
    {
      if (has_boundary)
        for (vtkIdType i = 1; i < n_tuples; ++i)
          if (boundary->GetValue(i) == 13)
            found = true;
      deallog << "2d found boundary face with id 13: " << found << std::endl;
    }
  else if constexpr (dim == 3)
    {
      if (has_manifold)
        for (vtkIdType i = 1; i < n_tuples; ++i)
          if (manifold->GetValue(i) == 5)
            found = true;
      deallog << "3d found manifold id 5 on subcell: " << found << std::endl;
    }
}

int
main()
{
  initlog();

  test<2, 2>();
  test<3, 3>();

  return 0;
}
