// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
// Copyright (C) 2023 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Detailed license information governing the source code and contributions
// can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
//
// -----------------------------------------------------------------------------

// Test reading of a vtk file into a deal.II triangulation using
// VTKWrappers::read_vtk

#include <deal.II/base/point.h>
#include <deal.II/base/utilities.h>

#include <deal.II/grid/tria.h>

#include <deal.II/vtk/utilities.h>

#include "../tests.h"

template <int dim, int spacedim>
void
test(const std::string &filename)
{
  Triangulation<dim, spacedim> triangulation;
  deallog << "Reading Triangulation<" << dim << ", " << spacedim
          << "> from file: " << filename << std::endl;
  VTKWrappers::read_tria(std::string(SOURCE_DIR) + "/" + filename,
                         triangulation);

  deallog << "Number of active cells: " << triangulation.n_active_cells()
          << std::endl;
}

int
main()
{
  initlog();
  test<1, 1>("data/data_1_1.vtk");
  test<1, 2>("data/data_1_1.vtk");
  test<1, 3>("data/data_1_1.vtk");

  test<2, 2>("data/data_2_2_simplex.vtk");
  test<2, 3>("data/data_2_2_simplex.vtk");

  test<2, 2>("data/data_2_2_quad.vtk");
  test<2, 3>("data/data_2_2_quad.vtk");

  test<3, 3>("data/data_3_3_tets.vtk");
  test<3, 3>("data/data_3_3_hexes.vtk");
  return 0;
}
