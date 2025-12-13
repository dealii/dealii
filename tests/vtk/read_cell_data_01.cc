// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

// Test reading of cell data contained in a vtk file into a deal.II
// Vector<double> using the VTKWrappers::read_cell_data function.

#include <deal.II/base/point.h>
#include <deal.II/base/utilities.h>

#include <deal.II/grid/tria.h>

#include <deal.II/lac/vector.h>

#include <deal.II/vtk/utilities.h>

#include "../tests.h"

template <int dim, int spacedim>
void
test(const std::string &filename)
{
  deallog << "Reading Triangulation<" << dim << ", " << spacedim
          << "> from file: " << filename << std::endl;

  Triangulation<dim, spacedim> triangulation;
  VTKWrappers::read_tria(std::string(SOURCE_DIR) + "/" + filename,
                         triangulation);

  Vector<double> data;
  VTKWrappers::read_cell_data(std::string(SOURCE_DIR) + "/" + filename,
                              "center_x",
                              data);

  Vector<double> vector_data;
  VTKWrappers::read_cell_data(std::string(SOURCE_DIR) + "/" + filename,
                              "center_xyz",
                              vector_data);

  deallog << "Scalar cell data size: " << data.size() << std::endl;
  deallog << "Vector cell data size: " << vector_data.size() << std::endl;

  for (const auto &cell : triangulation.active_cell_iterators())
    {
      const unsigned int cell_index = cell->active_cell_index();
      deallog << "Cell " << cell_index << ": "
              << " cell center: " << cell->center()
              << " data center_x = " << data(cell_index)
              << ", data center_xyz = (" << vector_data(cell_index * 3 + 0)
              << ", " << vector_data(cell_index * 3 + 1) << ", "
              << vector_data(cell_index * 3 + 2) << ")" << std::endl;
    }
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
