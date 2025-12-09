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
  deallog << "Extract FiniteElement<" << dim << ", " << spacedim
          << "> from file: " << filename << std::endl;
  auto [fe, names] = VTKWrappers::vtk_to_finite_element<dim, spacedim>(
    std::string(SOURCE_DIR) + "/" + filename);

  deallog << "FiniteElement: " << fe->get_name() << std::endl
          << "Number of components: " << fe->n_components() << std::endl
          << "Field names: " << Patterns::Tools::to_string(names) << std::endl;
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
