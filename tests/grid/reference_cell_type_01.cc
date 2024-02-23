// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2020 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// Test TriaAccessor::reference_cell() on a 1D/2D/3D HEX mesh.

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>

#include "../tests.h"


template <int dim>
void
test()
{
  Triangulation<dim> tria;
  GridGenerator::subdivided_hyper_cube(tria, 3);

  for (const auto &cell : tria.active_cell_iterators())
    deallog << static_cast<int>(cell->reference_cell()) << ' ';
  deallog << std::endl;

  if (dim != 1)
    for (const auto &face : tria.active_face_iterators())
      deallog << static_cast<int>(face->reference_cell()) << ' ';
  deallog << std::endl;

  for (const auto &cell : tria.active_cell_iterators())
    for (const auto &face : cell->face_iterators())
      deallog << static_cast<int>(face->reference_cell()) << ' ';
  deallog << std::endl;
}

int
main()
{
  initlog();

  {
    deallog.push("1D");
    test<1>();
    deallog.pop();
  }

  {
    deallog.push("2D");
    test<2>();
    deallog.pop();
  }

  {
    deallog.push("3D");
    test<3>();
    deallog.pop();
  }
}
