// ---------------------------------------------------------------------
//
// Copyright (C) 2020 - 2022 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE.md at
// the top level directory of deal.II.
//
// ---------------------------------------------------------------------


// Test TriaAccessor::reference_cell() on a 1D/2D/3D HEX mesh.

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>

#include "../tests.h"

using namespace dealii;

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
