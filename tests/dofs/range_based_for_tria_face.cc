// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
// Copyright (C) 2018 - 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Detailed license information governing the source code and contributions
// can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
//
// -----------------------------------------------------------------------------



// Check face range-based for loops for triangulations

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include <string>

#include "../tests.h"


template <int dim, int spacedim>
void
check()
{
  Triangulation<dim, spacedim> tr;
  GridGenerator::hyper_cube(tr);
  tr.refine_global(2);

  for (auto &face : tr.active_face_iterators())
    face->set_manifold_id(42);

  for (const auto &cell : tr.active_cell_iterators())
    for (const unsigned int face_n : GeometryInfo<dim>::face_indices())
      Assert(cell->face(face_n)->manifold_id() == 42, ExcInternalError());

  deallog << "OK" << std::endl;
}


int
main()
{
  initlog();

  check<2, 2>();
  check<2, 3>();
  check<3, 3>();
}
