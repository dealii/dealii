// ---------------------------------------------------------------------
//
// Copyright (C) 2018 - 2020 by the deal.II authors
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
