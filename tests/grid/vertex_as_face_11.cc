// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
// Copyright (C) 2010 - 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Detailed license information governing the source code and contributions
// can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
//
// -----------------------------------------------------------------------------


// verify that we can do things like cell->face() in 1d as well. here:
// getting boundary indicators


#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include "../tests.h"



template <int spacedim>
void
test()
{
  Triangulation<1, spacedim> tria;
  GridGenerator::hyper_cube(tria);

  tria.begin_active()->face(0)->set_boundary_id(2);
  tria.begin_active()->face(1)->set_boundary_id(4);

  deallog << (int)tria.begin_active()->face(0)->boundary_id() << std::endl;
  deallog << (int)tria.begin_active()->face(1)->boundary_id() << std::endl;
}



int
main()
{
  initlog();

  test<1>();
  test<2>();

  return 0;
}
