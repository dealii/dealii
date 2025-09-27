// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2010 - 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// verify that we can do things like cell->face() in 1d as well. here:
// setting boundary indicators


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

  std::vector<types::boundary_id> boundary_ids = tria.get_boundary_ids();

  for (unsigned int i = 0; i < boundary_ids.size(); ++i)
    deallog << (int)boundary_ids[i] << std::endl;
}



int
main()
{
  initlog();

  test<1>();
  test<2>();

  return 0;
}
