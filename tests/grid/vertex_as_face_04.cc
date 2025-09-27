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
// test vertex location


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

  deallog << "Coarse mesh:" << std::endl;
  deallog << "Left vertex=" << tria.begin_active()->face(0)->boundary_id()
          << std::endl;
  deallog << "Right vertex=" << tria.begin_active()->face(1)->boundary_id()
          << std::endl;

  tria.refine_global(2);

  for (typename Triangulation<1, spacedim>::active_cell_iterator cell =
         tria.begin_active();
       cell != tria.end();
       ++cell)
    {
      deallog << "Cell: " << cell << std::endl;
      deallog << "Left vertex=" << cell->face(0)->boundary_id() << std::endl;
      deallog << "Right vertex=" << cell->face(1)->boundary_id() << std::endl;
    }
}



int
main()
{
  initlog();

  test<1>();
  test<2>();

  return 0;
}
