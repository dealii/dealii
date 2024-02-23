// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2015 - 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// test TriaAccessor<0,dim,spacedim>


#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include "../tests.h"



void
test()
{
  Triangulation<2> tria;
  GridGenerator::hyper_cube(tria);
  tria.refine_global(2);

  for (Triangulation<2>::vertex_iterator vertex_it = tria.begin_vertex();
       vertex_it != tria.end_vertex();
       ++vertex_it)
    deallog << vertex_it->center() << std::endl;
  deallog << std::endl;

  for (Triangulation<2>::active_cell_iterator cell = tria.begin_active();
       cell != tria.end();
       ++cell)
    {
      for (unsigned int i = 0; i < 4; ++i)
        deallog << cell->vertex_iterator(i)->center() << std::endl;
      deallog << std::endl;
    }

  for (Triangulation<2>::active_cell_iterator cell = tria.begin_active();
       cell != tria.end();
       ++cell)
    for (unsigned int i = 0; i < 4; ++i)
      {
        for (unsigned int j = 0; j < 2; ++j)
          deallog << cell->line(i)->vertex_iterator(j)->center() << std::endl;
        deallog << std::endl;
      }
}



int
main()
{
  initlog();

  test();

  return 0;
}
