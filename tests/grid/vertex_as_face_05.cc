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
// test DoFHandler accessors for the same thing as the _01 test


#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>

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

  FESystem<1, spacedim>   fe(FE_Q<1, spacedim>(2), 2, FE_Q<1, spacedim>(1), 1);
  DoFHandler<1, spacedim> dof_handler(tria);
  dof_handler.distribute_dofs(fe);

  deallog << "Coarse mesh:" << std::endl;
  deallog << "Left vertex=" << dof_handler.begin_active()->face(0)->vertex(0)
          << std::endl;
  deallog << "Right vertex=" << dof_handler.begin_active()->face(1)->vertex(0)
          << std::endl;

  tria.refine_global(2);
  dof_handler.distribute_dofs(fe);

  for (typename DoFHandler<1, spacedim>::active_cell_iterator cell =
         dof_handler.begin_active();
       cell != dof_handler.end();
       ++cell)
    {
      deallog << "Cell: " << cell << std::endl;
      deallog << "Left vertex=" << cell->face(0)->vertex(0) << std::endl;
      deallog << "Right vertex=" << cell->face(1)->vertex(0) << std::endl;
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
