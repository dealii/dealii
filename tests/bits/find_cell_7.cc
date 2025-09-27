// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2003 - 2021 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// find_active_cell_around_point() stops working if some cells are
// refined in 3d.  this is caused by a bug in
// GridTools::find_cells_adjacent_to_vertex that got fixed in r 25562.
//
// the bug here is the same as in find_cell_6 but when calling the
// function with hp:: arguments

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/mapping_q1.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/hp/fe_collection.h>
#include <deal.II/hp/mapping_collection.h>

#include "../tests.h"

bool
inside(Triangulation<3> &tria, Point<3> &p)
{
  for (Triangulation<3>::cell_iterator cell = tria.begin(0);
       cell != tria.end(0);
       ++cell)
    if (cell->point_inside(p))
      return true;

  return false;
}

void
check2()
{
  Triangulation<3> tria;
  GridIn<3>        gridIn;
  gridIn.attach_triangulation(tria);
  std::ifstream inputFile(SOURCE_DIR "/grids/grid.inp");
  gridIn.read_ucd(inputFile);

  Point<3> p2(304.767, -57.0113, 254.766);

  int idx = 0;
  for (Triangulation<3>::active_cell_iterator cell = tria.begin_active();
       cell != tria.end();
       ++cell, ++idx)
    {
      if (idx == 21)
        cell->set_refine_flag();
    }
  tria.execute_coarsening_and_refinement();

  deallog << inside(tria, p2) << std::endl;

  hp::MappingCollection<3> mappings;
  mappings.push_back(MappingQ<3>(1));
  mappings.push_back(MappingQ<3>(1));

  hp::FECollection<3> fes;
  fes.push_back(FE_Q<3>(1));
  fes.push_back(FE_Q<3>(1));

  DoFHandler<3> dof_handler(tria);
  dof_handler.distribute_dofs(fes);

  GridTools::find_active_cell_around_point(mappings,
                                           dof_handler,
                                           p2); // triggered exception
}


int
main()
{
  initlog();

  try
    {
      check2();
    }
  catch (const std::exception &exc)
    {
      // we shouldn't get here...
      deallog << "Caught an error..." << std::endl;
      deallog << exc.what() << std::endl;
    }
}
