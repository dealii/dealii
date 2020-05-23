// ---------------------------------------------------------------------
//
// Copyright (C) 2003 - 2020 by the deal.II authors
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



// find_active_cell_around_point should throw an exception if the
// point is outside. Test that.

#include <deal.II/base/function.h>
#include <deal.II/base/quadrature_lib.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include <time.h>

#include <iostream>
#include <list>
#include <sstream>
#include <string>

#include "../tests.h"


void
test()
{
  Triangulation<2> tr;
  GridGenerator::hyper_cube(tr);

  Point<2> p;
  p(0) = -0.1;
  p(1) = 0.5;

  MappingQ<2> mapping(1);

  try
    {
      GridTools::find_active_cell_around_point(mapping, tr, p);
    }
  catch (GridTools::ExcPointNotFound<2> &e)
    {
      deallog << "outside" << std::endl;
    }
  deallog << "done" << std::endl;
}

int
main(int argc, char **argv)
{
  initlog();

  test();

  return 0;
}
