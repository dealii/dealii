// ---------------------------------------------------------------------
//
// Copyright (C) 2003 - 2013 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------



// find_active_cell_around_point should throw an exception if the
// point is outside. Test that.

#include "../tests.h"

#include <stdio.h>
#include <cstdlib>

#include <base/quadrature_lib.h>
#include <fe/mapping_q.h>
#include <base/function.h>
#include <base/logstream.h>
#include <grid/tria.h>
#include <grid/grid_generator.h>
#include <grid/tria_accessor.h>
#include <grid/tria_iterator.h>
#include <grid/grid_tools.h>
#include <dofs/dof_handler.h>
#include <fe/fe_q.h>
#include <fe/fe_values.h>
#include <grid/grid_in.h>


#include <iostream>
#include <fstream>
#include <list>
#include <string>
#include <sstream>
#include <time.h>

using namespace dealii;

void test()
{
  Triangulation<2> tr;
  GridGenerator::hyper_cube(tr);
  
  Point< 2 > p;
  p(0) = -0.1;
  p(1) = 0.5;

  MappingQ<2> mapping(1);
  
  try
    {
      GridTools::find_active_cell_around_point (mapping, tr, p);
    }
  catch (GridTools::ExcPointNotFound<2> &e)
    {
      deallog << "outside" << std::endl;
    }
  deallog << "done" << std::endl;
}

int main (int argc, char **argv)
{
  initlog();

  test();

  return 0;
}
