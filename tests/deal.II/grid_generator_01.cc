// ---------------------------------------------------------------------
//
// Copyright (C) 2007 - 2013 by the deal.II authors
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



#include "../tests.h"
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/base/logstream.h>
#include <cmath>
#include <cstdlib>

#include <fstream>
#include <iostream>
#include <iomanip>

std::ofstream logfile("output");


template<int dim>
void check_rect1 (unsigned int n, bool color, bool log)
{
  Point<dim> left;
  Point<dim> right;
  std::vector<unsigned int> subdivisions(dim);

  for (unsigned int d=0; d<dim; ++d)
    {
      left(d) = -1.;
      right(d) = d+2;
      subdivisions[d] = n*(d+3);
    }
  Triangulation<dim> tria;
  GridGenerator::subdivided_hyper_rectangle(tria, subdivisions, left, right, color);

  GridOut grid_out;
  if (dim == 2)
    {
      if (log)
        grid_out.write_xfig (tria, logfile);
      else
        grid_out.write_xfig (tria, std::cout);
    }
  else
    {
      if (log)
        grid_out.write_dx (tria, logfile);
      else
        grid_out.write_dx (tria, std::cout);
    }
}


int main()
{
  check_rect1<2> (1, true, true);
  check_rect1<2> (3, true, true);
  check_rect1<3> (1, true, true);
}
