//----------------------------------------------------------------------
//    $Id$
//    Version: $Name$ 
//
//    Copyright (C) 2007, 2008 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------------------------------------------------


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

std::ofstream logfile("grid_generator_01/output");


template<int dim>
void check_rect1 (unsigned int n, bool color, bool log)
{
  Point<dim> left;
  Point<dim> right;
  std::vector<unsigned int> subdivisions(dim);
  
  for (unsigned int d=0;d<dim;++d)
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
