// ---------------------------------------------------------------------
//
// Copyright (C) 2010 - 2013 by the deal.II authors
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


// verify that we can do things like cell->face() in 1d as well. here:
// test cell->face()->at_boundary()


#include "../tests.h"
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/grid_generator.h>

#include <fstream>


template <int spacedim>
void test ()
{
  Triangulation<1,spacedim> tria;
  GridGenerator::hyper_cube (tria);

  deallog << "Coarse mesh:" << std::endl;
  deallog << "Left vertex=" << tria.begin_active()->face(0)->at_boundary() << std::endl;
  deallog << "Right vertex=" << tria.begin_active()->face(1)->at_boundary() << std::endl;

  tria.refine_global (2);

  for (typename Triangulation<1,spacedim>::active_cell_iterator
       cell = tria.begin_active();
       cell != tria.end(); ++cell)
    {
      deallog << "Cell: " << cell << std::endl;
      deallog << "Left vertex=" << cell->face(0)->at_boundary() << std::endl;
      deallog << "Right vertex=" << cell->face(1)->at_boundary() << std::endl;
    }
}



int main ()
{
  std::ofstream logfile("output");
  deallog.attach(logfile);
  deallog.depth_console(0);

  test<1> ();
  test<2> ();

  return 0;
}
