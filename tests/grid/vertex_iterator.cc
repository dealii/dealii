// ---------------------------------------------------------------------
//
// Copyright (C) 2015, 2016 by the deal.II authors
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


// test TriaAccessor<0,dim,spacedim>


#include "../tests.h"
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/grid_generator.h>

#include <fstream>


void test ()
{
  Triangulation<2> tria;
  GridGenerator::hyper_cube (tria);
  tria.refine_global (2);

  for (Triangulation<2>::vertex_iterator
       vertex_it = tria.begin_vertex(); vertex_it != tria.end_vertex();
       ++vertex_it)
    deallog << vertex_it->center() <<std::endl;
  deallog << std::endl;

  for (Triangulation<2>::active_cell_iterator
       cell = tria.begin_active(); cell != tria.end(); ++cell)
    {
      for (unsigned int i=0; i<4; ++i)
        deallog << cell->vertex_iterator(i)->center() <<std::endl;
      deallog << std::endl;
    }

  for (Triangulation<2>::active_cell_iterator
       cell = tria.begin_active(); cell != tria.end(); ++cell)
    for (unsigned int i=0; i<4; ++i)
      {
        for (unsigned int j=0; j<2; ++j)
          deallog << cell->line(i)->vertex_iterator(j)->center() <<std::endl;
        deallog << std::endl;
      }
}



int main ()
{
  std::ofstream logfile("output");
  deallog.attach(logfile);

  test ();

  return 0;
}
