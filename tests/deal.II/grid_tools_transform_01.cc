// ---------------------------------------------------------------------
//
// Copyright (C) 2001 - 2013 by the deal.II authors
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
#include <deal.II/base/logstream.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/grid_out.h>

#include <fstream>
#include <iomanip>

template<int dim>
Point<dim> trans_func(Point<dim> & p);

template<>
Point<2> trans_func(Point<2> & p)
{
  Point<2> r(p(0)+p(1)*p(1),p(1));
  return r;
}
template<>
Point<3> trans_func(Point<3> & p)
{
  Point<3> r(p(0)+p(1)*p(1),p(1),p(2));
  return r;
}



template <int dim>
void test ()
{
  deallog << "dim = " << dim << std::endl;
  
  Triangulation<dim> tria;
  GridGenerator::hyper_cube(tria);
  tria.refine_global(1);
  tria.begin_active()->set_refine_flag();
  tria.execute_coarsening_and_refinement();
  

  deallog << "Unchanged grid:" << std::endl;
  GridOut().write_gnuplot (tria, deallog.get_file_stream());

  {
    std::ofstream f("grid1");
    GridOut().write_gnuplot (tria, f);
  }

  GridTools::transform(trans_func<dim>, tria);
  
  deallog << "transformed grid:" << std::endl;
  GridOut().write_gnuplot (tria, deallog.get_file_stream());
  {
    std::ofstream f("grid2");
    GridOut().write_gnuplot (tria, f);
  }
}


int main ()
{
  initlog();
  deallog.depth_console(0);

  test<2> ();
  test<3> ();

  return 0;
}

