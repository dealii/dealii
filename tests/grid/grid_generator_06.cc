// ---------------------------------------------------------------------
//
// Copyright (C) 2005 - 2013 by the deal.II authors
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



// Test GridGenerator::moebius

#include "../tests.h"
#include <deal.II/base/logstream.h>
#include <deal.II/base/tensor.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/fe/fe_q.h>

#include <fstream>
#include <iomanip>


template<int dim>
void test(std::ostream &out)
{
  Triangulation<2> triangulation;
  Triangulation<3> tr;
  GridGenerator::hyper_rectangle (triangulation, Point<2>(0,0), Point<2>(1,1), true);

  GridGenerator::extrude_triangulation(triangulation, 3, 2.0, tr);
  GridOut go;
  go.set_flags (GridOutFlags::Ucd(true));
  go.write_ucd(tr, out);
}


int main()
{
  std::ofstream logfile("output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  test<3>(logfile);
}
