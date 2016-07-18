// ---------------------------------------------------------------------
//
// Copyright (C) 2005 - 2015 by the deal.II authors
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



// Test if refining a quarter_hyper_shell using HyperBallBoundary yields
// correct results in 2d and 3d.

#include "../tests.h"
#include <deal.II/base/logstream.h>
#include <deal.II/base/tensor.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria_boundary_lib.h>
#include <deal.II/grid/grid_out.h>

#include <fstream>
#include <iomanip>


template<int dim>
void test(std::ostream &out)
{
  Point<dim> p1;
  p1[0] = 2.;

  HyperBallBoundary<dim> boundary_description (p1, 3);
  GridOut go;
  GridOut::OutputFormat format = GridOut::gnuplot;

  {
    deallog << "quarter_hyper_ball" << std::endl;
    Triangulation<dim> tr;
    GridGenerator::quarter_hyper_ball(tr, p1, 3.);
    tr.set_boundary (0, boundary_description);

    tr.refine_global (2);
    go.write(tr, out, format);
  }
}


int main()
{
  std::ofstream logfile("output");
  deallog.attach(logfile);
  deallog.threshold_double(1.e-10);

  deallog.push("2d");
  test<2>(logfile);
  deallog.pop();

  deallog.push("3d");
  test<3>(logfile);
  deallog.pop();
}
