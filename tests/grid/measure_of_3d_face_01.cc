// ---------------------------------------------------------------------
//
// Copyright (C) 2014 by the deal.II authors
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


// Test by Kevin Drzycimski: compute the measure of the faces of a 3d cell

#include "../tests.h"
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_tools.h>
#include <fstream>
#include <iomanip>


// move backward two adjacent vertices of the top face up by one
// unit. all faces remain flat this way
Point<3> distort_planar(Point<3> p) 
{
  if (p(1) > 0.5 && p(2) > 0.5) 
    {
      p(1) += 1;
    }
  return p;
}


// lift two opposite vertices of the top face up by one unit to create
// a saddle surface
Point<3> distort_twisted(Point<3> p) 
{
  if (p(2) > 0.5 && (p(0) > 0.5 ^ p(1) > 0.5)) 
    {
      p(2) += 1;
    }
  return p;
}



void test ()
{
  Triangulation<3> tria;
  GridOut gridout;

  deallog << "Planar\n";
  GridGenerator::hyper_cube(tria);
  GridTools::transform(&distort_planar, tria);
  gridout.write_eps(tria, deallog.get_file_stream());

  double measure_planar[] = {1.5, 1.5, 1, ::sqrt(2), 1, 2};
  deallog << "Face\tExact\tMeasure" << std::endl;
  Triangulation<3>::active_cell_iterator cell = tria.begin_active();
  for (int i = 0; i < 6; ++i)
    {
      double m = cell->face(i)->measure();
      deallog << i << '\t' << measure_planar[i] << '\t' << m << std::endl;
    }
 
  deallog << "Twisted\n";
  tria.clear();
  GridGenerator::hyper_cube(tria);
  GridTools::transform(&distort_twisted, tria);
  gridout.write_eps(tria, deallog.get_file_stream());

  double measure_twisted[] = {1.5, 1.5, 1.5, 1.5, 1, 5./3};
  deallog << "Face\tExact\tMeasure" << std::endl;
  cell = tria.begin_active();
  for (int i = 0; i < 6; ++i)
    {
      double m;
      try
	{
	  m = cell->face(i)->measure();
	}
      catch (...)
	{
	  m = std::numeric_limits<double>::quiet_NaN();
	}
      deallog << i << '\t' << measure_twisted[i] << '\t' << m << std::endl;
    }

  deallog << std::endl;
}



int main()
{
  std::ofstream logfile ("output");
  deallog << std::setprecision (5);
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  // run the tests but continue when finding an exception: we will try
  // out a distorted cell for which TriaAccessor::measure() will error
  // out, but we'd still like to print and check the data for the
  // remaining faces
  deal_II_exceptions::disable_abort_on_exception();
  
  test ();
}

