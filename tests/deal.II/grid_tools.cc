//----------------------------  grid_test.cc  ---------------------------
//    $Id$
//    Version: $Name$ 
//
//    Copyright (C) 2001, 2002, 2003, 2004 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  grid_test.cc  ---------------------------


#include "../tests.h"
#include <base/logstream.h>
#include <grid/tria.h>
#include <grid/grid_generator.h>
#include <grid/grid_tools.h>
#include <grid/grid_out.h>

#include <fstream>


std::ofstream logfile("grid_tools.output");



// check GridTools::diameter
template <int dim>
void test1 ()
{
				   // test 1: hypercube
  if (true)
    {
      Triangulation<dim> tria;
      GridGenerator::hyper_cube(tria);

      for (unsigned int i=0; i<2; ++i)
	{
	  tria.refine_global(2);
	  deallog << dim << "d, "
		  << "hypercube diameter, "
		  << i*2
		  << " refinements: "
		  << GridTools::diameter (tria)
		  << std::endl;
	};
    };

				   // test 2: hyperball
  if (dim == 2)
    {
      Triangulation<dim> tria;
      GridGenerator::hyper_ball(tria, Point<dim>(), 1);

      for (unsigned int i=0; i<2; ++i)
	{
	  tria.refine_global(2);
	  deallog << dim << "d, "
		  << "hyperball diameter, "
		  << i*2
		  << " refinements: "
		  << GridTools::diameter (tria)
		  << std::endl;
	};
    };
}


// GridTools::transform
void test2 ()
{
  Triangulation<2> tria;
  GridGenerator::hyper_cube(tria);

  logfile << "Unchanged grid:" << std::endl;
  GridOut().write_gnuplot (tria, logfile);
  
  logfile << "Shifted grid:" << std::endl;
  const Point<2> shift(1,2);
  GridTools::shift (shift, tria);
  GridOut().write_gnuplot (tria, logfile);

  logfile << "Rotated grid:" << std::endl;
  GridTools::rotate (3.14159265258/4, tria);
  GridOut().write_gnuplot (tria, logfile);
}


int main ()
{
  logfile.precision(4);
  deallog.attach(logfile);
  deallog.depth_console(0);

  test1<1> ();
  test1<2> ();
  test1<3> ();

  test2 ();
  
  return 0;
}

