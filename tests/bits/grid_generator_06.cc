//---------------------------------------------------------------------------
//    $Id: grid_generator_05.cc 23710 2011-05-17 04:50:10Z bangerth $
//    Version: $Name$ 
//
//    Copyright (C) 2005, 2006, 2007 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------


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
void test(std::ostream& out)
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
  std::ofstream logfile("grid_generator_06/output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);
  
  test<3>(logfile);
}
