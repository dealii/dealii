//---------------------------------------------------------------------------
//    $Id$
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
  Point<dim> p1;
  p1[0] = 2.;
  if (dim>1) p1[1] = -1.;
  if (dim>2) p1[2] = 0.;
  Point<dim> p2;
  p2[0] = 3.;
  if (dim>1) p2[1] = 2.;
  if (dim>2) p2[2] = 4.;
  Point<dim> p3;
  p3[0] = 2.;
  if (dim>1) p3[1] = 1.;
  if (dim>2) p3[2] = 4.;
  

                                   // loop without rotation 
  if (true)
    {
      deallog <<"moebius, no rotation" << std::endl;
      Triangulation<dim> tr;
      GridGenerator::moebius(tr, 20, 0, 10.0, 2.0);
      GridOut go;
      go.set_flags (GridOutFlags::Ucd(true));
      go.write_ucd(tr, out);
    }

                                   // loop with quarter rotation (1 * pi/2)
  if (true)
    {
      deallog << "---------------------------"<<std::endl
	      << "moebius, quarter rotation (1* Pi/2)" << std::endl;
      Triangulation<dim> tr;
      GridGenerator::moebius(tr, 20, 1, 10.0, 2.0);
      GridOut go;
      go.set_flags (GridOutFlags::Ucd(true));
      go.write_ucd(tr, out);
    }

                                   // loop with half rotation (2 * pi/2)
  if (true)
    {
      deallog << "---------------------------"<<std::endl
	      << "moebius, half rotation (2* Pi/2)" << std::endl;
      Triangulation<dim> tr;
      GridGenerator::moebius(tr, 20, 2, 10.0, 2.0);
      GridOut go;
      go.set_flags (GridOutFlags::Ucd(true));
      go.write_ucd(tr, out);
    }

                                   // loop with three quarter rotation (3 * pi/2)
  if (true)
    {
      deallog << "---------------------------"<<std::endl
	      << "moebius, three quarter rotation (3* Pi/2)" << std::endl;
      Triangulation<dim> tr;
      GridGenerator::moebius(tr, 20, 3, 10.0, 2.0);
      GridOut go;
      go.set_flags (GridOutFlags::Ucd(true));
      go.write_ucd(tr, out);
    }

                                   // loop with full rotation (1 * pi/2)
  if (true)
    {
      deallog << "---------------------------"<<std::endl
	      << "moebius, full rotation (2* Pi)" << std::endl;
      Triangulation<dim> tr;
      GridGenerator::moebius(tr, 20, 4, 10.0, 2.0);
      GridOut go;
      go.set_flags (GridOutFlags::Ucd(true));
      go.write_ucd(tr, out);
    }
}


int main()
{
  std::ofstream logfile("grid_generator_05/output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);
  
  test<3>(logfile);
}
