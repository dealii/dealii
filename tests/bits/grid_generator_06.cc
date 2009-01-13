//---------------------------------------------------------------------------
//    $Id$
//    Version: $Name$ 
//
//    Copyright (C) 2005, 2006, 2009 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------


// Test GridGenerator::half_hyper_shell in 2d and 3d, making sure that we get
// the boundary right. in 3d, we need to specify the inner and outer radii

#include "../tests.h"
#include <base/logstream.h>
#include <base/tensor.h>
#include <grid/tria.h>
#include <grid/grid_generator.h>
#include <grid/grid_out.h>
#include <grid/tria_boundary_lib.h>

#include <fstream>
#include <iomanip>


template<int dim>
void test_1(std::ostream& out)
{
  Point<dim> p1;
  p1[0] = 2.;
  if (dim>1) p1[1] = -1.;
  if (dim>2) p1[2] = 0.;
  
  GridOut go;
  deallog << "half_hyper_shell without specifying radii" << std::endl;
  Triangulation<dim> tr;
  GridGenerator::half_hyper_shell(tr, p1, 4., 6.);

  if (dim == 2)
    {
				       // test with now specifying the inner and
				       // outer radius
      static HalfHyperShellBoundary<dim> boundary (p1);
      tr.set_boundary (0, boundary);

      tr.refine_global (2);
  
      go.write_gnuplot(tr, out);
    }
}

template<int dim>
void test_2(std::ostream& out)
{
  Point<dim> p1;
  p1[0] = 2.;
  if (dim>1) p1[1] = -1.;
  if (dim>2) p1[2] = 0.;
  
  GridOut go;
  deallog << "half_hyper_shell with specifying radii" << std::endl;
  Triangulation<dim> tr;
  GridGenerator::half_hyper_shell(tr, p1, 4., 6.);

				   // test with now specifying the inner and
				   // outer radius
  static HalfHyperShellBoundary<dim> boundary (p1, 4., 6.);
  tr.set_boundary (0, boundary);

  tr.refine_global (1);
  
  go.write_gnuplot(tr, out);
}


int main()
{
  std::ofstream logfile("grid_generator_06/output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);
  
  deallog.push("2d");
  test_1<2>(logfile);
  test_2<2>(logfile);
  deallog.pop();
  deallog.push("3d");
  test_1<3>(logfile);
  test_2<3>(logfile);
  deallog.pop();
}
