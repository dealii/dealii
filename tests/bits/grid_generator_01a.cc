//---------------------------------------------------------------------------
//    $Id$
//    Version: $Name$ 
//
//    Copyright (C) 2005, 2006, 2013 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------


// HalfHyperSphereBoundary got refining the edges of the flat surface edges
// wrong in 3d.

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

  HalfHyperBallBoundary<dim> boundary_description (p1, 3);
  GridOut go;
  GridOut::OutputFormat format = GridOut::gnuplot;
  
    {
      deallog << "half_hyper_ball" << std::endl;
      Triangulation<dim> tr;
      GridGenerator::half_hyper_ball(tr, p1, 3.);
      tr.set_boundary (0, boundary_description);

      tr.refine_global (2);
      go.write(tr, out, format);
    }  
}


int main()
{
  std::ofstream logfile("grid_generator_01a/output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);
  
  deallog.push("3d");
  test<3>(logfile);
  deallog.pop();
}
