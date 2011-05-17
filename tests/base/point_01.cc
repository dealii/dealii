//----------------------------  point_01.cc  ---------------------------
//    $Id$
//    Version: $Name$ 
//
//    Copyright (C) 2010 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  point_01.cc  ---------------------------

// test for Point::operator()

#include "../tests.h"
#include <deal.II/base/point.h>
#include <deal.II/base/logstream.h>
#include <deal.II/lac/vector.h>
#include <fstream>
#include <iomanip>


template <int dim>
void check ()
{
  Point<dim> p;
  for (unsigned int i=0; i<dim; ++i)
    p[i] = i;

  for (unsigned int i=0; i<dim; ++i)
    deallog << p(i) << ' ';
  deallog << std::endl;
}

int main ()
{
  std::ofstream logfile("point_01/output");
  deallog << std::setprecision(3);
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  check<1>();
  check<2>();
  check<3>();
}
