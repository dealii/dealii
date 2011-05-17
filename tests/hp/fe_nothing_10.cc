//----------------------------  fe_nothing_10.cc  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2009, 2011 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  fe_nothing_10.cc  ---------------------------


// an extract of _09 that failed at the time of writing the test
// with an internal error


#include "../tests.h"
#include <deal.II/base/logstream.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/fe/fe_nothing.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/hp/fe_values.h>


#include <fstream>


template <int dim>
void test ()
{
  FESystem<dim> fe(FE_Nothing<dim>(), 2);
  FEValues<dim> fe_values (fe, QGauss<dim>(2), update_values);

  deallog << "OK" << std::endl;
}



int main ()
{
  std::ofstream logfile("fe_nothing_10/output");
  logfile.precision(2);

  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  test<1> ();
  test<2> ();
  test<3> ();

  deallog << "OK" << std::endl;
}
