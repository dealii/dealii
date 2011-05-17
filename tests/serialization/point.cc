//----------------------------  tensor_base.cc  ---------------------------
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
//----------------------------  tensor_base.cc  ---------------------------

// check serialization for Tensor<1,dim>

#include "serialization.h"
#include <deal.II/base/point.h>


void test ()
{
  const unsigned int dim=3;

  Point<dim> p1(1.,2.,3.);

  Point<dim> p2(4.,5.,6.);

  verify (p1, p2);
}


int main ()
{
  std::ofstream logfile("point/output");
  deallog << std::setprecision(3);
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  test ();

  deallog << "OK" << std::endl;
}
