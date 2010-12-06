//----------------------------  polynomial_lagrange_equidistant.cc  ---------------
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
//----------------------------  polynomial_lagrange_equidistant.cc  ---------------

// check serialization for LagrangeEquidistant

#include "serialization.h"
#include <base/polynomial.h>
#include <boost/serialization/vector.hpp>

void test ()
{
  unsigned int n1 = 3;
  unsigned int support_point1 = 1;
  
  Polynomials::LagrangeEquidistant p1(n1, support_point1);

  unsigned int n2 = 4;
  unsigned int support_point2 = 2;
  
  Polynomials::LagrangeEquidistant p2(n2, support_point2);

  verify (p1, p2);
}


int main ()
{
  std::ofstream logfile("polynomial_lagrange_equidistant/output");
  deallog << std::setprecision(3);
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  test ();

  deallog << "OK" << std::endl;
}
