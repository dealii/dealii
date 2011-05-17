//----------------------------  polynomial_hierarchical.cc  ---------------------------
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
//----------------------------  polynomial_hierarchical.cc  ---------------------------

// check serialization for Hierarchical

#include "serialization.h"
#include <deal.II/base/polynomial.h>
#include <boost/serialization/vector.hpp>

void test ()
{
  unsigned int degree1 = 3; 
  Polynomials::Hierarchical p1(degree1);
  
  unsigned int degree2 = 7;
  Polynomials::Hierarchical p2(degree2);
  
  verify (p1, p2);
}


int main ()
{
  std::ofstream logfile("polynomial_hierarchical/output");
  deallog << std::setprecision(3);
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  test ();

  deallog << "OK" << std::endl;
}
