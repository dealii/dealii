//----------------------------  polynomial.cc  ---------------------------
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
//----------------------------  polynomial.cc  ---------------------------

// check serialization for Polynomial

#include "serialization.h"
#include <boost/serialization/vector.hpp>
#include <deal.II/base/polynomial.h>


void test ()
{
  double c1[3] = {1., 2., 3.};
  std::vector<double> coefficients1(c1, &c1[3]);
  
  Polynomials::Polynomial<double> p1(coefficients1);

  double c2[3] = {4., 5., 6.};
  std::vector<double> coefficients2(c2, &c2[3]);
  
  Polynomials::Polynomial<double> p2(coefficients2);
  
  Polynomials::Polynomial<double> p3;

  verify (p1, p2);
  
  verify (p1, p3);
}


int main ()
{
  std::ofstream logfile("polynomial/output");
  deallog << std::setprecision(3);
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  test ();

  deallog << "OK" << std::endl;
}
