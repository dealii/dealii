//----------------------------  polynomial_monomial.cc  ---------------------------
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
//----------------------------  polynomial_monomial.cc  ---------------------------

// check serialization for Monomial

#include "serialization.h"
#include <deal.II/base/polynomial.h>
#include <boost/serialization/vector.hpp>


void test ()
{
  unsigned int n1 = 3;
  double coefficient1 = 5.;
  
  Polynomials::Monomial<double> m1(n1, coefficient1);

  unsigned int n2 = 3;
  double coefficient2 = 2.;
  
  Polynomials::Monomial<double> m2(n2, coefficient2);

  verify (m1, m2);
}


int main ()
{
  std::ofstream logfile("polynomial_monomial/output");
  deallog << std::setprecision(3);
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  test ();

  deallog << "OK" << std::endl;
}
