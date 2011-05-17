//----------------------------  vector.cc  -------------------------------
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
//----------------------------  vector.cc  -------------------------------

// check serialization for Vector

#include "serialization.h"
#include <deal.II/lac/vector.h>
#include <boost/serialization/vector.hpp>

void test ()
{
  unsigned int n = 5;
  
  Vector<double> v1(n);
  Vector<double> v2(n);
  
  Vector<double> v3;
  
  for (unsigned int i = 0; i<n; ++i)
  {
    v1(i) = i*1.;
    v2(i) = i*1. + n*1.;
  }
   

  verify (v1, v2);
  
  verify (v1, v3);
}


int main ()
{
  std::ofstream logfile("vector/output");
  deallog << std::setprecision(3);
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  test ();

  deallog << "OK" << std::endl;
}
