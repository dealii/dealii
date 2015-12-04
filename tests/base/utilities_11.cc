// ---------------------------------------------------------------------
//
// Copyright (C) 2005 - 2015 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------


// This tests the implementation of Utilities::to_string for different types.
// Note that the floating point number output might be depend on the system.

#include "../tests.h"
#include <iomanip>
#include <fstream>
#include <cmath>

#include <deal.II/base/utilities.h>


void test ()
{
  unsigned long long int i = std::pow(2,33);
  deallog << Utilities::to_string(i) << std::endl;
  deallog << Utilities::to_string(i,11) << std::endl;

  unsigned long long j = std::pow(2,31);
  deallog << Utilities::to_string (j) << std::endl;

  int k = - std::pow(2,30);
  deallog << Utilities::to_string (k) << std::endl;
  deallog << Utilities::to_string (k,12) << std::endl;

  long long int l = - std::pow(2,35);
  deallog << Utilities::to_string (l) << std::endl;
  deallog << Utilities::to_string (l,13) << std::endl;

  float f (-3.14159265358979323846264338327950288419716939937510);
  deallog << Utilities::to_string (f) << std::endl;
  deallog << Utilities::to_string (f,13) << std::endl;

  double d (-3.14159265358979323846264338327950288419716939937510);
  deallog << Utilities::to_string (d) << std::endl;
  deallog << Utilities::to_string (d,20) << std::endl;

  long double ld (-3.14159265358979323846264338327950288419716939937510);
  deallog << Utilities::to_string (ld) << std::endl;
  deallog << Utilities::to_string (ld,24) << std::endl;


}




int main()
{
  std::ofstream logfile("output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  test ();
}
