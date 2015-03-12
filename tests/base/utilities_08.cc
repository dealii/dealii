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


// Utilities::int_to_string produced wrong results with numbers larger
// than 10 digits (including padding)

#include "../tests.h"
#include <iomanip>
#include <iomanip>
#include <fstream>
#include <cmath>

#include <deal.II/base/utilities.h>


void test ()
{
  deallog << Utilities::int_to_string (9,10) << std::endl;
  deallog << Utilities::int_to_string (99,10) << std::endl;
  deallog << Utilities::int_to_string (999,10) << std::endl;
  deallog << Utilities::int_to_string (9999,10) << std::endl;
  deallog << Utilities::int_to_string (99999,10) << std::endl;
  deallog << Utilities::int_to_string (999999,10) << std::endl;
  deallog << Utilities::int_to_string (9999999,10) << std::endl;
  deallog << Utilities::int_to_string (99999999,10) << std::endl;
  deallog << Utilities::int_to_string (999999999,10) << std::endl;
}




int main()
{
  std::ofstream logfile("output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  test ();
}
