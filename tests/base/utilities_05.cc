// ---------------------------------------------------------------------
//
// Copyright (C) 2005 - 2013 by the deal.II authors
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


// Utilities::get_integer_at_position

#include "../tests.h"
#include <iomanip>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <sstream>

#include <deal.II/base/utilities.h>

using namespace dealii;




void test ()
{
  int number = 5;
  for (unsigned int i=0; i<7; ++i)
    {
      std::ostringstream s;
      s << "test test" << number << "test test";

      Assert (Utilities::get_integer_at_position (s.str(),
                                                  9).first
              == number,
              ExcInternalError());
      Assert (Utilities::get_integer_at_position (s.str(),
                                                  9).second
              == i+1,
              ExcInternalError());

      deallog << i << ' ' << Utilities::get_integer_at_position (s.str(),
                                                                 9).first
              << std::endl;

      number = number*10 + i;
    }
}




int main()
{
  std::ofstream logfile("output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  test ();
}
