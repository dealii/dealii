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


// test Utilities::trim

#include "../tests.h"
#include <iomanip>
#include <iomanip>
#include <fstream>
#include <cmath>

#include <deal.II/base/utilities.h>


void check(const std::string &input, const std::string &expected)
{
  deallog << "trim(\"" << input << "\") = \"" << Utilities::trim(input) << "\"" << std::endl;
  Assert(Utilities::trim(input) == expected, ExcInternalError());
}



void test ()
{
  check("Hello World", "Hello World");
  check("", "");
  check(" ", "");
  check("    ", "");
  check("  middle   ", "middle");
  check("left   ", "left");
  check("  right", "right");
  check("  multiple  words with spaces  ", "multiple  words with spaces");
}


int main()
{
  std::ofstream logfile("output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  test ();
}
