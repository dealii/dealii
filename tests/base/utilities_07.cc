// ---------------------------------------------------------------------
//
// Copyright (C) 2005 - 2014 by the deal.II authors
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


// verify that Utilities::string_to_double actually catches errors

#include "../tests.h"
#include <iomanip>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <sstream>

#include <deal.II/base/utilities.h>

using namespace dealii;




void verify (const std::string &s)
{
  bool exception_caught = false;
  try
    {
      Utilities::string_to_double(s);
    }
  catch (...)
    {
      exception_caught = true;
    }
  Assert (exception_caught == true, ExcMessage ("Function is broken!"));

  deallog << "Done correctly: " << s << std::endl;
}




int main()
{
  std::ofstream logfile("output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  verify ("abc");
  verify ("1.23.4");
  verify ("1 23 4");
  verify ("123abc");
}
