// ---------------------------------------------------------------------
//
// Copyright (C) 1998 - 2013 by the deal.II authors
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



// test the functions of ConditionalOStream


#include "../tests.h"
#include <deal.II/base/conditional_ostream.h>
#include <fstream>
#include <iomanip>
#include <limits>


int main()
{
  std::ofstream logfile("output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  ConditionalOStream o(logfile, true);
  o << "Yes" << std::endl;
  deallog << o.is_active() << std::endl;

  o.set_condition (false);
  o << "No" << std::endl;
  deallog << o.is_active() << std::endl;
}
