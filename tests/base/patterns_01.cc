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


#include "../tests.h"
#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/logstream.h>
#include <memory>

int main()
{
  std::ofstream logfile("output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  // create a pattern and let it
  // output its description
  Patterns::Integer pattern(-1,42);
  const std::string desc = pattern.description();

  // now let the same class re-create
  // a pattern object from the
  // description and verify that the
  // result is the same as what we
  // started out with
  Patterns::Integer *pattern2 = Patterns::Integer::create (desc);

  Assert (pattern2 != 0, ExcInternalError());
  Assert (desc == pattern2->description(), ExcInternalError());

  deallog << desc << std::endl;

  delete pattern2;
}
