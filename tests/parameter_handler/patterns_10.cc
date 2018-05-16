// ---------------------------------------------------------------------
//
// Copyright (C) 2005 - 2017 by the deal.II authors
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
#include <deal.II/base/std_cxx14/memory.h>
#include <memory>

int
main()
{
  initlog();

  // create a pattern and match a string
  const auto &pattern = Patterns::Tuple(";", Patterns::Integer(), Patterns::Double(), Patterns::Anything());
  const std::string desc = pattern.description();

  deallog << desc << std::endl;

  std::string test = "5; 3.14; Ciao";

  if (pattern.match(test))
    deallog << "OK" << std::endl;
  else
    deallog << "Not OK" << std::endl;
}
