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


#include "../tests.h"
#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/logstream.h>
#include <memory>

int main()
{
  initlog();

  // create a pattern and try to convert to and from string
  Patterns::Selection pattern("one|two|three");

  std::string a = "one";

  deallog << "Selection to string: "
          << pattern.to_string(a) << std::endl;

  pattern.to_value("two", a);
  deallog << "To value string: "
          << a << std::endl;

  return 0;
}
