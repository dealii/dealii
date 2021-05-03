// ---------------------------------------------------------------------
//
// Copyright (C) 2005 - 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE.md at
// the top level directory of deal.II.
//
// ---------------------------------------------------------------------

// Check to_string and to_value

#include <deal.II/base/parameter_handler.h>

#include <memory>

#include "../tests.h"

using dealii::Patterns::Tools::to_string;
using dealii::Patterns::Tools::to_value;

int
main()
{
  initlog();

  auto a = std::make_tuple(1, std::string("ciao"));

  auto s = to_string(a);
  to_value("2 : mondo", a);

  deallog << "From: " << s << " to " << to_string(a) << std::endl;
}
