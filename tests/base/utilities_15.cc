// ---------------------------------------------------------------------
//
// Copyright (C) 2019 by the deal.II authors
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


// test Utilities::type_to_string()

#include <deal.II/base/utilities.h>

#include "../tests.h"

int
main()
{
  initlog();

  double   a = 1.0;
  int      b = 3;
  Point<2> c;

  deallog << Utilities::type_to_string(a) << std::endl
          << Utilities::type_to_string(b) << std::endl
          << Utilities::type_to_string(c) << std::endl;
}
