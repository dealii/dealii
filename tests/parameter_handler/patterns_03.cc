// ---------------------------------------------------------------------
//
// Copyright (C) 2017 by the deal.II authors
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

// verify that Patterns::List objects can be copied

#include "../tests.h"
#include <deal.II/base/parameter_handler.h>
#include <memory>

int
main()
{
  initlog();

  {
    // create one pattern and a copy of it
    Patterns::Map map(
      Patterns::Integer(-1, 42), Patterns::Double(-1e9, 1e9), 2, 3);
    Patterns::Map map2(map);

    // both now go out of scope -- ensure that their destruction does
    // not lead to memory corruption
  }
  deallog << "OK" << std::endl;
}
