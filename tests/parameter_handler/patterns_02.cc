// ---------------------------------------------------------------------
//
// Copyright (C) 2017 - 2018 by the deal.II authors
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

// verify that Patterns::List objects can be copied

#include <deal.II/base/parameter_handler.h>

#include <memory>

#include "../tests.h"

int
main()
{
  initlog();

  {
    // create one pattern and a copy of it
    Patterns::List list(Patterns::Integer(-1, 42), 2, 3);
    Patterns::List list2(list);

    // both now go out of scope -- ensure that their destruction does
    // not lead to memory corruption
  }
  deallog << "OK" << std::endl;
}
