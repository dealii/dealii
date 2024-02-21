// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2017 - 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

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
