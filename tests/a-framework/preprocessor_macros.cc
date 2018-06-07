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
// The full text of the license can be found in the file LICENSE.md at
// the top level directory of deal.II.
//
// ---------------------------------------------------------------------



// test the testsuite framework.

// This test checks whether the cmake configuration correctly exports
// preprocessor definitions, in particular whether the DEBUG maro is
// defined in debug mode.

#include "../tests.h"

int
main()
{
  initlog();
#ifdef DEBUG
  deallog << "DEBUG defined!" << std::endl;
#else
  deallog << "DEBUG undefined!" << std::endl;
#endif
}
