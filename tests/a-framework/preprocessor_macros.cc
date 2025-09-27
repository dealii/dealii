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
