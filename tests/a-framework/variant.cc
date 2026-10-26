// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
// Copyright (C) 2026 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Detailed license information governing the source code and contributions
// can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
//
// -----------------------------------------------------------------------------



// test the testsuite framework. this test is supposed to miscompare with
// the main comparison file but to match the variant comparison file.
// Furthermore, the output is normalized prior to the comparison.

#include "../tests.h"


int
main()
{
  initlog();

  deallog << "OK" << std::endl;
  deallog << "-0.000" << std::endl;
  deallog << "JobId 2026" << std::endl;
}
