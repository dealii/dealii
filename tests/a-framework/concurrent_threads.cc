// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
// Copyright (C) 2018 - 2022 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Detailed license information governing the source code and contributions
// can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
//
// -----------------------------------------------------------------------------

// test the testsuite framework.

// Verify that the "threads=N" test annotation actually works.

#include "../tests.h"

int
main()
{
  std::cout << "Number of threads " << dealii::MultithreadInfo::n_threads()
            << std::endl;
  return 0;
}
