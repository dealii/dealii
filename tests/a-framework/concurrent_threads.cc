// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2018 - 2022 by the deal.II authors
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

// Verify that the "threads=N" test annotation actually works.

#include "../tests.h"

int
main()
{
  std::cout << "Number of threads " << dealii::MultithreadInfo::n_threads()
            << std::endl;
  return 0;
}
