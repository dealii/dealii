// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2010 - 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// Test the testsuite framework. This test is supposed to segfault when on
// first try and succeed when invoked a second time.

#include <filesystem>

#include "../tests.h"


int
main()
{
  initlog();

  if (!std::filesystem::exists("segfault_encountered"))
    {
      // Dereference a null pointer and cause a segfault. Note that we need
      // to cast to a volatile int pointer here to prevent the compiler from
      // deleting the dereference.
      *static_cast<volatile int *>(nullptr) = 42;
    }

  deallog << "OK" << std::endl;
}
