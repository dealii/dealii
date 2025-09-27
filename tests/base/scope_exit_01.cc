// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// Check that ScopeExit works in both regular return statements and
// via `throw` statements.

#include <deal.II/base/scope_exit.h>

#include "../tests.h"


void
f1()
{
  ScopeExit s([]() { deallog << "Exiting function regularly." << std::endl; });

  return;
}


void
f2()
{
  ScopeExit s([]() { deallog << "Exiting function via throw." << std::endl; });

  throw 123;
}


void
f3()
{
  ScopeExit s([]() { deallog << "Exiting function regularly." << std::endl; });

  try
    {
      throw 123;
    }
  catch (...)
    {
      // just eat the exception and return regularly
      deallog << "Caught exception. Making it go away." << std::endl;
    }
}



int
main()
{
  initlog();

  f1();

  try
    {
      f2();
    }
  catch (int i)
    {
      deallog << "Caught int=" << i << std::endl;
    }

  f3();
}
