// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
// Copyright (C) 2000 - 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Detailed license information governing the source code and contributions
// can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
//
// -----------------------------------------------------------------------------

// test the functions of JobIdentifier::base_name

#include "../tests.h"

int
main()
{
  initlog();

  deallog << JobIdentifier::get_dealjobid().base_name("mypath/test.cc")
          << std::endl;
  deallog << JobIdentifier::get_dealjobid().base_name("/foo.bar/mypath/test.cc")
          << std::endl;

  return 0;
}
