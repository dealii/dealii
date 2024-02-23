// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2015 - 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// test Utilities::split_string_list with an empty string

#include <deal.II/base/utilities.h>

#include "../tests.h"


void
test()
{
  // test an empty string -- should yield a list of zero elements with
  // any delimiter
  {
    const char *p = "";
    deallog << Utilities::split_string_list(p).size() << std::endl;
    deallog << Utilities::split_string_list(p, ' ').size() << std::endl;
  }

  // also test a string that consists only of whitespace. this should
  // yield a list of zero elements even if (maybe not very usefully)
  // the delimiter is chosen as a whitespace itself
  {
    const char *p = "  ";
    deallog << Utilities::split_string_list(p).size() << std::endl;
    deallog << Utilities::split_string_list(p, ' ').size() << std::endl;
  }
}



int
main()
{
  initlog();

  test();
}
