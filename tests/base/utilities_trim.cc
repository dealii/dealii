// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2005 - 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// test Utilities::trim. Note that deallog does not like being given '\n's in
// the middle of strings, so the output file is almost empty.

#include <deal.II/base/utilities.h>

#include "../tests.h"



void
check(const std::string &input, const std::string &expected)
{
  AssertThrow(Utilities::trim(input) == expected, ExcInternalError());
}



void
test()
{
  check("\r\nHello World\r\n\r", "Hello World");
  check("", "");
  check(" ", "");
  check("    ", "");
  check("  \r\n\v\t\r\n\f  ", "");
  check("  \rmiddle\r\n   ", "middle");
  check("left\v\t\r\n   ", "left");
  check("  \n\v\f\r\nright", "right");
  check(" \t\v\f\t\r\n\r multiple  words with spaces  \v\f\n",
        "multiple  words with spaces");

  deallog << "OK" << std::endl;
}


int
main()
{
  initlog();

  test();
}
