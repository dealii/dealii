// ---------------------------------------------------------------------
//
// Copyright (C) 2005 - 2020 by the deal.II authors
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
