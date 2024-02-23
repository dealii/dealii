// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2012 - 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// test Utilities::replace_in_string

#include <deal.II/base/utilities.h>

#include "../tests.h"

void
check(const std::string in,
      const std::string from,
      const std::string to,
      std::string       out)
{
  std::string result = Utilities::replace_in_string(in, from, to);
  if (result != out)
    {
      deallog << "in='" << in << "' from='" << from << "' to='" << to
              << "' result='" << result << "' != '" << out << "'" << std::endl;
    }
}


void
test()
{
  check("wie geht es dir?", "dir", "euch", "wie geht es euch?");
  check("empty from", "", "abc", "empty from");
  check("eins zwei drei", "ei", "", "ns zw dr");
  check("eins zwei drei", "zwei", "zweiundvierzig", "eins zweiundvierzig drei");
  check("wer das liest ist doof", "das liest ", "", "wer ist doof");
  check("string string", "string", "", " ");
  check(" same is same", " same", " same", " same is same");
  check("  ", " ", "", "");
  check("", " ", "abc", "");
  check("small SMALL", "LL", "ll", "small SMAll");
  deallog << "OK" << std::endl;
}


int
main()
{
  initlog();

  test();
}
