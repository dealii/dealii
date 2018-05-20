// ---------------------------------------------------------------------
//
// Copyright (C) 2005 - 2017 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------

// test Utilities::replace_in_string

#include "../tests.h"

#include <deal.II/base/utilities.h>

void
check(const std::string in,
      const std::string from,
      const std::string to,
      std::string       out)
{
  std::string result = Utilities::replace_in_string(in, from, to);
  if(result != out)
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
