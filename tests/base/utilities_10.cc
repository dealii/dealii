// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
// Copyright (C) 2015 - 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Detailed license information governing the source code and contributions
// can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
//
// -----------------------------------------------------------------------------


// test Utilities::split_string_list with a string that only contains
// the delimiter and, possibly, spaces

#include <deal.II/base/utilities.h>

#include "../tests.h"


void
test()
{
  // verify the documented behavior of eating trailing delimiters
  {
    deallog << Utilities::split_string_list(",").size() << std::endl;
    deallog << Utilities::split_string_list(" , ").size() << std::endl;
  }

  {
    deallog << Utilities::split_string_list(",,").size() << std::endl;
    deallog << Utilities::split_string_list(" , , ").size() << std::endl;
  }

  // try some more esoteric cases:
  {
    deallog << Utilities::split_string_list(" , , ", ' ').size() << std::endl;
  }

  {
    deallog << Utilities::split_string_list(" ", ' ').size() << std::endl;
    deallog << Utilities::split_string_list("   ", ' ').size() << std::endl;
  }

  Assert(Utilities::split_string_list(" ; ", ';').size() == 1,
         ExcInternalError());
  Assert(Utilities::split_string_list(" ; ", ';')[0] == "", ExcInternalError());
}



int
main()
{
  initlog();

  test();
}
