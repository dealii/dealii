// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2012 - 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// Utilities::get_integer_at_position

#include <deal.II/base/utilities.h>

#include <sstream>

#include "../tests.h"



void
test()
{
  int number = 5;
  for (unsigned int i = 0; i < 7; ++i)
    {
      std::ostringstream s;
      s << "test test" << number << "test test";

      AssertThrow(Utilities::get_integer_at_position(s.str(), 9).first ==
                    number,
                  ExcInternalError());
      AssertThrow(Utilities::get_integer_at_position(s.str(), 9).second ==
                    i + 1,
                  ExcInternalError());

      deallog << i << ' '
              << Utilities::get_integer_at_position(s.str(), 9).first
              << std::endl;

      number = number * 10 + i;
    }
}



int
main()
{
  initlog();

  test();
}
