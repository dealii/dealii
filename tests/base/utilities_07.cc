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


// verify that Utilities::string_to_double actually catches errors

#include <deal.II/base/utilities.h>

#include <sstream>

#include "../tests.h"



void
verify(const std::string &s)
{
  bool exception_caught = false;
  try
    {
      Utilities::string_to_double(s);
    }
  catch (...)
    {
      exception_caught = true;
    }
  Assert(exception_caught == true, ExcMessage("Function is broken!"));

  deallog << "Done correctly: " << s << std::endl;
}



int
main()
{
  initlog();

  verify("abc");
  verify("1.23.4");
  verify("1 23 4");
  verify("123abc");
}
