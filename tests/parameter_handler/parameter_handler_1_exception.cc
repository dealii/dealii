// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2002 - 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// ParameterHandler::declare_entry throws an exception if the default
// value of an entry doesn't match the pattern; but it should still
// yield a properly declared entry

#include <deal.II/base/parameter_handler.h>

#include "../tests.h"

void
check(const char *p)
{
  ParameterHandler prm;
  try
    {
      prm.declare_entry("test_1",
                        "abc",
                        Patterns::List(Patterns::Integer(-1, 1), 2, 3));
    }
  catch (const ParameterHandler::ExcValueDoesNotMatchPattern &)
    {
      deallog << "Exception caught as expected." << std::endl;
    }

  std::ifstream in(p);
  prm.parse_input(in);

  deallog << "test_1=" << prm.get("test_1") << std::endl;
}


int
main()
{
  initlog();

  check(SOURCE_DIR "/prm/parameter_handler_1_exception.prm");

  return 0;
}
