// ---------------------------------------------------------------------
//
// Copyright (C) 2002 - 2017 by the deal.II authors
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
