// ---------------------------------------------------------------------
//
// Copyright (C) 2002 - 2018 by the deal.II authors
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



// test that using ParameterHandler::set with a parameter that doesn't conform
// to the specs leads to an error

#include <deal.II/base/parameter_handler.h>

#include "../tests.h"

void
check()
{
  ParameterHandler prm;
  prm.declare_entry("test_1", "3", Patterns::Integer());

  try
    {
      prm.set("test_1", "3.1415");
    }
  catch (const ParameterHandler::ExcValueDoesNotMatchPattern &)
    {
      deallog << "OK" << std::endl;
    }
  deallog << "test_1=" << prm.get("test_1") << std::endl;
}


int
main()
{
  initlog();

  check();

  return 0;
}
