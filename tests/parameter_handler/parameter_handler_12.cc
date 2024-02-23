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
