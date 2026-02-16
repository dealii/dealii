// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
// Copyright (C) 2002 - 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Detailed license information governing the source code and contributions
// can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
//
// -----------------------------------------------------------------------------



// One can (accidentally) call ParameterHandler::get_int/double on parameters
// that are really strings. There used to be a bug in these functions in that
// they didn't throw an error when the string wasn't actually convertible to a
// number

#include <deal.II/base/parameter_handler.h>

#include "../tests.h"

void
check()
{
  ParameterHandler prm;
  prm.declare_entry("a", "this that and the other", Patterns::Anything(), "");
  try
    {
      prm.get_double("a");
    }
  catch (...)
    {
      deallog << "get_double() detected the mistake" << std::endl;
    }

  try
    {
      prm.get_integer("a");
    }
  catch (...)
    {
      deallog << "get_integer() detected the mistake" << std::endl;
    }

  try
    {
      prm.get_bool("a");
    }
  catch (...)
    {
      deallog << "get_bool() detected the mistake" << std::endl;
    }
}


int
main()
{
  initlog();

  check();

  return 0;
}
