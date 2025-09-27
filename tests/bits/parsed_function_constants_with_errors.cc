// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2007 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// This program tests the functionality of the function parser
// wrapper.

#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/parsed_function.h>
#include <deal.II/base/utilities.h>

#include <map>

#include "../tests.h"

void
test()
{
  // A parameter handler
  ParameterHandler prm;

  Functions::ParsedFunction<1>::declare_parameters(prm, 1);

  // Setting a constant to an expression should fail:
  prm.set("Function constants", "f=100.0*60*60*24*365.2425");

  prm.set("Function expression", "f*f+2*x");

  Functions::ParsedFunction<1> function(1);

  bool caught_error = false;
  try
    {
      function.parse_parameters(prm);
    }
  catch (...)
    {
      deallog << "Caught error about invalid set for function constant."
              << std::endl;
      caught_error = true;
      // Ignore the actual error.
    }

  Assert(caught_error == true, ExcInternalError());
}

int
main()
{
  initlog();

  test();
}
