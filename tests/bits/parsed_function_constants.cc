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
  prm.set("Function constants", "f=2");

  prm.set("Function expression", "f*f+2*x");

  Functions::ParsedFunction<1> function(1);
  function.parse_parameters(prm);

  const double value = function.value(Point<1>(3.0));
  deallog << "Value at x=3: " << value << std::endl;

  prm.log_parameters(deallog);
}

int
main()
{
  initlog();

  test();
}
