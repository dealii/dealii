// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
// Copyright (C) 2003 - 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Detailed license information governing the source code and contributions
// can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
//
// -----------------------------------------------------------------------------



// test ParameterHandler::Double description of limits

#include <deal.II/base/parameter_handler.h>

#include "../tests.h"


int
main()
{
  initlog();

  ParameterHandler prm;
  prm.declare_entry("a", "1.2", Patterns::Double(), "no limit");
  prm.declare_entry("b", "1.2", Patterns::Double(-2.13), "lower limit");
  prm.declare_entry("c",
                    "1.2",
                    Patterns::Double(Patterns::Double::min_double_value, 42.0),
                    "upper limit");
  prm.declare_entry("d", "1.2", Patterns::Double(0.2, 42.0), "both limits");
  prm.declare_entry("e", "1.2", Patterns::Double(1.0, -1.0), "no limits");

  prm.print_parameters(deallog.get_file_stream(), ParameterHandler::LaTeX);

  return 0;
}
