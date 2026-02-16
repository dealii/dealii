// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
// Copyright (C) 2019 - 2022 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Detailed license information governing the source code and contributions
// can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
//
// -----------------------------------------------------------------------------

// Output standard parameters for ParsedConvergenceTable class

#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/parsed_convergence_table.h>

#include "../tests.h"

int
main()
{
  initlog();

  ParsedConvergenceTable table;

  ParameterHandler prm;
  table.add_parameters(prm);

  prm.log_parameters(deallog);
}
