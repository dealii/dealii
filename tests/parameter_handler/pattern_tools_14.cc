// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
// Copyright (C) 2010 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Detailed license information governing the source code and contributions
// can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
//
// -----------------------------------------------------------------------------

// Check Patterns::Tools::Convert for enum types
// Test the conversion of UpdateFlags to string and back

#include <deal.II/fe/fe_values.h>

#include <memory>

#include "../tests.h"

int
main()
{
  initlog();

  UpdateFlags flags = update_values;

  deallog << Patterns::Tools::Convert<UpdateFlags>::to_string(flags)
          << std::endl;

  flags = Patterns::Tools::Convert<UpdateFlags>::to_value(
    "update_values|update_gradients");

  deallog << Patterns::Tools::Convert<UpdateFlags>::to_string(flags)
          << std::endl;
}
