// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2019 - 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// Check that scalar differentiation can be performed using the utility
// functions

#include <deal.II/differentiation/sd.h>

#include "../tests.h"

namespace SD = Differentiation::SD;


int
main()
{
  initlog();

  using SD_number_t = SD::Expression;

  const SD_number_t x     = SD::make_symbol("x");
  const SD_number_t f     = 2 * x * x;
  const SD_number_t df_dx = SD::differentiate(f, x);

  deallog << "x: " << x << std::endl;
  deallog << "f: " << f << std::endl;
  deallog << "df_dx: " << df_dx << std::endl;

  deallog << "OK" << std::endl;
}
