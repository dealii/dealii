// ---------------------------------------------------------------------
//
// Copyright (C) 2019 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------


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
