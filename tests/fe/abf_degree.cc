// ---------------------------------------------------------------------
//
// Copyright (C) 2017 - 2020 by the deal.II authors
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



// Verify that FE_ABF(r) returns that its polynomial degree is r+2, not r+1.

#include <deal.II/fe/fe_abf.h>

#include "../tests.h"



int
main()
{
  initlog();

  deallog << "dim=2:" << std::endl;
  for (unsigned int degree = 0; degree < 2; ++degree)
    deallog << degree << ": " << FE_ABF<2>(degree).degree << std::endl;
}
