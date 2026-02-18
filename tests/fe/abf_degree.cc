// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
// Copyright (C) 2017 - 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Detailed license information governing the source code and contributions
// can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
//
// -----------------------------------------------------------------------------



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
