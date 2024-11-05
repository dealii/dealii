// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// Putting a cubic FE_SimplexP into an FESystem triggered an error
// because we didn't yet know how to deal with hanging node
// constraints on simplices in FESystem.

#include <deal.II/fe/fe_simplex_p.h>
#include <deal.II/fe/fe_system.h>

#include "../tests.h"



int
main()
{
  initlog();

  FESystem<2> fe_system_2(FE_SimplexP<2>(3), 1);

  deallog << "Success!" << std::endl;
}
