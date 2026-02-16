// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
// Copyright (C) 2007 - 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Detailed license information governing the source code and contributions
// can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
//
// -----------------------------------------------------------------------------


#include "../tests.h"

#include "fe_prolongation_common.h"



int
main()
{
  initlog();
  deallog << std::setprecision(9);

  CHECK_ALL(DGQ, 0, 2);
  CHECK_ALL(DGQ, 1, 2);
  CHECK_ALL(DGQ, 2, 2);
  CHECK_ALL(DGQ, 3, 2);
  CHECK_ALL(DGQ, 4, 2);

  CHECK_ALL(DGQ, 0, 3);
  CHECK_ALL(DGQ, 1, 3);
  CHECK_ALL(DGQ, 2, 3);
}
