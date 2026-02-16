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

#include "fe_restriction_common.h"



int
main()
{
  initlog();

  CHECK_ALL(DGP, 0, 2);
  CHECK_ALL(DGP, 1, 2);
  CHECK_ALL(DGP, 2, 2);
  CHECK_ALL(DGP, 3, 2);
  CHECK_ALL(DGP, 4, 2);

  CHECK_ALL(DGP, 0, 3);
  CHECK_ALL(DGP, 1, 3);
  CHECK_ALL(DGP, 2, 3);
}
