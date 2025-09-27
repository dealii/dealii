// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2007 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


#include "../tests.h"

#include "fe_prolongation_common.h"



int
main()
{
  initlog();

  CHECK_ALL(NedelecSZ, 0, 2);
  CHECK_ALL(NedelecSZ, 0, 3);
  CHECK_ALL(NedelecSZ, 1, 2);
  CHECK_ALL(NedelecSZ, 1, 3);
}
