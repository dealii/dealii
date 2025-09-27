// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2007 - 2020 by the deal.II authors
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

#include "fe_support_points_common.h"


int
main()
{
  initlog();

  CHECK_ALL(Q_DG0, 1, 2);
  CHECK_ALL(Q_DG0, 2, 2);
  CHECK_ALL(Q_DG0, 3, 2);

  CHECK_ALL(Q_DG0, 1, 3);
  CHECK_ALL(Q_DG0, 2, 3);
}
