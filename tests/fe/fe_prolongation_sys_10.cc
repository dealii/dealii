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

  CHECK_SYS3(FE_Nedelec<2>(1),
             1,
             FESystem<2>(FE_DGQ<2>(3), 3),
             1,
             FESystem<2>(FE_Q<2>(2), 3, FE_Nedelec<2>(1), 2),
             2,
             2);
}
