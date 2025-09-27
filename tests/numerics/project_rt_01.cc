// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2006 - 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// check that VectorTools::project works for RT elements correctly

#include "project_common.h"


template <int dim>
void
test()
{
  if (dim != 1)
    // this is interesting also in 3d, but is
    // exceedingly slow there. limit to the
    // case of RT(0) elements in 3d
    for (unsigned int p = 0; p < (dim == 2 ? 3 : 1); ++p)
      test_no_hanging_nodes(FE_RaviartThomas<dim>(p), p + 1, 1);
}
