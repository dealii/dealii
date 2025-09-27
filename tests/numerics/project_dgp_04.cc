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



// check that VectorTools::project works for DGP elements correctly

// note that for mapped DGP(k) spaces, P(k) is not in the range of the
// mapping. However, P(k/2) is. Therefore, we have a gap for
// convergence, which we have to specify in the last argument to the
// call below

#include "project_common.h"


template <int dim>
void
test()
{
  for (unsigned int p = 0; p < 6 - dim; ++p)
    test_with_2d_deformed_mesh(FE_DGP<dim>(p), p, (p + 1) / 2);
}
