// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
// Copyright (C) 2006 - 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Detailed license information governing the source code and contributions
// can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
//
// -----------------------------------------------------------------------------


/*
 * Project the function [1,1] onto a deformed grid and see whether the RT
 * elements can represent it exactly.
 */



char logname[] = "output";
#include "deformed_projection.h"


void
test()
{
  FE_RaviartThomas<2>               fe(0);
  const std::array<unsigned int, 3> min_convergence_steps = {{3, 3, 3}};
  check(fe, min_convergence_steps);
}
