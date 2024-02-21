// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2006 - 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


/*
 * Project the function [1,1] onto a deformed grid and see whether the
 * FESystem elements can represent it exactly. This shouldn't be a surprise,
 * but it is nice to compare with the RT and ABF elements
 */



char logname[] = "output";
#include "deformed_projection.h"


void
test()
{
  FESystem<2> fe(FE_Q_Bubbles<2>(QIterated<1>(QTrapezoid<1>(), 3)), 2);
  const std::array<unsigned int, 3> min_convergence_steps = {{15, 15, 15}};
  check(fe, min_convergence_steps);
}
