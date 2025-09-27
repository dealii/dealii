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



// this is reduced from hp_constraints_q_system_x_01, where we test that we
// can deal with FESystem(FE_Q(p),FE_DGQ(q)) for different p,q. note that for
// fixed p but varying q, neither of the two elements will dominate the other

char logname[] = "output";


#include "hp_constraints_common.h"


template <int dim>
void
test()
{
  hp::FECollection<dim> fe;
  fe.push_back(FESystem<dim>(FE_Q<dim>(1), 1, FE_DGQ<dim>(0), 1));
  fe.push_back(FESystem<dim>(FE_Q<dim>(1), 1, FE_DGQ<dim>(1), 1));

  test_no_hanging_nodes(fe);
}
