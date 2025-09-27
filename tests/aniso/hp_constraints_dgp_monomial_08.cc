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



// check that computation of hp-constraints works for DGP_Monomial elements
// correctly

char logname[] = "output";


#include "../hp/hp_constraints_common.h"


template <int dim>
void
test()
{
  hp::FECollection<dim> fe;
  for (unsigned int i = 0; i < 4; ++i)
    fe.push_back(FE_DGPMonomial<dim>(i));

  test_with_hanging_nodes_random_aniso(fe);
}
