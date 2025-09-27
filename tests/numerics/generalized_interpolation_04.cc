// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2017 - 2021 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

// Check projection property of VectorTools::interpolate for a complex,
// staggered system of Hdiv / Hcurl / L2 conforming spaces.

#include <deal.II/fe/fe_nedelec.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_raviart_thomas.h>
#include <deal.II/fe/fe_system.h>

#include "../tests.h"

#include "generalized_interpolation.h"

int
main()
{
  initlog();

  static const int dim = 2;

  FESystem<dim> fe(FE_RaviartThomas<dim>(1),
                   1,
                   FE_Q<dim>(1),
                   1,
                   FE_Nedelec<dim>(1),
                   2,
                   FESystem<dim>(FE_Q<dim>(1), dim),
                   1);

  const unsigned int n_comp = fe.n_components();

  test<2>(fe, F<2>(n_comp, 1), 1, false);
}
