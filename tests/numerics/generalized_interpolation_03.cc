// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2017 - 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

// Check projection property of VectorTools::interpolate for
// Hdiv and Hcurl conforming spaces on something nontrivial.

#include <deal.II/fe/fe_nedelec.h>
#include <deal.II/fe/fe_raviart_thomas.h>

#include "../tests.h"

#include "generalized_interpolation.h"

int
main()
{
  deallog.depth_console(3);

  test<2>(FE_RaviartThomas<2>(0), F<2>(2, 1), 1, false);
  test<2>(FE_RaviartThomas<2>(1), F<2>(2, 0), 2, false);
  test<2>(FE_RaviartThomas<2>(1), F<2>(2, 2), 2, false);

  test<3>(FE_RaviartThomas<3>(0), F<3>(3, 0), 1, false);
  test<3>(FE_RaviartThomas<3>(1), F<3>(3, 0), 2, false);
  test<3>(FE_RaviartThomas<3>(1), F<3>(3, 2), 2, false);
}
