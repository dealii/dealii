// ---------------------------------------------------------------------
//
// Copyright (C) 2017 - 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE.md at
// the top level directory of deal.II.
//
// ---------------------------------------------------------------------

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
