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
                   FE_Nedelec<dim>(2),
                   2,
                   FESystem<dim>(FE_Q<dim>(1), dim),
                   1);

  const unsigned int n_comp = fe.n_components();

  test<2>(fe, F<2>(n_comp, 1), 1, false);
}
