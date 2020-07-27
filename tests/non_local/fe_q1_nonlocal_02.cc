// ---------------------------------------------------------------------
//
// Copyright (C) 1998 - 2018 by the deal.II authors
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

// build an FE_Q1_Nonlocal finite element, and check identities.

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>

#include <iostream>

#include "../tests.h"

#include "../non_local/fe_q1_nonlocal.h"

template <int dim>
void
test()
{
  deallog << "DIM: " << dim << std::endl << std::endl;
  FE_Q1_Nonlocal<dim> fe;
  // Check that the basis functions are what we expect
  const auto &quad = fe.get_unit_support_points();

  for (unsigned int i = 0; i < fe.n_dofs_per_cell(); ++i)
    {
      deallog << "Shape[" << i << "]: ";
      for (unsigned int q = 0; q < quad.size(); ++q)
        {
          const auto p = quad[q];
          deallog << int(fe.shape_value(i, p)) << " ";
        }
      deallog << std::endl;
    }
}

int
main()
{
  initlog();

  test<1>();
  test<2>();
  test<3>();
}
