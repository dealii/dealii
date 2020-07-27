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

// build an FE_Q1_Nonlocal finite element, where dofs are not distributed
// on vertices (with dpo[0]=1) but on cells (with dpo[dim+2]=vertices_per_cell)

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/tensor_product_polynomials.h>

#include <deal.II/fe/fe_q_base.h>

#include <iostream>

#include "../tests.h"

#include "../non_local/fe_q1_nonlocal.h"

template <int dim>
void
test()
{
  FE_Q1_Nonlocal<dim> fe;
  if (fe.n_dofs_per_cell() != GeometryInfo<dim>::vertices_per_cell)
    deallog << "Not OK" << std::endl;
  else
    deallog << "OK" << std::endl;
}

int
main()
{
  initlog();

  test<1>();
  test<2>();
  test<3>();
}
