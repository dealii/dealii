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

// build an FE_Q1_Nonlocal finite element, distribute dofs on a simple
// Triangulation, and interpolate a scalar function.

#include <deal.II/base/function_lib.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>

#include <deal.II/numerics/vector_tools.h>

#include <iostream>

#include "../tests.h"

#include "../non_local/fe_q1_nonlocal.h"

template <int dim>
void
test()
{
  Triangulation<dim> tria;
  GridGenerator::hyper_cube(tria);
  FE_Q1_Nonlocal<dim> fe(tria);
  tria.refine_global(1);
  DoFHandler<dim> dh(tria);
  dh.distribute_dofs(fe);

  Tensor<1, dim> ones;
  for (unsigned int d = 0; d < dim; ++d)
    ones[d] = 1;

  Vector<double> solution(dh.n_dofs());

  VectorTools::interpolate(dh, Functions::Monomial<dim>(ones), solution);

  deallog << solution << std::endl;
}

int
main()
{
  initlog();

  test<1>();
  test<2>();
  test<3>();
}
