// ---------------------------------------------------------------------
//
// Copyright (C) 2015 - 2017 by the deal.II Authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------

// Call VectorTools::interpolate for a function in FE_TraceQ. The purpose is
// to check that all functions needed for interpolate to work exist.

#include "../tests.h"

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/fe/fe_trace.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>
#include <deal.II/numerics/vector_tools.h>

using namespace dealii;

template <int dim>
void
test()
{
  Triangulation<dim> triangulation;
  FE_TraceQ<dim>     fe(2);
  DoFHandler<dim>    dof_handler(triangulation);

  GridGenerator::hyper_cube(triangulation, 0, 1);
  triangulation.refine_global(6);

  dof_handler.distribute_dofs(fe);
  Vector<double> solution(dof_handler.n_dofs());

  VectorTools::interpolate(dof_handler, ZeroFunction<dim>(), solution);
  deallog << "Success, dim = " << dim << std::endl;
}

int
main()
{
  initlog();

  test<2>();
  test<3>();
  return 0;
}
