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

// Call VectorTools::interpolate for a function in FE_FaceQ. The purpose is
// to check that all functions needed for interpolate to work exist.

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_face.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>

#include <deal.II/numerics/vector_tools.h>

#include "../tests.h"


template <int dim>
void
test()
{
  Triangulation<dim> triangulation;
  FE_FaceQ<dim>      fe(2);
  DoFHandler<dim>    dof_handler(triangulation);

  GridGenerator::hyper_cube(triangulation, 0, 1);
  triangulation.refine_global(6);

  dof_handler.distribute_dofs(fe);
  Vector<double> solution(dof_handler.n_dofs());

  VectorTools::interpolate(dof_handler,
                           Functions::ZeroFunction<dim>(),
                           solution);
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
