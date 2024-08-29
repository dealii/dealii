// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2020 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

// Check mesh smoothing for a simplex mesh

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_simplex_p.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>

#include "../tests.h"


template <int dim>
void
test()
{
  Triangulation<dim> tria(Triangulation<dim>::maximum_smoothing);
  Triangulation<dim> hex_tria;
  GridGenerator::hyper_cube(hex_tria);
  hex_tria.refine_global(3);
  GridGenerator::convert_hypercube_to_simplex_mesh(hex_tria, tria);
  tria.refine_global();
  deallog << "cells: " << tria.n_active_cells() << std::endl;
  deallog << "OK" << std::endl;
}

int
main()
{
  initlog();

  test<2>();
}
