// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2021 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#include <deal.II/base/quadrature_lib.h>

#include <deal.II/fe/fe_pyramid_p.h>
#include <deal.II/fe/fe_simplex_p.h>
#include <deal.II/fe/fe_simplex_p_bubbles.h>
#include <deal.II/fe/fe_tools.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/fe_wedge_p.h>
#include <deal.II/fe/mapping_fe.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/grid_tools_cache.h>

#include "../tests.h"

// Check that Triangulation::get_boundary_ids() works for
// simplices. This used to crash.


template <int dim>
void
test_in_unit_cube()
{
  Triangulation<dim> tria;
  GridGenerator::subdivided_hyper_cube_with_simplices(tria, 1);

  for (const auto b : tria.get_boundary_ids())
    deallog << b << ' ';
  deallog << std::endl;
}

int
main()
{
  initlog();

  test_in_unit_cube<2>();
  test_in_unit_cube<3>();
}
