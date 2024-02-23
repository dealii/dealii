// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2022 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

// Check that has hanging nodes works correctly for anysotropic refinement

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>

#include "../tests.h"

int
main()
{
  initlog();
  Triangulation<2> tria;
  GridGenerator::hyper_cube(tria);
  tria.begin_active()->set_refine_flag(RefinementCase<2>::cut_x);
  tria.execute_coarsening_and_refinement();
  tria.begin_active()->set_refine_flag(RefinementCase<2>::cut_x);
  tria.execute_coarsening_and_refinement();

  if (tria.has_hanging_nodes() == false)
    {
      deallog << "OK" << std::endl;
    }
  else
    {
      deallog << "NOT OK" << std::endl;
    }
}
