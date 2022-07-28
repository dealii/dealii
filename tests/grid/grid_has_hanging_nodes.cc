// ---------------------------------------------------------------------
//
// Copyright (C) 2022 by the deal.II authors
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
