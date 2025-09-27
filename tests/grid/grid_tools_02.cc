// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2001 - 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// check GridTools::diameter for codim-1 meshes


#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria.h>

#include "../tests.h"

template <int dim>
void
test1()
{
  // test 1: hypercube
  if (true)
    {
      Triangulation<dim, dim + 1> tria;
      GridGenerator::hyper_cube(tria);

      for (unsigned int i = 0; i < 2; ++i)
        {
          tria.refine_global(2);
          deallog << dim << "d, "
                  << "hypercube diameter, " << i * 2
                  << " refinements: " << GridTools::diameter(tria) << std::endl;
        }
    }
}



int
main()
{
  initlog();
  deallog << std::setprecision(4);

  test1<1>();
  test1<2>();

  return 0;
}
