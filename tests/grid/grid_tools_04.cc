// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
// Copyright (C) 2001 - 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Detailed license information governing the source code and contributions
// can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
//
// -----------------------------------------------------------------------------


// check GridTools::volume for codim-one


#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/manifold_lib.h>
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
                  << "hypercube volume, " << i * 2
                  << " refinements: " << GridTools::volume(tria) << std::endl;
        };
    };
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
