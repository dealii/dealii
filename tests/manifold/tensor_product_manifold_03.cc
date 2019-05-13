// ---------------------------------------------------------------------
//
// Copyright (C) 2016 - 2019 by the deal.II authors
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


// Test the correct handling of data members within the TensorProductManifold
// class. The test generates a general hypercube grid and, in an inner scope,
// attempts to assign a TensorProductManifold to it. The test is successful if
// it terminates without subscriptor errors or memory leaks.

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/grid/tensor_product_manifold.h>

#include "../tests.h"


bool
test()
{
  Triangulation<3, 3> tria;
  GridGenerator::hyper_cube(tria);
  {
    FunctionManifold<1, 1> F("x", "x");
    PolarManifold<2, 2>    G(Point<2>(0.5, 0.5));
    try
      {
        TensorProductManifold<3, 2, 2, 2, 1, 1, 1> manifold(G, F);
        tria.set_all_manifold_ids(0);
        tria.set_manifold(0, manifold);
      }
    catch (std::bad_cast &e)
      {
        (void)e;
        return false;
      }
  }
  tria.refine_global(1);
  return true;
}



int
main()
{
  initlog();

  if (test())
    {
      deallog << "OK" << std::endl;
    }
  else
    {
      deallog << "FAILED" << std::endl;
    }

  return 0;
}
