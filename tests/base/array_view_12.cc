// ---------------------------------------------------------------------
//
// Copyright (C) 2017 by the deal.II authors
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


// test the ArrayView constructor that converts from std::initializer_list
#include "../tests.h"

#include <deal.II/base/array_view.h>
#include <deal.II/grid/manifold_lib.h>


void test ()
{
  PolarManifold<2> polar_manifold;
  polar_manifold.get_new_point({Point<2>(1.0, 0.0), Point<2>(0.0, 1.0)},
  {0.5, 0.5});

  deallog << "OK" << std::endl;
}



int main()
{
  initlog();

  test ();
}
