// ---------------------------------------------------------------------
//
// Copyright (C) 2018 by the deal.II authors
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


// Test: Triangulation<dim, spacedim>::reverse_vertices

#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include "../tests.h"

void
test()
{
  // 1d test
  std::vector<Point<1>> p_1d   = {Point<1>(1.), Point<1>(-1.)};
  std::vector<Point<1>> ref_1d = {Point<1>(-1.), Point<1>(1.)};
  Triangulation<1, 1>::reorder_vertices(p_1d);
  for (int i = 0; i < 2; ++i)
    AssertThrow(p_1d[i].distance(ref_1d[i]) < 1.0e-15,
                ExcMessage("1d test fail"));
  deallog << "1D test okay" << std::endl;

  // 2d test
  // corner case test with points on axis
  std::vector<Point<2>> p1 = {Point<2>(-1., 0.),
                              Point<2>(1., 0.),
                              Point<2>(0., -1.),
                              Point<2>(0., 1.)};
  std::vector<Point<2>> r1 = {Point<2>(0., -1.),
                              Point<2>(1., 0.),
                              Point<2>(-1., 0.),
                              Point<2>(0., 1.)};

  // regular case test
  std::vector<Point<2>> p2 = {Point<2>(-1.5, 0.1),
                              Point<2>(1., 0.3),
                              Point<2>(-0.2, -1.),
                              Point<2>(0., 1.)};
  std::vector<Point<2>> r2 = {Point<2>(-0.2, -1.),
                              Point<2>(1., 0.3),
                              Point<2>(-1.5, 0.1),
                              Point<2>(0., 1.)};
  // reorder vertices
  Triangulation<2, 2>::reorder_vertices(p1);
  Triangulation<2, 2>::reorder_vertices(p2);

  for (int i = 0; i < 4; ++i)
    {
      AssertThrow(p1[i].distance(r1[i]) < 1.0e-15,
                  ExcMessage("2d corner case test fail"));
      AssertThrow(p2[i].distance(r2[i]) < 1.0e-15,
                  ExcMessage("2d regular test fail"));
    }
  deallog << "2D test okay" << std::endl;
}


int
main()
{
  init_log();

  test();

  return 0;
}
