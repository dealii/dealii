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


// test for BoundingBox<unsigned int spacedim> which tests the functions
// volume and merge_with of the class

#include <deal.II/base/bounding_box.h>
#include <deal.II/base/point.h>

#include "../tests.h"

template <int spacedim>
void
test_bounding_box()
{
  BoundingBox<spacedim> a;
  deallog << "Empty BoundingBox: " << std::endl;
  deallog << a.volume() << std::endl;


  std::pair<Point<spacedim>, Point<spacedim>> unit;
  for (int i = 0; i < spacedim; i++)
    {
      unit.first[i]  = 0.0;
      unit.second[i] = 1.0;
    }

  BoundingBox<spacedim> b(unit);
  deallog << "Bounding box boundaries: " << std::endl;
  deallog << b.get_boundary_points().first << std::endl;
  deallog << b.get_boundary_points().second << std::endl;
  deallog << " Boundary Box volume: " << b.volume() << std::endl;
  b.merge_with(a);
  deallog << "Merging it with 0: volume " << b.volume() << std::endl;
  deallog << "Boundary points:" << std::endl;
  deallog << b.get_boundary_points().first << std::endl;
  deallog << b.get_boundary_points().second << std::endl;

  std::pair<Point<spacedim>, Point<spacedim>> boundaries;
  for (int i = 0; i < spacedim; i++)
    {
      unit.first[i]  = 1.0;
      unit.second[i] = 2.0 + i;
    }

  BoundingBox<spacedim> c(unit);
  deallog << "Bounding box boundaries: " << std::endl;
  deallog << c.get_boundary_points().first << std::endl;
  deallog << c.get_boundary_points().second << std::endl;
  deallog << c.volume() << std::endl;
  c.merge_with(b);
  deallog << "Merging it with previous: volume " << c.volume() << std::endl;
  deallog << "Boundary points:" << std::endl;
  deallog << c.get_boundary_points().first << std::endl;
  deallog << c.get_boundary_points().second << std::endl;
  deallog << "End test for dimension " << spacedim << std::endl;
}

int
main()
{
  initlog();

  deallog << "Test: Bounding Box class merge and volume functions" << std::endl;
  deallog << std::endl << "Test for dimension 1" << std::endl;
  test_bounding_box<1>();

  deallog << std::endl << "Test for dimension 2" << std::endl;
  test_bounding_box<2>();

  deallog << std::endl << "Test for dimension 3" << std::endl;
  test_bounding_box<3>();
}
