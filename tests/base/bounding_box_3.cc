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


// test for BoundingBox<unsigned int spacedim> which tests the function
// is_neighbour

#include "../tests.h"

#include <deal.II/base/point.h>
#include <deal.II/base/bounding_box.h>

template <int spacedim>
void test_bounding_box()
{

  std::pair<Point<spacedim>,Point<spacedim>> unit;
  for (int i=0; i<spacedim; i++)
    {
      unit.first[i] = 0.0;
      unit.second[i] = 1.0;
    }

  BoundingBox<spacedim> a(unit);

  deallog << "Bounding box boundaries A: " << std::endl;
  deallog << a.get_boundary_points().first << std::endl;
  deallog << a.get_boundary_points().second << std::endl;
  deallog << "Is neighbour of itself: " << a.is_neighbour(a) << std::endl;

  std::pair<Point<spacedim>,Point<spacedim>> second;
  for (int i=0; i<spacedim; i++)
    {
      second.first[i] = 1.0;
      second.second[i] = 2.0;
    }

  BoundingBox<spacedim> b(second);

  deallog << "Bounding box boundaries B: " << std::endl;
  deallog << b.get_boundary_points().first << std::endl;
  deallog << b.get_boundary_points().second << std::endl;

  deallog << "Is neighbour A with B: " << a.is_neighbour(b) << std::endl;
  deallog << "Is neighbour B with A: " << b.is_neighbour(a) << std::endl;
  deallog << "Is neighbour B with B: " << b.is_neighbour(b) << std::endl;

  BoundingBox<spacedim> c(std::make_pair(unit.first,second.second));
  deallog << "Bounding box boundaries C: " << std::endl;
  deallog << c.get_boundary_points().first << std::endl;
  deallog << c.get_boundary_points().second << std::endl;

  deallog << "Is neighbour C with B: " << c.is_neighbour(b) << std::endl;
  deallog << "Is neighbour B with C: " << b.is_neighbour(c) << std::endl;
  deallog << "Is neighbour A with C: " << a.is_neighbour(c) << std::endl;
  deallog << "Is neighbour C with A: " << c.is_neighbour(a) << std::endl;
  deallog << "Is neighbour C with C: " << c.is_neighbour(c) << std::endl;


  unit.second *= 1.4;
  BoundingBox<spacedim> d(unit);

  deallog << "Bounding box boundaries D: " << std::endl;
  deallog << d.get_boundary_points().first << std::endl;
  deallog << d.get_boundary_points().second << std::endl;

  deallog << "Is neighbour D with A: " << d.is_neighbour(a) << std::endl;
  deallog << "Vice-versa: " << a.is_neighbour(d) << std::endl;
  deallog << "Is neighbour D with B: " << d.is_neighbour(b) << std::endl;
  deallog << "Vice-versa: " << b.is_neighbour(d) << std::endl;
  deallog << "Is neighbour D with C: " << d.is_neighbour(c) << std::endl;
  deallog << "Vice-versa: " << c.is_neighbour(d) << std::endl;
  deallog << "Is neighbour D with D: " << d.is_neighbour(d) << std::endl;

  for (int i=0; i<spacedim; i++)
    {
      second.first[i] = -10.0;
      second.second[i] = -8.0;
    }

  BoundingBox<spacedim> e(second);

  deallog << "Bounding box boundaries E: " << std::endl;
  deallog << e.get_boundary_points().first << std::endl;
  deallog << e.get_boundary_points().second << std::endl;

  deallog << "Is neighbour E with A: " << e.is_neighbour(a) << std::endl;
  deallog << "Vice-versa: " << a.is_neighbour(e) << std::endl;
  deallog << "Is neighbour E with B: " << e.is_neighbour(b) << std::endl;
  deallog << "Vice-versa: " << b.is_neighbour(e) << std::endl;
  deallog << "Is neighbour E with C: " << e.is_neighbour(c) << std::endl;
  deallog << "Vice-versa: " << c.is_neighbour(e) << std::endl;
  deallog << "Is neighbour E with D: " << e.is_neighbour(d) << std::endl;
  deallog << "Vice-versa: " << d.is_neighbour(e) << std::endl;
  deallog << "Is neighbour E with E: " << e.is_neighbour(e) << std::endl;

  deallog << "End test for dimension " << spacedim << std::endl;
}

int main()
{
  initlog();

  deallog << "Test: Bounding Box class Is neighbour and volume functions" << std::endl;
  deallog << std::endl << "Test for dimension 1" << std::endl;
  test_bounding_box<1>();

  deallog << std::endl << "Test for dimension 2" << std::endl;
  test_bounding_box<2>();

  deallog << std::endl << "Test for dimension 3" << std::endl;
  test_bounding_box<3>();
}
