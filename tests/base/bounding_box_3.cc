// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2017 - 2021 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// test for BoundingBox<unsigned int spacedim> which tests the function
// get_neighbor_type

#include <deal.II/base/bounding_box.h>
#include <deal.II/base/point.h>

#include "../tests.h"

template <int spacedim>
void
test_bounding_box()
{
  std::pair<Point<spacedim>, Point<spacedim>> unit;
  for (int i = 0; i < spacedim; ++i)
    {
      unit.first[i]  = 0.0;
      unit.second[i] = 1.0;
    }

  BoundingBox<spacedim> a(unit);

  deallog << "Bounding box boundaries A: " << std::endl;
  deallog << a.get_boundary_points().first << std::endl;
  deallog << a.get_boundary_points().second << std::endl;
  deallog << "Is neighbor of itself: " << (int)a.get_neighbor_type(a)
          << std::endl;

  // This is a simple neighbor
  std::pair<Point<spacedim>, Point<spacedim>> second;
  for (int i = 0; i < spacedim; ++i)
    {
      second.first[i]  = 1.0;
      second.second[i] = 2.0;
    }

  BoundingBox<spacedim> b(second);

  deallog << "Bounding box boundaries B: " << std::endl;
  deallog << b.get_boundary_points().first << std::endl;
  deallog << b.get_boundary_points().second << std::endl;

  deallog << "Is neighbor A with B: " << (int)a.get_neighbor_type(b)
          << std::endl;
  deallog << "Is neighbor B with A: " << (int)b.get_neighbor_type(a)
          << std::endl;
  deallog << "Is neighbor B with B: " << (int)b.get_neighbor_type(b)
          << std::endl;

  // Testing the two 3d only cases
  if (spacedim == 3)
    {
      std::pair<Point<spacedim>, Point<spacedim>> second_1;
      second_1.first[0]  = 0.0;
      second_1.second[0] = 1.0;
      for (int i = 1; i < spacedim; ++i)
        {
          second_1.first[i]  = 1.0;
          second_1.second[i] = 2.0;
        }
      BoundingBox<spacedim> b_1(second_1);

      deallog << "Bounding box boundaries B1: " << std::endl;
      deallog << b_1.get_boundary_points().first << std::endl;
      deallog << b_1.get_boundary_points().second << std::endl;

      deallog << "Is neighbor A with B1: " << (int)a.get_neighbor_type(b_1)
              << std::endl;
      deallog << "Is neighbor B1 with A: " << (int)b_1.get_neighbor_type(a)
              << std::endl;
      deallog << "Is neighbor B1 with B: " << (int)b_1.get_neighbor_type(b_1)
              << std::endl;

      // Case of partial edge which is still simple neighbor
      // Note: in dimension 1 everything is only either not neighbor or
      // mergeable neighbor
      second_1.first[0] = 0.2;

      BoundingBox<spacedim> b_2(second_1);

      deallog << "Bounding box boundaries B2: " << std::endl;
      deallog << b_2.get_boundary_points().first << std::endl;
      deallog << b_2.get_boundary_points().second << std::endl;

      deallog << "Is neighbor A with B2: " << (int)a.get_neighbor_type(b_2)
              << std::endl;
      deallog << "Is neighbor B2 with A: " << (int)b_2.get_neighbor_type(a)
              << std::endl;
      deallog << "Is neighbor B2 with B: " << (int)b_2.get_neighbor_type(b_2)
              << std::endl;

      // Case which is attached_neighbor

      second_1.first[0] = 0.0;
      second_1.first[2] = 0.8;

      BoundingBox<spacedim> b_3(second_1);

      deallog << "Bounding box boundaries B3: " << std::endl;
      deallog << b_3.get_boundary_points().first << std::endl;
      deallog << b_3.get_boundary_points().second << std::endl;

      deallog << "Is neighbor A with B3: " << (int)a.get_neighbor_type(b_3)
              << std::endl;
      deallog << "Is neighbor B3 with A: " << (int)b_3.get_neighbor_type(a)
              << std::endl;
      deallog << "Is neighbor B3 with B: " << (int)b_3.get_neighbor_type(b_3)
              << std::endl;
    }

  // C contains both A and B
  BoundingBox<spacedim> c(std::make_pair(unit.first, second.second));
  deallog << "Bounding box boundaries C: " << std::endl;
  deallog << c.get_boundary_points().first << std::endl;
  deallog << c.get_boundary_points().second << std::endl;

  deallog << "Is neighbor C with B: " << (int)c.get_neighbor_type(b)
          << std::endl;
  deallog << "Is neighbor B with C: " << (int)b.get_neighbor_type(c)
          << std::endl;
  deallog << "Is neighbor A with C: " << (int)a.get_neighbor_type(c)
          << std::endl;
  deallog << "Is neighbor C with A: " << (int)c.get_neighbor_type(a)
          << std::endl;
  deallog << "Is neighbor C with C: " << (int)c.get_neighbor_type(c)
          << std::endl;

  // D contains A, intersects with B
  unit.second *= 1.4;
  BoundingBox<spacedim> d(unit);

  deallog << "Bounding box boundaries D: " << std::endl;
  deallog << d.get_boundary_points().first << std::endl;
  deallog << d.get_boundary_points().second << std::endl;

  deallog << "Is neighbor D with A: " << (int)d.get_neighbor_type(a)
          << std::endl;
  deallog << "Vice-versa: " << (int)a.get_neighbor_type(d) << std::endl;
  deallog << "Is neighbor D with B: " << (int)d.get_neighbor_type(b)
          << std::endl;
  deallog << "Vice-versa: " << (int)b.get_neighbor_type(d) << std::endl;
  deallog << "Is neighbor D with C: " << (int)d.get_neighbor_type(c)
          << std::endl;
  deallog << "Vice-versa: " << (int)c.get_neighbor_type(d) << std::endl;
  deallog << "Is neighbor D with D: " << (int)d.get_neighbor_type(d)
          << std::endl;

  for (int i = 0; i < spacedim; ++i)
    {
      second.first[i]  = -10.0;
      second.second[i] = -8.0;
    }

  // E is separated from all others
  BoundingBox<spacedim> e(second);

  deallog << "Bounding box boundaries E: " << std::endl;
  deallog << e.get_boundary_points().first << std::endl;
  deallog << e.get_boundary_points().second << std::endl;

  deallog << "Is neighbor E with A: " << (int)e.get_neighbor_type(a)
          << std::endl;
  deallog << "Vice-versa: " << (int)a.get_neighbor_type(e) << std::endl;
  deallog << "Is neighbor E with B: " << (int)e.get_neighbor_type(b)
          << std::endl;
  deallog << "Vice-versa: " << (int)b.get_neighbor_type(e) << std::endl;
  deallog << "Is neighbor E with C: " << (int)e.get_neighbor_type(c)
          << std::endl;
  deallog << "Vice-versa: " << (int)c.get_neighbor_type(e) << std::endl;
  deallog << "Is neighbor E with D: " << (int)e.get_neighbor_type(d)
          << std::endl;
  deallog << "Vice-versa: " << (int)d.get_neighbor_type(e) << std::endl;
  deallog << "Is neighbor E with E: " << (int)e.get_neighbor_type(e)
          << std::endl;

  // A and F are mergeable because one next to the other
  for (int i = 0; i < spacedim; ++i)
    {
      second.first[i]  = 0.0;
      second.second[i] = 1.0;
    }
  second.first[0] += 1.0;
  second.second[0] += 1.0;

  BoundingBox<spacedim> f(second);

  deallog << "Bounding box boundaries F: " << std::endl;
  deallog << f.get_boundary_points().first << std::endl;
  deallog << f.get_boundary_points().second << std::endl;

  deallog << "Is neighbor F with A: " << (int)f.get_neighbor_type(a)
          << std::endl;
  deallog << "Vice-versa: " << (int)a.get_neighbor_type(f) << std::endl;
  deallog << "Is neighbor F with B: " << (int)f.get_neighbor_type(b)
          << std::endl;
  deallog << "Vice-versa: " << (int)b.get_neighbor_type(f) << std::endl;
  deallog << "Is neighbor F with C: " << (int)f.get_neighbor_type(c)
          << std::endl;
  deallog << "Vice-versa: " << (int)c.get_neighbor_type(f) << std::endl;
  deallog << "Is neighbor F with D: " << (int)f.get_neighbor_type(d)
          << std::endl;
  deallog << "Vice-versa: " << (int)d.get_neighbor_type(f) << std::endl;
  deallog << "Is neighbor F with E: " << (int)f.get_neighbor_type(e)
          << std::endl;
  deallog << "Vice-versa: " << (int)e.get_neighbor_type(f) << std::endl;
  deallog << "Is neighbor F with F: " << (int)f.get_neighbor_type(f)
          << std::endl;

  deallog << "End test for dimension " << spacedim << std::endl;
  deallog << std::endl;

  // G is mergeable with A and F and C..
  second.first[0] -= 0.2;
  second.second[0] -= 0.2;

  BoundingBox<spacedim> g(second);

  deallog << "Bounding box boundaries G: " << std::endl;
  deallog << g.get_boundary_points().first << std::endl;
  deallog << g.get_boundary_points().second << std::endl;

  deallog << "Is neighbor G with A: " << (int)g.get_neighbor_type(a)
          << std::endl;
  deallog << "Vice-versa: " << (int)a.get_neighbor_type(g) << std::endl;
  deallog << "Is neighbor G with B: " << (int)g.get_neighbor_type(b)
          << std::endl;
  deallog << "Vice-versa: " << (int)b.get_neighbor_type(g) << std::endl;
  deallog << "Is neighbor G with C: " << (int)g.get_neighbor_type(c)
          << std::endl;
  deallog << "Vice-versa: " << (int)c.get_neighbor_type(g) << std::endl;
  deallog << "Is neighbor G with D: " << (int)g.get_neighbor_type(d)
          << std::endl;
  deallog << "Vice-versa: " << (int)d.get_neighbor_type(g) << std::endl;
  deallog << "Is neighbor G with E: " << (int)g.get_neighbor_type(e)
          << std::endl;
  deallog << "Vice-versa: " << (int)e.get_neighbor_type(g) << std::endl;
  deallog << "Is neighbor G with F: " << (int)g.get_neighbor_type(f)
          << std::endl;
  deallog << "Vice-versa: " << (int)f.get_neighbor_type(g) << std::endl;
  deallog << "Is neighbor G with G: " << (int)g.get_neighbor_type(g)
          << std::endl;

  deallog << "End test for dimension " << spacedim << std::endl;
  deallog << std::endl;
}

int
main()
{
  initlog();

  deallog << "Test: Bounding Box class Is neighbor and volume functions"
          << std::endl;
  deallog << std::endl << "Test for dimension 1" << std::endl;
  test_bounding_box<1>();

  deallog << std::endl << "Test for dimension 2" << std::endl;
  test_bounding_box<2>();

  deallog << std::endl << "Test for dimension 3" << std::endl;
  test_bounding_box<3>();
}
