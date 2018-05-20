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

// testing serialize function for class BoundingBox

#include "../tests.h"

#include <deal.II/base/bounding_box.h>
#include <deal.II/base/point.h>
#include <deal.II/base/utilities.h>

template <int spacedim>
void
test(const unsigned int& size)
{
  std::vector<BoundingBox<spacedim>> b_boxes(size);

  for(auto& b : b_boxes)
    {
      auto p1 = random_point<spacedim>();
      auto p2 = p1;
      for(unsigned int i = 0; i < spacedim; ++i)
        p2[i] += i;
      b = BoundingBox<spacedim>(std::make_pair(p1, p2));
    }

  auto buffer = Utilities::pack(b_boxes);

  auto unpacked = Utilities::unpack<std::vector<BoundingBox<spacedim>>>(buffer);

  unsigned int i  = 0;
  bool         ok = true;
  for(auto& b : b_boxes)
    {
      const auto& b_points = b.get_boundary_points();
      const auto& u_points = unpacked[i++].get_boundary_points();
      if(b_points.first.distance(u_points.first) > 1e-12)
        {
          deallog << "NOT OK: " << b_points.first << " != " << u_points.first
                  << std::endl;
          ok = false;
        }
      else if(b_points.second.distance(u_points.second) > 1e-12)
        {
          deallog << "NOT OK: " << b_points.second << " != " << u_points.second
                  << std::endl;
          ok = false;
        }
    }

  if(ok)
    deallog << "OK!" << std::endl;
}

int
main()
{
  initlog();

  test<1>(10);
  test<2>(20);
  test<3>(30);
}
