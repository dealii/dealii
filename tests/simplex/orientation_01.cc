// ---------------------------------------------------------------------
//
// Copyright (C) 2021 by the deal.II authors
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



// Test ReferenceCell::permute_according_orientation and compute_orientation
// by reverting the permutation due to a given orientation.

#include <deal.II/grid/reference_cell.h>

#include "../tests.h"

template <unsigned int n_points>
void
test(const ReferenceCell type, const unsigned int n_orientations)
{
  for (unsigned int o = 0; o < n_orientations; ++o)
    {
      std::array<unsigned int, n_points> origin;

      for (unsigned int i = 0; i < n_points; ++i)
        origin[i] = i;

      const auto permuted = type.permute_according_orientation(origin, o);
      const unsigned int origin_back =
        type.compute_orientation(permuted, origin);

      AssertDimension(o, origin_back);
    }
}

int
main()
{
  initlog();

  test<2>(ReferenceCells::Line, 2);
  test<3>(ReferenceCells::Triangle, 3);
  test<4>(ReferenceCells::Quadrilateral, 4);

  deallog << "OK!" << std::endl;
}
