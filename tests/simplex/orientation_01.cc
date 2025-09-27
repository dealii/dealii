// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2021 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



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
