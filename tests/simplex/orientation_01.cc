// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
// Copyright (C) 2021 - 2026 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Detailed license information governing the source code and contributions
// can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
//
// -----------------------------------------------------------------------------



// Test ReferenceCell::permute_according_orientation and compute_orientation
// by reverting the permutation due to a given orientation.

#include <deal.II/grid/reference_cell.h>

#include "../tests.h"

template <unsigned int n_points, int dim>
void
test(const ReferenceCell<dim> type, const unsigned int n_orientations)
{
  for (unsigned int o = 0; o < n_orientations; ++o)
    {
      std::array<unsigned int, n_points> origin;

      for (unsigned int i = 0; i < n_points; ++i)
        origin[i] = i;

      ArrayView<const unsigned int> origin_view(origin.cbegin(), n_points);
      const auto                    permuted =
        type.permute_by_combined_orientation(origin_view, o);
      const unsigned int origin_back = type.get_combined_orientation(
        make_array_view(permuted.cbegin(), permuted.cend()), origin_view);

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
