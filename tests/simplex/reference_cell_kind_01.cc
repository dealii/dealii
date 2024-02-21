// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2021 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// Test ReferenceCell::Kind::faces_for_given_vertex().


#include <deal.II/base/quadrature_lib.h>

#include <deal.II/grid/reference_cell.h>

#include "../tests.h"


template <int dim>
void
test(const ReferenceCell &reference_cell)
{
  const auto kind = ReferenceCell(reference_cell);

  for (const auto v : reference_cell.vertex_indices())
    {
      deallog << v << ": ";
      for (const auto i : kind.faces_for_given_vertex(v))
        deallog << i << ' ';
      deallog << std::endl;
    }
  deallog << std::endl;
}



int
main()
{
  initlog();

  test<2>(ReferenceCells::Line);
  test<2>(ReferenceCells::Triangle);
  test<2>(ReferenceCells::Quadrilateral);
  test<3>(ReferenceCells::Tetrahedron);
  test<3>(ReferenceCells::Pyramid);
  test<3>(ReferenceCells::Wedge);
  test<3>(ReferenceCells::Hexahedron);
}
