// ---------------------------------------------------------------------
//
// Copyright (C) 2021 - 2022 by the deal.II authors
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


// Test ReferenceCell::Kind::faces_for_given_vertex().


#include <deal.II/base/quadrature_lib.h>

#include <deal.II/grid/reference_cell.h>

#include "../tests.h"

using namespace dealii;

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
