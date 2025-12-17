// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// Test ReferenceCell::opposite_face_index()
#include <deal.II/grid/reference_cell.h>

#include "../tests.h"


void
test(const ReferenceCell reference_cell)
{
  for (const auto &face_no : reference_cell.face_indices())
    {
      const unsigned int face_guess =
        reference_cell.opposite_face_index(face_no);
      deallog << "Face number: " << face_no
              << " opposite face guess: " << face_guess << std::endl;
    }
}

int
main()
{
  initlog();

  deallog.push("Quadrilateral");
  test(ReferenceCells::Quadrilateral);
  deallog.pop();

  deallog.push("Triangle");
  test(ReferenceCells::Triangle);
  deallog.pop();

  deallog.push("Hexahedron");
  test(ReferenceCells::Hexahedron);
  deallog.pop();

  deallog.push("Tetrahedron");
  test(ReferenceCells::Tetrahedron);
  deallog.pop();

  deallog.push("Pyramid");
  test(ReferenceCells::Pyramid);
  deallog.pop();

  deallog.push("Wedge");
  test(ReferenceCells::Wedge);
  deallog.pop();
}
