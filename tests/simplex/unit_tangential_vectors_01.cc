// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2020 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// Test ReferenceCell::face_tangent_vector() and ::face_normal_vector()
// for all reference-cell types.


#include <deal.II/base/quadrature_lib.h>

#include <deal.II/grid/reference_cell.h>

#include "../tests.h"


template <int dim>
void
test(const ReferenceCell &reference_cell)
{
  for (const auto face_no : reference_cell.face_indices())
    {
      deallog << reference_cell.template face_normal_vector<dim>(face_no)
              << std::endl;

      for (unsigned int i = 0; i < dim - 1; ++i)
        deallog << reference_cell.template face_tangent_vector<dim>(face_no, i)
                << std::endl;
    }
  deallog << std::endl;
}



int
main()
{
  initlog();

  test<2>(ReferenceCells::Triangle);
  test<2>(ReferenceCells::Quadrilateral);
  test<3>(ReferenceCells::Tetrahedron);
  test<3>(ReferenceCells::Pyramid);
  test<3>(ReferenceCells::Wedge);
  test<3>(ReferenceCells::Hexahedron);
}
