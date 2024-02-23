// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2023 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

// Test ReferenceCell::get_nodal_type_quadrature()

#include <deal.II/base/quadrature.h>

#include <deal.II/grid/reference_cell.h>

#include "../tests.h"


template <int dim>
void
test(const ReferenceCell reference_cell)
{
  const Quadrature<dim> &quad =
    reference_cell.template get_nodal_type_quadrature<dim>();

  for (unsigned int q = 0; q < quad.size(); ++q)
    {
      // There is nothing to print for dim = 0 in the point
      if (dim > 0)
        deallog << quad.point(q) << ' ';
      deallog << quad.weight(q) << ' ';
      deallog << std::endl;
    }
}

int
main()
{
  initlog();

  {
    deallog.push("0d");
    deallog.pop();
  }
  {
    deallog.push("1d");
    test<1>(ReferenceCells::Line);
    deallog.pop();
  }
  {
    deallog.push("2d-1");
    test<2>(ReferenceCells::Triangle);
    deallog.pop();
  }
  {
    deallog.push("2d-2");
    test<2>(ReferenceCells::Quadrilateral);
    deallog.pop();
  }
  {
    deallog.push("3d-1");
    test<3>(ReferenceCells::Tetrahedron);
    deallog.pop();
  }
  {
    deallog.push("3d-2");
    test<3>(ReferenceCells::Pyramid);
    deallog.pop();
  }
  {
    deallog.push("3d-2");
    test<3>(ReferenceCells::Wedge);
    deallog.pop();
  }
  {
    deallog.push("3d-3");
    test<3>(ReferenceCells::Hexahedron);
    deallog.pop();
  }
}
