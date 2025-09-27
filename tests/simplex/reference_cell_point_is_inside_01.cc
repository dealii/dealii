// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2022 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// check ReferenceCell::contains_point()

#include <deal.II/grid/reference_cell.h>

#include "../tests.h"


template <int dim>
void
test(const ReferenceCell &reference_cell)
{
  deallog << "Reference cell is: " << reference_cell.to_string() << std::endl;

  deallog << std::boolalpha;

  // Generate equispaced points and check whether they are inside the box
  for (double x = -0.1; x <= 1.1; x += 0.1)
    for (double y = -0.1; y <= (dim > 1 ? 1.1 : -0.1); y += 0.1)
      for (double z = -0.1; z <= (dim > 2 ? 1.1 : -0.1); z += 0.1)
        {
          Point<dim> p =
            (dim == 1 ? Point<dim>(x) :
                        (dim == 2 ? Point<dim>(x, y) : Point<dim>(x, y, z)));

          deallog << p << ' ' << reference_cell.contains_point(p) << std::endl;
        }

  // Make sure that all vertices are inside:
  for (unsigned int v = 0; v < reference_cell.n_vertices(); ++v)
    Assert(reference_cell.contains_point(reference_cell.vertex<dim>(v)),
           ExcInternalError());

  // Make sure that all vertices are outside with a negative tolerance:
  for (unsigned int v = 0; v < reference_cell.n_vertices(); ++v)
    Assert(reference_cell.contains_point(reference_cell.vertex<dim>(v),
                                         -1e-12) == false,
           ExcInternalError());
}


int
main()
{
  initlog();

  test<1>(ReferenceCells::Line);
  test<2>(ReferenceCells::Quadrilateral);
  test<2>(ReferenceCells::Triangle);
  test<3>(ReferenceCells::Hexahedron);
  test<3>(ReferenceCells::Tetrahedron);

  // TODO: wedge and pyramid are not currently implemented

  return 0;
}
