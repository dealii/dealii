// ---------------------------------------------------------------------
//
// Copyright (C) 2003 - 2022 by the deal.II authors
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
