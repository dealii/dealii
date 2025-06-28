// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2022 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// Test consistency of ReferenceCell::face_to_cell_vertices(),
// ReferenceCell::vertex(), and ReferenceCell:face_vertex_location().

#include <deal.II/grid/reference_cell.h>

#include <random>

#include "../tests.h"



template <int dim>
void
test(const ReferenceCell &reference_cell)
{
  deallog << "ReferenceCell: " << reference_cell.to_string() << std::endl;
  for (unsigned int f = 0; f < reference_cell.n_faces(); ++f)
    for (unsigned int v = 0;
         v < reference_cell.face_reference_cell(f).n_vertices();
         ++v)
      {
        deallog << "Face=" << f << ", vertex=" << v << ", location="
                << reference_cell.face_vertex_location<dim>(f, v) << std::endl;
        Assert(reference_cell.face_vertex_location<dim>(f, v) ==
                 reference_cell.vertex<dim>(
                   reference_cell.face_to_cell_vertices(
                     f, v, numbers::default_geometric_orientation)),
               ExcInternalError());
      }
}

int
main()
{
  initlog();

  {
    deallog.push("1D");
    test<1>(ReferenceCells::Line);
    deallog.pop();
  }

  {
    deallog.push("2D");
    test<2>(ReferenceCells::Quadrilateral);
    test<2>(ReferenceCells::Triangle);
    deallog.pop();
  }

  {
    deallog.push("3D");
    test<3>(ReferenceCells::Tetrahedron);
    test<3>(ReferenceCells::Pyramid);
    test<3>(ReferenceCells::Wedge);
    test<3>(ReferenceCells::Hexahedron);
    deallog.pop();
  }
}
